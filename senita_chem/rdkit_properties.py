import os
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from math import ceil
from typing import Dict, List, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]

_taut_opts = rdMolStandardize.CleanupParameters()
_taut_opts.tautomerRemoveSp3Stereo = False  # type: ignore[assignment]
_taut_opts.tautomerRemoveBondStereo = False  # type: ignore[assignment]
TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator(_taut_opts)

_CACHE_MAXSIZE = 65536
_cache: "OrderedDict[str, Optional[Dict]]" = OrderedDict()


def _cache_put(smiles: str, props: Optional[Dict]) -> None:
    """
    Insert a SMILES → properties mapping into the shared cache, evicting the
    least recently used entry when the cache is full.

    Args:
        smiles (str): The SMILES string to cache.
        props (Optional[Dict]): The computed properties (or None if invalid).
    """
    _cache[smiles] = props
    if len(_cache) > _CACHE_MAXSIZE:
        _cache.popitem(last=False)


def cache_clear() -> None:
    """Clear the module-level LRU cache."""
    _cache.clear()


def _compute_rdkit_properties_impl(smiles: str) -> Optional[Dict]:
    """
    Compute all RDKit physicochemical properties from a SMILES string.

    This is the uncached implementation used by both ``compute_rdkit_properties``
    (single-call, cached) and ``compute_rdkit_properties_batch`` (parallel,
    cache-aware).

    Args:
        smiles (str): A SMILES string to process.

    Returns:
        Optional[Dict]: Property dict if the SMILES is valid, otherwise None.
    """
    if not smiles or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Stereocenters must be captured before canonicalization strips them
    num_stereocenters = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    num_undefined_stereocenters = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(
        mol
    )

    # Capture pre-tautomer canonical SMILES and InChIKey
    original_canonical_smiles = Chem.MolToSmiles(
        mol, isomericSmiles=True, canonical=True
    )
    original_inchi = MolToInchi(mol)
    original_inchikey = InchiToInchiKey(original_inchi) if original_inchi else None

    mol = TAUTOMER_ENUMERATOR.Canonicalize(mol)

    inchi = MolToInchi(mol)
    inchikey = InchiToInchiKey(inchi) if inchi else None
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    is_multi_fragment = "." in canonical_smiles

    return {
        "canonical_smiles": canonical_smiles,
        "original_canonical_smiles": original_canonical_smiles,
        "original_inchikey": original_inchikey,
        "tautomer_smiles": canonical_smiles,
        "tautomer_inchikey": inchikey,
        "inchi": inchi,
        "inchikey": inchikey,
        "is_multi_fragment": is_multi_fragment,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 4),  # type: ignore[attr-defined]
        "mw_exact": round(Descriptors.ExactMolWt(mol), 4),  # type: ignore[attr-defined]
        "logp": round(Descriptors.MolLogP(mol), 4),  # type: ignore[attr-defined]
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 4),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "num_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "num_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "num_heterocycles": rdMolDescriptors.CalcNumHeterocycles(mol),
        "frac_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 4),
        "num_stereocenters": num_stereocenters,
        "num_defined_stereocenters": num_stereocenters,
        "num_undefined_stereocenters": num_undefined_stereocenters,
        "formal_charge": Chem.GetFormalCharge(mol),
    }


def compute_rdkit_properties(smiles: str) -> Optional[Dict]:
    """
    Compute all RDKit physicochemical properties from a SMILES string, with
    module-level LRU caching.

    Args:
        smiles (str): A SMILES string to process.

    Returns:
        Optional[Dict]: Property dict if the SMILES is valid, otherwise None.
        Both valid results and None (invalid SMILES) are cached.
    """
    if smiles in _cache:
        _cache.move_to_end(smiles)
        return _cache[smiles]
    result = _compute_rdkit_properties_impl(smiles)
    _cache_put(smiles, result)
    return result


def compute_rdkit_properties_batch(
    smiles_list: List[str],
    max_workers: Optional[int] = None,
    chunksize: Optional[int] = None,
) -> List[Optional[Dict]]:
    """
    Compute RDKit physicochemical properties for a list of SMILES strings in
    parallel, with cache-aware dispatch.

    Uses a process pool to bypass the GIL — RDKit's Python wrappers hold the GIL
    for most operations, so threads provide little real parallelism. Separate
    processes give true CPU-level parallelism at the cost of IPC serialization.

    Cache hits are served instantly from the module-level LRU cache; only cache
    misses are dispatched to worker processes.  Results from workers are stored
    back in the cache so subsequent calls (single or batch) benefit from prior
    computations.

    Args:
        smiles_list (List[str]): SMILES strings to process. Order is preserved
            in the output.
        max_workers (Optional[int]): Maximum number of worker processes. Defaults
            to ``min(len(misses), os.cpu_count() or 4)``.
        chunksize (Optional[int]): Number of SMILES per task chunk. Larger values
            reduce IPC overhead. Defaults to
            ``ceil(len(misses) / (4 * max_workers))``.

    Returns:
        List[Optional[Dict]]: One result per input SMILES, in the same order.
            Each element is the dict from ``compute_rdkit_properties`` or ``None``
            if the SMILES was invalid.
    """
    if not smiles_list:
        return []

    results: List[Optional[Dict]] = [None] * len(smiles_list)
    misses: List[str] = []
    miss_indices: List[int] = []

    for i, smiles in enumerate(smiles_list):
        if smiles in _cache:
            _cache.move_to_end(smiles)
            results[i] = _cache[smiles]
        else:
            misses.append(smiles)
            miss_indices.append(i)

    if not misses:
        return results

    if max_workers is None:
        max_workers = min(len(misses), os.cpu_count() or 4)

    if chunksize is None:
        chunksize = max(1, ceil(len(misses) / (4 * max_workers)))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        computed = list(
            executor.map(_compute_rdkit_properties_impl, misses, chunksize=chunksize)
        )

    for idx, smiles, props in zip(miss_indices, misses, computed):
        results[idx] = props
        _cache_put(smiles, props)

    return results
