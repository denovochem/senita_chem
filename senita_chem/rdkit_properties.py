import os
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


def compute_rdkit_properties(smiles: str) -> Optional[Dict]:
    """
    Computes all RDKit physicochemical properties from a SMILES string.

    Returns None if the SMILES is invalid.
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


def compute_rdkit_properties_batch(
    smiles_list: List[str],
    max_workers: Optional[int] = None,
    chunksize: Optional[int] = None,
) -> List[Optional[Dict]]:
    """
    Compute RDKit physicochemical properties for a list of SMILES strings in parallel.

    Uses a process pool to bypass the GIL — RDKit's Python wrappers hold the GIL
    for most operations, so threads provide little real parallelism. Separate
    processes give true CPU-level parallelism at the cost of IPC serialization.

    Args:
        smiles_list (List[str]): SMILES strings to process. Order is preserved in the output.
        max_workers (Optional[int]): Maximum number of worker processes. Defaults to
            ``min(len(smiles_list), os.cpu_count())``.
        chunksize (Optional[int]): Number of SMILES per task chunk. Larger values
            reduce IPC overhead. Defaults to ``ceil(len(smiles_list) / (4 * max_workers))``.

    Returns:
        List[Optional[Dict]]: One result per input SMILES, in the same order.
            Each element is the dict from ``compute_rdkit_properties`` or ``None``
            if the SMILES was invalid.
    """
    if not smiles_list:
        return []

    if max_workers is None:
        max_workers = min(len(smiles_list), os.cpu_count() or 4)

    if chunksize is None:
        chunksize = max(1, ceil(len(smiles_list) / (4 * max_workers)))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        return list(
            executor.map(compute_rdkit_properties, smiles_list, chunksize=chunksize)
        )
