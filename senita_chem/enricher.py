import json
import logging
import re
from functools import lru_cache
from importlib import resources
from typing import Dict, List, Optional

from openclatura import name_smiles

from senita_chem.local_pubchem_batch import batch_lookup_by_inchikeys_sqlite
from senita_chem.pubchem_batch import batch_lookup_by_inchikeys
from senita_chem.rdkit_properties import compute_rdkit_properties
from senita_chem.synonym_cleaning import (
    clean_synonyms_list,
    get_cas_nos_from_synonyms_list,
)

# InChIKey pattern: 14 uppercase letters, hyphen, 10 uppercase letters, hyphen, 1 uppercase letter
INCHIKEY_PATTERN = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

logger = logging.getLogger(__name__)


@lru_cache(maxsize=65536)
def _cached_compute_rdkit_properties(smiles: str) -> Optional[Dict]:
    """Compute RDKit properties with module-level LRU caching."""
    return compute_rdkit_properties(smiles)


def enrich_compounds(
    compounds: Optional[List[Dict]] = None,
    inchikeys: Optional[List[str]] = None,
    max_synonyms: int = 75,
    pubchem_method: str = "api",
    db_path: Optional[str] = None,
    precomputed_rdkit: Optional[Dict[str, Dict]] = None,
) -> Dict[str, Dict]:
    """
    Enrich compounds with RDKit properties, PubChem data, and OpenClatura IUPAC names.

    Accepts either:
      - compounds: List[{smiles, name}]
      - inchikeys: List[str]

    For each compound, RDKit properties are computed and a PubChem batch lookup is
    performed. When PubChem returns data, synonyms, CAS numbers, and PubChem metadata
    are merged into the record. When PubChem has no data, OpenClatura's ``name_smiles``
    is called on the canonical SMILES to generate an IUPAC name, which is used as the
    ``iupac_name``, ``preferred_name``, and sole entry in ``synonyms``. Single-fragment
    compounds without PubChem data are marked ``enrichment_source='failed'``; multi-
    fragment compounds are marked ``enrichment_source='rdkit_only'``.

    Args:
        compounds (Optional[List[Dict]]): List of dicts with ``smiles`` and ``name`` keys.
            Entries may also contain an ``inchikey`` key and a ``precomputed_rdkit`` key
            (used when ``precomputed_rdkit`` is provided).
        inchikeys (Optional[List[str]]): List of InChIKey strings to look up directly.
        max_synonyms (int): Maximum number of synonyms to retain after cleaning.
        pubchem_method (str): Either ``'api'`` for PubChem REST API or ``'local_db'``
            for a local SQLite database.
        db_path (Optional[str]): Path to the local PubChem SQLite database. Required
            when ``pubchem_method='local_db'``.
        precomputed_rdkit (Optional[Dict[str, Dict]]): Mapping from InChIKey to
            pre-computed RDKit properties. When provided, the RDKit computation pass
            is skipped for any compound whose InChIKey is found in this dict.

    Returns:
        Dict[str, Dict]: Results keyed by InChIKey. Each value is a dict containing
        RDKit properties, PubChem data (if found), OpenClatura IUPAC name (if PubChem
        had no data), synonyms, CAS numbers, common names, and an ``enrichment_source``
        field indicating the provenance (``'pubchem'``, ``'rdkit_only'``, or ``'failed'``).

    Raises:
        ValueError: If neither ``compounds`` nor ``inchikeys`` is provided, if
            ``pubchem_method`` is invalid, or if ``db_path`` is missing when
            ``pubchem_method='local_db'``.
    """
    if compounds is None and inchikeys is None:
        raise ValueError("Provide either compounds or inchikeys.")

    with (
        resources.files("senita_chem.data")
        .joinpath("enrichment_dict.json")
        .open("r") as fh
    ):
        enrichment_dict: Dict[str, List[str]] = json.load(fh)

    results: Dict[str, Dict] = {}

    if inchikeys is not None:
        compounds = [
            {"smiles": "", "name": "", "_inchikey_only": ik} for ik in inchikeys
        ]

    if compounds is None:
        raise ValueError("compounds must not be None after normalization.")

    # --- Step 1: RDKit pass (or use pre-computed properties) ---
    rdkit_batch: Dict[str, Dict] = {}
    inchikey_to_inputs: Dict[str, List[Dict]] = {}

    for compound in compounds:
        smiles = compound.get("smiles", "")
        name = compound.get("name", "")
        inchikey_only = compound.get("_inchikey_only")

        if inchikey_only:
            rdkit_batch[inchikey_only] = {
                "inchikey": inchikey_only,
                "is_multi_fragment": False,
            }
            inchikey_to_inputs.setdefault(inchikey_only, []).append(compound)
            continue

        # Try pre-computed properties first (keyed by inchikey)
        if precomputed_rdkit:
            ik = compound.get("inchikey")
            if ik and ik in precomputed_rdkit:
                props = precomputed_rdkit[ik].copy()
                props.setdefault("enrichment_source", "pubchem")
                props["input_smiles"] = smiles
                props["input_name"] = name
                rdkit_batch[ik] = props
                inchikey_to_inputs.setdefault(ik, []).append(compound)
                continue

        # Fall back to computing RDKit properties
        props = _cached_compute_rdkit_properties(smiles)
        if props is None:
            logger.warning(f"Invalid SMILES, skipping: {smiles[:60]}")
            continue

        inchikey = props["inchikey"]

        # Set enrichment source and input tracking
        props.setdefault("enrichment_source", "pubchem")
        props["input_smiles"] = smiles
        props["input_name"] = name

        rdkit_batch[inchikey] = props
        inchikey_to_inputs.setdefault(inchikey, []).append(compound)

    # --- Step 2: PubChem batch lookup ---
    all_inchikeys = list(rdkit_batch.keys())

    if pubchem_method == "local_db":
        if db_path is None:
            raise ValueError("db_path is required when pubchem_method='local_db'.")
        pubchem_results = batch_lookup_by_inchikeys_sqlite(
            all_inchikeys, db_path=db_path
        )
    elif pubchem_method == "api":
        pubchem_results = batch_lookup_by_inchikeys(all_inchikeys)
    else:
        raise ValueError(f"Invalid pubchem_method: {pubchem_method}")

    # --- Step 3: Merge ---
    for inchikey, rdkit_props in rdkit_batch.items():
        record = rdkit_props.copy()
        pubchem = pubchem_results.get(inchikey, {})

        if pubchem:
            raw_synonyms = pubchem.pop("raw_synonyms", [])
            record.update(pubchem)
            record["synonyms"] = clean_synonyms_list(raw_synonyms, max_synonyms)
            record["cas"] = get_cas_nos_from_synonyms_list(raw_synonyms)
        else:
            try:
                generated_iupac_name = name_smiles(rdkit_props["canonical_smiles"])
            except Exception:
                generated_iupac_name = ""

            record["synonyms"] = [generated_iupac_name] if generated_iupac_name else []
            record["cas"] = []
            record["pubchem_cid"] = None
            record["iupac_name"] = generated_iupac_name
            record["preferred_name"] = generated_iupac_name
            if record.get("is_multi_fragment"):
                record["enrichment_source"] = "rdkit_only"
            else:
                record["enrichment_source"] = "failed"

        common_names = enrichment_dict.get(inchikey, [])
        preferred = record.get("preferred_name", "")
        if preferred:
            norm_preferred = preferred.lower().replace(" ", "")
            common_names = [
                cn
                for cn in common_names
                if cn.lower().replace(" ", "") != norm_preferred
            ]
        record["common_names"] = common_names

        results[inchikey] = record

    return results
