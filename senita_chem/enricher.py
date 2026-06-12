import logging
import re
from functools import lru_cache
from typing import Dict, List, Optional

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


@lru_cache(maxsize=16384)
def _cached_compute_rdkit_properties(smiles: str) -> Optional[Dict]:
    """Compute RDKit properties with module-level LRU caching."""
    return compute_rdkit_properties(smiles)


def enrich_compounds(
    compounds: Optional[List[Dict]] = None,
    inchikeys: Optional[List[str]] = None,
    max_synonyms: int = 150,
    pubchem_method: str = "local_db",
    db_path: Optional[str] = None,
) -> Dict[str, Dict]:
    """
    Main entry point. Accepts either:
      - compounds: List[{smiles, name}]
      - inchikeys: List[str]

    Returns dict keyed by InChIKey.
    """
    if compounds is None and inchikeys is None:
        raise ValueError("Provide either compounds or inchikeys.")

    results: Dict[str, Dict] = {}

    if inchikeys is not None:
        compounds = [
            {"smiles": "", "name": "", "_inchikey_only": ik} for ik in inchikeys
        ]

    if compounds is None:
        raise ValueError("compounds must not be None after normalization.")

    # --- Step 1: RDKit pass ---
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
            record["synonyms"] = []
            record["cas"] = []
            record["pubchem_cid"] = None
            record["iupac_name"] = ""
            record["preferred_name"] = ""
            if record.get("is_multi_fragment"):
                record["enrichment_source"] = "rdkit_only"
            else:
                record["enrichment_source"] = "failed"

        results[inchikey] = record

    return results
