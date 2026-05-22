import logging
import re
from typing import Dict, List, Optional

# InChIKey pattern: 14 uppercase letters, hyphen, 10 uppercase letters, hyphen, 1 uppercase letter
INCHIKEY_PATTERN = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

from senita_chem.rdkit_properties import compute_rdkit_properties
from senita_chem.pubchem_batch import batch_lookup_by_inchikeys
from senita_chem.synonym_cleaning import clean_synonyms_list, get_cas_nos_from_synonyms_list

logger = logging.getLogger(__name__)


def enrich_compounds(
    compounds: Optional[List[Dict]] = None,
    inchikeys: Optional[List[str]] = None,
    max_synonyms: int = 75,
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

    if inchikeys:
        compounds = [{"smiles": "", "name": "", "_inchikey_only": ik} for ik in inchikeys]

    # --- Step 1: RDKit pass ---
    rdkit_cache: Dict[str, Dict] = {}
    inchikey_to_inputs: Dict[str, List[Dict]] = {}
    # Collect multi-fragment compounds for batch cholla_chem resolution
    multi_fragment_names: List[str] = []
    multi_fragment_compounds: List[Dict] = []

    for compound in compounds:
        smiles = compound.get("smiles", "")
        name = compound.get("name", "")
        inchikey_only = compound.get("_inchikey_only")

        if inchikey_only:
            rdkit_cache[inchikey_only] = {"inchikey": inchikey_only, "is_multi_fragment": False}
            inchikey_to_inputs.setdefault(inchikey_only, []).append(compound)
            continue

        props = compute_rdkit_properties(smiles)
        if props is None:
            logger.warning(f"Invalid SMILES, skipping: {smiles[:60]}")
            continue

        inchikey = props["inchikey"]

        # Set enrichment source and input tracking
        props.setdefault("enrichment_source", "pubchem")
        props["input_smiles"] = smiles
        props["input_name"] = name

        # Collect multi-fragment compounds for batch resolution
        if props["is_multi_fragment"] and name:
            multi_fragment_names.append(name)
            multi_fragment_compounds.append({
                "original_inchikey": inchikey,
                "original_smiles": smiles,
                "original_name": name,
            })

        rdkit_cache[inchikey] = props
        inchikey_to_inputs.setdefault(inchikey, []).append(compound)

    # --- Step 1b: Batch cholla_chem resolution for multi-fragment compounds ---
    if multi_fragment_names:
        logger.info(f"Batch resolving {len(multi_fragment_names)} multi-fragment compounds with cholla_chem")
        try:
            from cholla_chem import resolve_compounds_to_smiles
            resolved = resolve_compounds_to_smiles(compounds_list=multi_fragment_names)

            for compound_info in multi_fragment_compounds:
                name = compound_info["original_name"]
                resolved_smiles = resolved.get(name)
                # Guard: skip if resolved value looks like an InChIKey instead of SMILES
                if resolved_smiles and not INCHIKEY_PATTERN.match(resolved_smiles):
                    resolved_props = compute_rdkit_properties(resolved_smiles)
                    if resolved_props and not resolved_props["is_multi_fragment"]:
                        # Update the cache with resolved SMILES
                        original_inchikey = compound_info["original_inchikey"]

                        resolved_props["input_smiles"] = compound_info["original_smiles"]
                        resolved_props["input_name"] = name
                        resolved_props["enrichment_source"] = "cholla_chem+pubchem"

                        # Remove old multi-fragment entry and add resolved entry
                        del rdkit_cache[original_inchikey]
                        rdkit_cache[resolved_props["inchikey"]] = resolved_props
        except Exception as e:
            logger.warning(f"cholla_chem batch resolution failed: {e}")

    # --- Step 2: PubChem batch lookup ---
    all_inchikeys = list(rdkit_cache.keys())

    pubchem_results = batch_lookup_by_inchikeys(all_inchikeys)

    # --- Step 3: Merge ---
    for inchikey, rdkit_props in rdkit_cache.items():
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
