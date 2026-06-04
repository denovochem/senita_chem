import logging
from typing import Dict, List

import requests

logger = logging.getLogger(__name__)

DEFAULT_REST_PROPERTIES = [
    "IsomericSMILES",
    "CanonicalSMILES",
    "InChI",
    "InChIKey",
    "IUPACName",
    "MolecularFormula",
    "Title",
]


def direct_rest_lookup_by_inchikeys(inchikeys: List[str]) -> Dict[str, Dict]:
    """
    Direct REST lookup for any size batch, chunked at 50 (REST endpoint limit).
    Always uses REST - bypasses PUG entirely for reliability.
    """
    if not inchikeys:
        return {}

    all_results: Dict[str, Dict] = {}
    property_list = ",".join(DEFAULT_REST_PROPERTIES)

    # Chunk at 50 (REST endpoint limit)
    for i in range(0, len(inchikeys), 50):
        chunk = inchikeys[i : i + 50]
        chunk_num = (i // 50) + 1
        logger.info(f"REST chunk {chunk_num}: {len(chunk)} InChIKeys")

        inchikey_list = ",".join(chunk)
        props_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey_list}/property/{property_list}/JSON"
        synonyms_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey_list}/synonyms/JSON"

        try:
            # Fetch properties
            props_response = requests.get(props_url, timeout=60)
            props_response.raise_for_status()
            props_data = props_response.json()

            # Fetch synonyms
            synonyms_response = requests.get(synonyms_url, timeout=60)
            synonyms_response.raise_for_status()
            synonyms_data = synonyms_response.json()

            # Process properties for this chunk
            chunk_results: Dict[str, Dict] = {}
            for prop in props_data.get("PropertyTable", {}).get("Properties", []):
                inchikey = prop.get("InChIKey", "")
                if inchikey:
                    chunk_results[inchikey] = {
                        "pubchem_cid": str(prop.get("CID", "")),
                        "iupac_name": prop.get("IUPACName", ""),
                        "preferred_name": prop.get("Title", ""),
                        "canonical_smiles_pubchem": prop.get("CanonicalSMILES", ""),
                        "isomeric_smiles_pubchem": prop.get("IsomericSMILES", ""),
                        "inchi_pubchem": prop.get("InChI", ""),
                        "inchikey_pubchem": inchikey,
                        "formula_pubchem": prop.get("MolecularFormula", ""),
                    }

            # Build CID -> InChIKey mapping from chunk results
            cid_to_inchikey = {
                chunk_results[ik]["pubchem_cid"]: ik
                for ik in chunk_results
                if chunk_results[ik].get("pubchem_cid")
            }

            # Add synonyms (synonyms endpoint returns CID, not InChIKey)
            for info in synonyms_data.get("InformationList", {}).get("Information", []):
                cid = str(info.get("CID", ""))
                if cid in cid_to_inchikey:
                    inchikey = cid_to_inchikey[cid]
                    chunk_results[inchikey]["raw_synonyms"] = info.get("Synonym", [])

            # Ensure all results have raw_synonyms
            for inchikey in chunk_results:
                if "raw_synonyms" not in chunk_results[inchikey]:
                    chunk_results[inchikey]["raw_synonyms"] = []

            # Merge chunk results into all_results
            all_results.update(chunk_results)
            logger.info(f"REST chunk {chunk_num} done: {len(chunk_results)} results")

        except Exception as e:
            logger.error(f"Error in REST chunk {chunk_num}: {e}")
            continue

    return all_results


def batch_lookup_by_inchikeys(
    inchikeys: List[str],
) -> Dict[str, Dict]:
    """
    Main batch lookup: InChIKey list → PubChem properties + raw synonyms.

    Always uses direct REST (chunked at 50) - PUG path is disabled for reliability.

    Returns dict keyed by InChIKey:
    {
        inchikey: {
            pubchem_cid, iupac_name, preferred_name, canonical_smiles,
            isomeric_smiles, inchi, inchikey, formula,
            raw_synonyms: [...]
        }
    }
    """
    unique_keys = list(dict.fromkeys(inchikeys))
    logger.info(f"Using direct REST for {len(unique_keys)} InChIKeys (PUG disabled)")
    return direct_rest_lookup_by_inchikeys(unique_keys)
