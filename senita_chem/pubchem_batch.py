import logging
import time
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

logger = logging.getLogger(__name__)

PUBCHEM_PUG_URL = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"

DEFAULT_REST_PROPERTIES = [
    "IsomericSMILES",
    "CanonicalSMILES",
    "InChI",
    "InChIKey",
    "IUPACName",
    "MolecularFormula",
    "Title",
]


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=10),
    retry=retry_if_exception_type(requests.exceptions.RequestException),
)
def send_xml_to_pubchem(xml_data: str) -> Optional[str]:
    headers = {"Content-Type": "application/xml"}
    response = requests.post(
        PUBCHEM_PUG_URL, data=xml_data, headers=headers, timeout=30
    )
    response.raise_for_status()
    return response.text


def create_batch_cid_request_xml(
    identifiers: List[str],
    identifier_type: str,
) -> str:
    """
    Creates PUG XML to exchange a list of identifiers for CIDs.
    Supports identifier_type: 'name', 'smiles', 'inchikey'.
    """
    root = ET.Element("PCT-Data")
    input_data = ET.SubElement(ET.SubElement(root, "PCT-Data_input"), "PCT-InputData")
    query = ET.SubElement(
        ET.SubElement(ET.SubElement(input_data, "PCT-InputData_query"), "PCT-Query"),
        "PCT-Query_type",
    )
    id_exchange = ET.SubElement(
        ET.SubElement(query, "PCT-QueryType"), "PCT-QueryType_id-exchange"
    )
    id_exchange = ET.SubElement(id_exchange, "PCT-QueryIDExchange")

    input_elem = ET.SubElement(id_exchange, "PCT-QueryIDExchange_input")
    query_uids = ET.SubElement(input_elem, "PCT-QueryUids")

    id_type_lower = identifier_type.lower()
    if id_type_lower == "smiles":
        container = ET.SubElement(query_uids, "PCT-QueryUids_smiles")
        tag = "PCT-QueryUids_smiles_E"
    elif id_type_lower == "inchikey":
        # InChIKeys must be passed as synonyms - PUG has no native InChIKey input type
        container = ET.SubElement(query_uids, "PCT-QueryUids_synonyms")
        tag = "PCT-QueryUids_synonyms_E"
    else:
        container = ET.SubElement(query_uids, "PCT-QueryUids_synonyms")
        tag = "PCT-QueryUids_synonyms_E"

    for identifier in identifiers:
        ET.SubElement(container, tag).text = str(identifier)

    ET.SubElement(id_exchange, "PCT-QueryIDExchange_operation-type").set("value", "same")
    ET.SubElement(id_exchange, "PCT-QueryIDExchange_output-type").set("value", "cid")
    ET.SubElement(id_exchange, "PCT-QueryIDExchange_output-method").set("value", "file-pair")
    ET.SubElement(id_exchange, "PCT-QueryIDExchange_compression").set("value", "none")

    return ET.tostring(root, encoding="unicode")


def create_status_request_xml(req_id: str) -> str:
    root = ET.Element("PCT-Data")
    input_data = ET.SubElement(ET.SubElement(root, "PCT-Data_input"), "PCT-InputData")
    request = ET.SubElement(
        ET.SubElement(input_data, "PCT-InputData_request"), "PCT-Request"
    )
    ET.SubElement(request, "PCT-Request_reqid").text = req_id
    ET.SubElement(request, "PCT-Request_type").set("value", "status")
    return ET.tostring(root, encoding="unicode")


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def download_file_from_pubchem(url: str) -> Optional[str]:
    https_url = url.replace("ftp://", "https://")
    response = requests.get(https_url, timeout=60)
    response.raise_for_status()
    return response.text


def poll_request_status(
    req_id: str,
    check_interval: int = 3,
    timeout: int = 120,
) -> Optional[str]:
    start_time = time.time()
    while time.time() - start_time < timeout:
        status_xml = create_status_request_xml(req_id)
        try:
            status_response = send_xml_to_pubchem(status_xml)
        except Exception as e:
            logger.error(f"Error checking status: {e}")
            return None

        root = ET.fromstring(status_response)
        status_elem = root.find(".//PCT-Status")
        if status_elem is None:
            logger.error("Could not find status element in response")
            return None

        status = status_elem.attrib.get("value", "unknown")
        if status == "success":
            download_url_elem = root.find(".//PCT-Download-URL_url")
            return download_url_elem.text if download_url_elem is not None else None
        elif status == "error":
            error_elem = root.find(".//PCT-Status-Message_message")
            logger.error(f"PubChem error: {error_elem.text if error_elem is not None else 'unknown'}")
            return None
        else:
            logger.debug(f"Request {req_id} status: {status}, waiting...")
            time.sleep(check_interval)

    logger.error(f"Timeout after {timeout}s for request {req_id}")
    return None


def parse_cid_file(file_content: str) -> Dict[str, str]:
    identifier_to_cid = {}
    for line in file_content.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            identifier, cid = parts[0].strip(), parts[1].strip()
            if cid.isdigit():
                identifier_to_cid[identifier] = cid
    return identifier_to_cid


def fetch_synonyms_for_cids(cids: List[str]) -> Dict[str, List[str]]:
    """Fetches all synonyms for a list of CIDs (max 1000) via PubChem REST."""
    if not cids:
        return {}
    cid_list = ",".join(str(cid) for cid in cids[:1000])
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_list}/synonyms/JSON"
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        data = response.json()
        results = {}
        for info in data.get("InformationList", {}).get("Information", []):
            cid = str(info.get("CID"))
            results[cid] = info.get("Synonym", [])
        return results
    except Exception as e:
        logger.error(f"Error fetching synonyms: {e}")
        return {}


def fetch_properties_for_cids(cids: List[str]) -> Dict[str, Dict]:
    """Fetches properties for a list of CIDs (max 1000) via PubChem REST."""
    if not cids:
        return {}
    cid_list = ",".join(str(cid) for cid in cids[:1000])
    property_list = ",".join(DEFAULT_REST_PROPERTIES)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_list}/property/{property_list}/JSON"
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        data = response.json()
        results = {}
        for prop in data.get("PropertyTable", {}).get("Properties", []):
            cid = str(prop.get("CID"))
            results[cid] = {
                "pubchem_cid": cid,
                "iupac_name": prop.get("IUPACName", ""),
                "preferred_name": prop.get("Title", ""),
                "canonical_smiles": prop.get("CanonicalSMILES", ""),
                "isomeric_smiles": prop.get("IsomericSMILES", ""),
                "inchi": prop.get("InChI", ""),
                "inchikey": prop.get("InChIKey", ""),
                "formula": prop.get("MolecularFormula", ""),
            }
        return results
    except Exception as e:
        logger.error(f"Error fetching properties: {e}")
        return {}


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
                        "canonical_smiles": prop.get("CanonicalSMILES", ""),
                        "isomeric_smiles": prop.get("IsomericSMILES", ""),
                        "inchi": prop.get("InChI", ""),
                        "inchikey": inchikey,
                        "formula": prop.get("MolecularFormula", ""),
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
    chunk_size: int = 50,
    check_interval: int = 3,
    timeout: int = 120,
    rest_threshold: int = 999999,
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
