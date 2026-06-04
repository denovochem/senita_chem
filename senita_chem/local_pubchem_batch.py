import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Union


def batch_lookup_by_inchikeys_sqlite(
    inchikeys: List[str], db_path: Union[str, Path]
) -> Dict[str, Dict[str, Any]]:
    """
    Look up InChIKeys in a local PubChem SQLite database.

    Queries the compounds and synonyms tables to return CID, IUPAC name,
    preferred name (label), and all synonyms for each matched InChIKey.

    Args:
        inchikeys: List of InChIKeys to look up
        db_path: Path to the SQLite database file

    Returns:
        Dictionary mapping each InChIKey to a dict with keys:
        - ``pubchem_cid`` (int): PubChem CID
        - ``iupac_name`` (str | None): Preferred IUPAC name
        - ``preferred_name`` (str | None): Preferred name / title
        - ``raw_synonyms`` (List[str]): List of synonym strings
    """
    if not inchikeys:
        return {}

    conn = sqlite3.connect(str(db_path))
    try:
        placeholders = ",".join("?" * len(inchikeys))

        # Fetch core compound data
        compound_rows = conn.execute(
            f"SELECT cid, inchikey, iupac_name, title FROM compounds "
            f"WHERE inchikey IN ({placeholders})",
            inchikeys,
        ).fetchall()

        results: Dict[str, Dict[str, Any]] = {}
        for cid, inchikey, iupac_name, title in compound_rows:
            results[inchikey] = {
                "pubchem_cid": cid,
                "iupac_name": iupac_name,
                "preferred_name": title,
                "raw_synonyms": [],
            }

        # Fetch synonyms for the matched InChIKeys
        synonym_rows = conn.execute(
            f"SELECT c.inchikey, s.synonym_text FROM compounds c "
            f"JOIN synonyms s ON c.cid = s.cid "
            f"WHERE c.inchikey IN ({placeholders})",
            inchikeys,
        ).fetchall()

        for inchikey, synonym_text in synonym_rows:
            if inchikey in results:
                results[inchikey]["raw_synonyms"].append(synonym_text)

        return results
    finally:
        conn.close()
