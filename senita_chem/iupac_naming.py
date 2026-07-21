from functools import lru_cache
from typing import Optional

from openclatura import name_smiles


@lru_cache(maxsize=65536)
def cached_name_smiles(smiles: str) -> Optional[str]:
    """
    Generate an IUPAC name from a SMILES string using OpenClatura, with LRU caching.

    Args:
        smiles (str): A canonical SMILES string to convert to an IUPAC name.

    Returns:
        Optional[str]: The IUPAC name if conversion succeeds, otherwise None.

    Note:
        Results are cached at the module level so repeated calls with the same
        SMILES string avoid redundant calls to the underlying ``name_smiles``
        function.  Calls that raise are not cached by ``lru_cache``.
    """
    return name_smiles(smiles)
