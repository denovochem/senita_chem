import os
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from math import ceil
from typing import Dict, List, Optional

from openclatura import name_smiles

_CACHE_MAXSIZE = 65536
_cache: "OrderedDict[str, Optional[str]]" = OrderedDict()


def _cache_put(smiles: str, name: Optional[str]) -> None:
    """
    Insert a SMILES → name mapping into the shared cache, evicting the least
    recently used entry when the cache is full.

    Args:
        smiles (str): The SMILES string to cache.
        name (Optional[str]): The IUPAC name (or None if conversion failed).
    """
    _cache[smiles] = name
    if len(_cache) > _CACHE_MAXSIZE:
        _cache.popitem(last=False)


def cached_name_smiles(smiles: str) -> Optional[str]:
    """
    Generate an IUPAC name from a SMILES string using OpenClatura, with LRU caching.

    Args:
        smiles (str): A canonical SMILES string to convert to an IUPAC name.

    Returns:
        Optional[str]: The IUPAC name if conversion succeeds, otherwise None.

    Note:
        Results are cached in a module-level LRU cache (max 65 536 entries) so
        repeated calls with the same SMILES string avoid redundant calls to the
        underlying ``name_smiles`` function.  Both successful results and
        failures (None) are cached.
    """
    if smiles in _cache:
        _cache.move_to_end(smiles)
        return _cache[smiles]
    try:
        result = name_smiles(smiles)
    except Exception:
        result = None
    _cache_put(smiles, result)
    return result


def cache_clear() -> None:
    """Clear the module-level LRU cache."""
    _cache.clear()


def _name_smiles_safe(smiles: str) -> Optional[str]:
    """
    Wrapper around ``name_smiles`` that returns None on any exception.

    Args:
        smiles (str): A canonical SMILES string to convert to an IUPAC name.

    Returns:
        Optional[str]: The IUPAC name if conversion succeeds, otherwise None.
    """
    try:
        return name_smiles(smiles)
    except Exception:
        return None


def batch_name_smiles(
    smiles_list: List[str],
    max_workers: Optional[int] = None,
    chunksize: Optional[int] = None,
) -> Dict[str, Optional[str]]:
    """
    Generate IUPAC names for a list of SMILES strings in parallel using processes.

    Uses ``ProcessPoolExecutor`` because ``name_smiles`` calls RDKit under the
    hood, making it CPU-bound.  Cache hits are served instantly from the
    module-level LRU cache; only cache misses are dispatched to worker
    processes.  Results from workers are stored back in the cache so subsequent
    calls (single or batch) benefit from prior computations.  Duplicate SMILES
    in the input are collapsed to a single lookup.

    Args:
        smiles_list (List[str]): List of canonical SMILES strings to convert.
        max_workers (Optional[int]): Maximum number of worker processes. When
            None, defaults to ``min(len(misses), os.cpu_count() or 4)``.
        chunksize (Optional[int]): Number of SMILES per chunk submitted to each
            worker.  When None, a reasonable default is computed from the number
            of misses and worker count.

    Returns:
        Dict[str, Optional[str]]: Mapping from each unique SMILES string to its
        IUPAC name (or None if conversion failed).  Duplicate SMILES in the input
        are collapsed to a single key.
    """
    unique_smiles = list(dict.fromkeys(smiles_list))

    if not unique_smiles:
        return {}

    results: Dict[str, Optional[str]] = {}
    misses: List[str] = []

    for smiles in unique_smiles:
        if smiles in _cache:
            _cache.move_to_end(smiles)
            results[smiles] = _cache[smiles]
        else:
            misses.append(smiles)

    if not misses:
        return results

    if max_workers is None:
        max_workers = min(len(misses), os.cpu_count() or 4)

    if chunksize is None:
        chunksize = max(1, ceil(len(misses) / (4 * max_workers)))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        names = list(executor.map(_name_smiles_safe, misses, chunksize=chunksize))

    for smiles, name in zip(misses, names):
        results[smiles] = name
        _cache_put(smiles, name)

    return results
