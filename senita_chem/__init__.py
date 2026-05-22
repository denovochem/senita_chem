"""
senita_chem — Batch enrichment of chemical compounds

Given SMILES strings, InChIKeys, or chemical names, returns:
- RDKit-computed physicochemical properties (local, zero network)
- PubChem-sourced identity metadata: names, cleaned synonyms, CAS numbers

Designed for scale — millions of compounds — using batched async PubChem requests,
deduplication by InChIKey, local caching, and resumable processing.
"""

from .enricher import enrich_compounds

__version__ = "0.1.0"
__all__ = ["enrich_compounds"]
