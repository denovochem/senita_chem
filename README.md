# senita_chem

Batch enrichment of chemical compounds from SMILES, InChIKeys, or chemical names.

## Description

`senita_chem` is a pip-installable Python library and CLI tool for batch enrichment of chemical compounds. Given a list of SMILES strings, InChIKeys, or chemical names, it returns:

- RDKit-computed physicochemical properties (local, zero network)
- PubChem-sourced identity metadata: names, cleaned synonyms, CAS numbers (batched via PUG API)

It is designed for scale — millions of compounds — using batched async PubChem requests, deduplication by InChIKey, local caching, and resumable processing.

## Installation

```bash
pip install git+https://github.com/denovochem/senita_chem.git
```

## Quick Start

### Python API

```python
from senita_chem import enrich_compounds

results = enrich_compounds(
    compounds=[
        {"smiles": "COC(=O)c1ccc(Cc2ccccc2)cc1", "name": ""},
        {"smiles": "[Br-].[Br-].[Mg+2]", "name": "MgBr2"},
        {"smiles": "Fc1nc(F)c(F)c(SCc2ccccc2)c1F", "name": "4-(benzylthio)-2,3,5,6-tetrafluoropyridine"},
    ]
)

print(results)
```

### Command Line

```bash
# Single compound
senita-chem "COC(=O)c1ccc(Cc2ccccc2)cc1"

# From file (one SMILES per line, optional tab-separated name)
senita-chem --input smiles.txt --output results.json

# From InChIKey list
senita-chem --inchikeys keys.txt --output results.json
```

## Input Formats

### Compounds List
```python
compounds = [
    {"smiles": "COC(=O)c1ccc(Cc2ccccc2)cc1", "name": ""},
    {"smiles": "[Br-].[Br-].[Mg+2]", "name": "MgBr2"},
]
```

### InChIKeys List
```python
results = enrich_compounds(inchikeys=["UHOVQNZJYSORNB-UHFFFAOYSA-N"])
```

### File Formats

**SMILES file (smiles.txt):**
```
COC(=O)c1ccc(Cc2ccccc2)cc1	methyl 4-benzylbenzoate
[Br-].[Br-].[Mg+2]	MgBr2
Fc1nc(F)c(F)c(SCc2ccccc2)c1F	4-(benzylthio)-2,3,5,6-tetrafluoropyridine
```

**InChIKey file (keys.txt):**
```
XBDQKXXYIPTUBI-UHFFFAOYSA-N
UHOVQNZJYSORNB-UHFFFAOYSA-N
```

## Output Format

Results are returned as a dictionary keyed by InChIKey:

```python
{
    "XBDQKXXYIPTUBI-UHFFFAOYSA-N": {
        # --- Identity (PubChem) ---
        "pubchem_cid": "7005",
        "iupac_name": "methyl 4-benzylbenzoate",
        "preferred_name": "Methyl 4-benzylbenzoate",
        "synonyms": ["methyl 4-benzylbenzoate", "4-benzylbenzoic acid methyl ester", ...],
        "cas": ["7148-03-0"],

        # --- Structure ---
        "inchikey": "XBDQKXXYIPTUBI-UHFFFAOYSA-N",
        "inchi": "InChI=1S/C15H14O2/c1-17-15(16)13-9-7-12(8-10-13)11-14-5-3-2-4-6-14/h2-10H,11H2,1H3",
        "canonical_smiles": "COC(=O)c1ccc(Cc2ccccc2)cc1",

        # --- RDKit properties ---
        "formula": "C15H14O2",
        "mw": 226.27,
        "mw_exact": 226.0994,
        "logp": 3.42,
        "tpsa": 26.3,
        "hba": 2,
        "hbd": 0,
        "num_heavy_atoms": 17,
        "num_rotatable_bonds": 4,
        "num_rings": 2,
        "num_aromatic_rings": 2,
        "num_aliphatic_rings": 0,
        "num_heterocycles": 0,
        "frac_csp3": 0.13,
        "num_stereocenters": 0,
        "num_defined_stereocenters": 0,
        "num_undefined_stereocenters": 0,
        "formal_charge": 0,

        # --- Metadata ---
        "is_multi_fragment": False,
        "enrichment_source": "pubchem",
        "input_smiles": "COC(=O)c1ccc(Cc2ccccc2)cc1",
        "input_name": "",
    }
}
```

## Processing Pipeline

1. **RDKit pass** - Compute properties and detect multi-fragment compounds
2. **Deduplication** - Group by InChIKey for efficient processing
3. **PubChem lookup** - Batch async requests for single-fragment compounds
4. **Multi-fragment resolution** - Use cholla_chem for named multi-fragment compounds
5. **Merge results** - Combine RDKit and PubChem data

## Dependencies

- `rdkit>=2023.9.1` - Chemical informatics and property calculation
- `requests>=2.31.0` - HTTP client for PubChem API
- `tenacity>=8.2.0` - Retry logic for resilient API calls
- `cholla_chem` - Name to SMILES resolution (for multi-fragment compounds)

## CLI Options

```
usage: senita-chem [-h] [--input INPUT] [--inchikeys INCHIKEYS] 
                   [--output OUTPUT] [--max-synonyms MAX_SYNONYMS] 
                   [--verbose] [smiles]

Batch enrichment of chemical compounds using senita_chem

positional arguments:
  smiles                Single SMILES string to enrich

optional arguments:
  -h, --help            show this help message and exit
  --input, -i INPUT     Input file with SMILES (one per line, optional tab-separated name)
  --inchikeys INCHIKEYS Input file with InChIKeys (one per line)
  --output, -o OUTPUT   Output JSON file (default: stdout)
  --max-synonyms MAX_SYNONYMS
                        Maximum number of synonyms to keep per compound (default: 75)
  --verbose, -v         Enable verbose logging
```

## Relationship to cholla_chem

| Library | Direction | Input | Output |
|---|---|---|---|
| `cholla_chem` | Name → Structure | Chemical name | SMILES |
| `senita_chem` | Structure → Identity | SMILES / InChIKey / Name | Properties + names + CAS |

For multi-fragment SMILES (ionic compounds, mixtures) where the structure alone is ambiguous, `senita_chem` accepts an optional `name` alongside the SMILES and uses `cholla_chem` to resolve the name to a canonical single-fragment SMILES before performing the PubChem lookup.

## License

MIT License - see LICENSE file for details.
