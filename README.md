# senita_chem
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://gitHub.com/denovochem/senita_chem/graphs/commit-activity)
[![License](https://img.shields.io/github/license/denovochem/senita_chem)](https://github.com/denovochem/senita_chem/blob/main/LICENSE)
[![Run Tests](https://img.shields.io/github/actions/workflow/status/denovochem/senita_chem/tests.yml?logo=github&logoColor=%23ffffff&label=tests)](https://github.com/denovochem/senita_chem/actions/workflows/tests.yml)
[![Build Docs](https://img.shields.io/github/actions/workflow/status/denovochem/senita_chem/docs.yml?logo=github&logoColor=%23ffffff&label=docs)](https://denovochem.github.io/senita_chem/)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/denovochem/senita_chem/blob/main/examples/example_notebook.ipynb)

Batch enrichment of chemical compounds from SMILES, InChIKeys, or chemical names.

## Description

`senita_chem` is a pip-installable Python library and CLI tool for batch enrichment of chemical compounds. Given a list of SMILES strings, InChIKeys, or chemical names, it returns:

- RDKit-computed physicochemical properties (local, zero network)
- PubChem-sourced identity metadata: names, cleaned synonyms, CAS numbers

PubChem data can be retrieved via the public REST API or queried from a local SQLite database. Input compounds are deduplicated by InChIKey for efficient batch processing.

## Installation

```bash
pip install git+https://github.com/denovochem/senita_chem.git
```

## Quick Start

### Python API

```python
from senita_chem import enrich_compounds

# Via PubChem REST API
results = enrich_compounds(
    compounds=[
        {"smiles": "COC(=O)c1ccc(Cc2ccccc2)cc1", "name": ""},
        {"smiles": "[Br-].[Br-].[Mg+2]", "name": "MgBr2"},
        {"smiles": "Fc1nc(F)c(F)c(SCc2ccccc2)c1F", "name": "4-(benzylthio)-2,3,5,6-tetrafluoropyridine"},
    ],
    pubchem_method="api",
)

# Via local SQLite database
results = enrich_compounds(
    compounds=[{"smiles": "CCO", "name": "ethanol"}],
    pubchem_method="local_db",
    db_path="/path/to/pubchem.sqlite",
)

print(results)
```

### Command Line

```bash
# Single compound via API
senita-chem "COC(=O)c1ccc(Cc2ccccc2)cc1" --pubchem-method api

# From file via local SQLite database (default)
senita-chem --input smiles.txt --output results.json --db-path /path/to/pubchem.sqlite

# From InChIKey list via API
senita-chem --inchikeys keys.txt --output results.json --pubchem-method api
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
results = enrich_compounds(
    inchikeys=["UHOVQNZJYSORNB-UHFFFAOYSA-N"],
    pubchem_method="api",
)
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
3. **PubChem lookup** - Batch lookup via REST API or local SQLite database
4. **Merge results** - Combine RDKit and PubChem data

## Dependencies

- `rdkit>=2023.9.1` - Chemical informatics and property calculation
- `requests>=2.31.0` - HTTP client for PubChem REST API
- `sqlite3` - Local PubChem SQLite database queries

## CLI Options

```
usage: senita-chem [-h] [--input INPUT] [--inchikeys INCHIKEYS]
                   [--output OUTPUT] [--max-synonyms MAX_SYNONYMS]
                   [--pubchem-method {local_db,api}] [--db-path DB_PATH]
                   [--verbose]
                   [smiles]

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
  --pubchem-method {local_db,api}
                        PubChem lookup method: local_db or api (default: local_db)
  --db-path DB_PATH     Path to local PubChem SQLite database (required for local_db)
  --verbose, -v         Enable verbose logging
```

## License

MIT License - see LICENSE file for details.
