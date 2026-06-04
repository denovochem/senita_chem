#!/usr/bin/env python3
"""
Command-line interface for senita_chem
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List

from senita_chem import enrich_compounds


def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_smiles_file(input_file: Path) -> List[Dict]:
    """
    Read SMILES from file.

    Format: one SMILES per line, optional tab-separated name
    Example: "COC(=O)c1ccc(Cc2ccccc2)cc1\tmethyl 4-benzylbenzoate"
    """
    compounds = []
    with open(input_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) == 1:
                smiles = parts[0].strip()
                name = ""
            elif len(parts) == 2:
                smiles, name = parts[0].strip(), parts[1].strip()
            else:
                logging.warning(
                    f"Line {line_num}: Too many tabs, using first two parts only"
                )
                smiles, name = parts[0].strip(), parts[1].strip()

            if smiles:
                compounds.append({"smiles": smiles, "name": name})
            else:
                logging.warning(f"Line {line_num}: Empty SMILES, skipping")

    return compounds


def read_inchikey_file(input_file: Path) -> List[str]:
    """Read InChIKeys from file (one per line)."""
    inchikeys = []
    with open(input_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            inchikeys.append(line)
    return inchikeys


def write_results(results: Dict, output_file: Path) -> None:
    """Write results to JSON file."""
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)


def main() -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Batch enrichment of chemical compounds using senita_chem",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single compound
  senita-chem "COC(=O)c1ccc(Cc2ccccc2)cc1"
  
  # From file (SMILES + optional names)
  senita-chem --input smiles.txt --output results.json
  
  # From InChIKey list
  senita-chem --inchikeys keys.txt --output results.json
        """,
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("smiles", nargs="?", help="Single SMILES string to enrich")
    input_group.add_argument(
        "--input",
        "-i",
        type=Path,
        help="Input file with SMILES (one per line, optional tab-separated name)",
    )
    input_group.add_argument(
        "--inchikeys", type=Path, help="Input file with InChIKeys (one per line)"
    )

    # Output options
    parser.add_argument(
        "--output", "-o", type=Path, help="Output JSON file (default: stdout)"
    )

    # Processing options
    parser.add_argument(
        "--max-synonyms",
        type=int,
        default=75,
        help="Maximum number of synonyms to keep per compound (default: 75)",
    )

    # Logging options
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Setup logging
    setup_logging(args.verbose)

    try:
        # Get input data
        if args.smiles:
            compounds = [{"smiles": args.smiles, "name": ""}]
            inchikeys = None
            logging.info(f"Processing single SMILES: {args.smiles}")
        elif args.input:
            if not args.input.exists():
                logging.error(f"Input file not found: {args.input}")
                sys.exit(1)
            compounds = read_smiles_file(args.input)
            inchikeys = None
            logging.info(f"Read {len(compounds)} compounds from {args.input}")
        elif args.inchikeys:
            if not args.inchikeys.exists():
                logging.error(f"InChIKey file not found: {args.inchikeys}")
                sys.exit(1)
            inchikeys = read_inchikey_file(args.inchikeys)
            compounds = None
            logging.info(f"Read {len(inchikeys)} InChIKeys from {args.inchikeys}")

        # Enrich compounds
        logging.info("Starting enrichment...")
        results = enrich_compounds(
            compounds=compounds, inchikeys=inchikeys, max_synonyms=args.max_synonyms
        )

        # Output results
        if args.output:
            write_results(results, args.output)
            logging.info(f"Results written to {args.output}")
        else:
            print(json.dumps(results, indent=2, ensure_ascii=False))

        # Summary
        successful = sum(1 for r in results.values() if r.get("pubchem_cid"))
        rdkit_only = sum(
            1 for r in results.values() if r.get("enrichment_source") == "rdkit_only"
        )
        failed = sum(
            1 for r in results.values() if r.get("enrichment_source") == "failed"
        )

        logging.info("Processing complete:")
        logging.info(f"  Total compounds: {len(results)}")
        logging.info(f"  PubChem enriched: {successful}")
        logging.info(f"  RDKit only: {rdkit_only}")
        logging.info(f"  Failed: {failed}")

    except KeyboardInterrupt:
        logging.info("Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error: {e}")
        if args.verbose:
            logging.exception("Full traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()
