"""Unit tests for senita_chem.cli module."""

import json
import logging
import sys
from pathlib import Path
from typing import Dict
from unittest.mock import MagicMock, patch

import pytest

from senita_chem.cli import (
    main,
    read_inchikey_file,
    read_smiles_file,
    setup_logging,
    write_results,
)


class TestSetupLogging:
    """Tests for setup_logging function."""

    def test_default_logging_level(self) -> None:
        """setup_logging configures INFO level by default."""
        with patch("senita_chem.cli.logging.basicConfig") as mock_basic_config:
            setup_logging()
            mock_basic_config.assert_called_once()
            assert mock_basic_config.call_args.kwargs["level"] == logging.INFO

    def test_verbose_logging_level(self) -> None:
        """setup_logging configures DEBUG level when verbose=True."""
        with patch("senita_chem.cli.logging.basicConfig") as mock_basic_config:
            setup_logging(verbose=True)
            mock_basic_config.assert_called_once()
            assert mock_basic_config.call_args.kwargs["level"] == logging.DEBUG


class TestReadSmilesFile:
    """Tests for read_smiles_file function."""

    def test_reads_smiles_only(self, tmp_path: Path) -> None:
        """read_smiles_file parses lines with only SMILES."""
        input_file = tmp_path / "smiles.txt"
        input_file.write_text("CCO\nCCC\n")
        result = read_smiles_file(input_file)
        assert result == [
            {"smiles": "CCO", "name": ""},
            {"smiles": "CCC", "name": ""},
        ]

    def test_reads_smiles_with_names(self, tmp_path: Path) -> None:
        """read_smiles_file parses tab-separated SMILES and names."""
        input_file = tmp_path / "smiles.txt"
        input_file.write_text("CCO\tethanol\nCCC\tpropane\n")
        result = read_smiles_file(input_file)
        assert result == [
            {"smiles": "CCO", "name": "ethanol"},
            {"smiles": "CCC", "name": "propane"},
        ]

    def test_skips_empty_lines_and_comments(self, tmp_path: Path) -> None:
        """read_smiles_file skips empty lines and comments."""
        input_file = tmp_path / "smiles.txt"
        input_file.write_text("\n# comment\nCCO\tethanol\n\n")
        result = read_smiles_file(input_file)
        assert result == [{"smiles": "CCO", "name": "ethanol"}]

    def test_warns_on_too_many_tabs(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """read_smiles_file warns when more than one tab is present."""
        input_file = tmp_path / "smiles.txt"
        input_file.write_text("CCO\tethanol\textra\n")
        with caplog.at_level(logging.WARNING):
            result = read_smiles_file(input_file)
        assert result == [{"smiles": "CCO", "name": "ethanol"}]
        assert "Too many tabs" in caplog.text


class TestReadInchikeyFile:
    """Tests for read_inchikey_file function."""

    def test_reads_inchikeys(self, tmp_path: Path) -> None:
        """read_inchikey_file reads one InChIKey per line."""
        input_file = tmp_path / "inchikeys.txt"
        input_file.write_text(
            "ABCDEFGHIJKLMNO-PQRSTUVWXY-Z\n12345678901234-5678901234-A\n"
        )
        result = read_inchikey_file(input_file)
        assert result == [
            "ABCDEFGHIJKLMNO-PQRSTUVWXY-Z",
            "12345678901234-5678901234-A",
        ]

    def test_skips_empty_lines_and_comments(self, tmp_path: Path) -> None:
        """read_inchikey_file skips empty lines and comments."""
        input_file = tmp_path / "inchikeys.txt"
        input_file.write_text("\n# comment\nABCDEFGHIJKLMNO-PQRSTUVWXY-Z\n\n")
        result = read_inchikey_file(input_file)
        assert result == ["ABCDEFGHIJKLMNO-PQRSTUVWXY-Z"]


class TestWriteResults:
    """Tests for write_results function."""

    def test_writes_json(self, tmp_path: Path) -> None:
        """write_results serializes results to JSON file."""
        output_file = tmp_path / "results.json"
        results: Dict = {"key": {"value": 1}}
        write_results(results, output_file)
        assert output_file.exists()
        with open(output_file, "r", encoding="utf-8") as f:
            loaded = json.load(f)
        assert loaded == results


class TestMain:
    """Tests for main CLI entry point."""

    @patch("senita_chem.cli.enrich_compounds")
    @patch(
        "senita_chem.cli.sys.argv", ["senita-chem", "CCO", "--pubchem-method", "api"]
    )
    def test_single_smiles(
        self, mock_enrich: MagicMock, capsys: pytest.CaptureFixture
    ) -> None:
        """main processes a single SMILES string and prints JSON to stdout."""
        mock_enrich.return_value = {
            "inchikey1": {"pubchem_cid": 123, "enrichment_source": "pubchem"}
        }
        main()
        captured = capsys.readouterr()
        assert "inchikey1" in captured.out
        mock_enrich.assert_called_once_with(
            compounds=[{"smiles": "CCO", "name": ""}],
            inchikeys=None,
            max_synonyms=75,
            pubchem_method="api",
            db_path=None,
        )

    @patch("senita_chem.cli.enrich_compounds")
    @patch("senita_chem.cli.sys.argv", ["senita-chem", "--input", "/tmp/smiles.txt"])
    def test_input_file_not_found(
        self, mock_enrich: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """main exits with error when input file does not exist."""
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
        assert "not found" in caplog.text

    @patch("senita_chem.cli.enrich_compounds")
    @patch("senita_chem.cli.sys.argv", ["senita-chem", "--inchikeys", "/tmp/keys.txt"])
    def test_inchikey_file_not_found(
        self, mock_enrich: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """main exits with error when InChIKey file does not exist."""
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
        assert "not found" in caplog.text

    @patch("senita_chem.cli.enrich_compounds")
    @patch(
        "senita_chem.cli.sys.argv",
        ["senita-chem", "CCO", "--pubchem-method", "local_db"],
    )
    def test_local_db_requires_db_path(
        self, mock_enrich: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """main exits with error when local_db is chosen without --db-path."""
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
        assert "--db-path is required" in caplog.text

    @patch("senita_chem.cli.enrich_compounds")
    def test_writes_to_output_file(
        self, mock_enrich: MagicMock, tmp_path: Path
    ) -> None:
        """main writes results to the specified output file."""
        output_file = tmp_path / "results.json"
        mock_enrich.return_value = {
            "inchikey1": {"pubchem_cid": 123, "enrichment_source": "pubchem"}
        }
        with patch.object(
            sys,
            "argv",
            [
                "senita-chem",
                "CCO",
                "--output",
                str(output_file),
                "--pubchem-method",
                "api",
            ],
        ):
            main()
        assert output_file.exists()
        with open(output_file, "r", encoding="utf-8") as f:
            loaded = json.load(f)
        assert loaded == mock_enrich.return_value

    @patch("senita_chem.cli.enrich_compounds")
    @patch(
        "senita_chem.cli.sys.argv",
        ["senita-chem", "CCO", "--verbose", "--pubchem-method", "api"],
    )
    def test_verbose_mode(
        self, mock_enrich: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """main runs successfully in verbose mode."""
        mock_enrich.return_value = {
            "inchikey1": {"pubchem_cid": 123, "enrichment_source": "pubchem"}
        }
        with caplog.at_level(logging.DEBUG):
            main()
        assert "Processing single SMILES" in caplog.text

    @patch("senita_chem.cli.enrich_compounds")
    @patch(
        "senita_chem.cli.sys.argv", ["senita-chem", "CCO", "--pubchem-method", "api"]
    )
    def test_keyboard_interrupt(
        self, mock_enrich: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """main exits gracefully on KeyboardInterrupt."""
        mock_enrich.side_effect = KeyboardInterrupt()
        with caplog.at_level(logging.INFO):
            with pytest.raises(SystemExit) as exc_info:
                main()
        assert exc_info.value.code == 1
        assert "Interrupted by user" in caplog.text

    @patch("senita_chem.cli.enrich_compounds")
    @patch(
        "senita_chem.cli.sys.argv", ["senita-chem", "CCO", "--pubchem-method", "api"]
    )
    def test_unexpected_error(
        self, mock_enrich: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """main exits with error on unexpected exceptions."""
        mock_enrich.side_effect = RuntimeError("boom")
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
        assert "Error: boom" in caplog.text

    @patch("senita_chem.cli.enrich_compounds")
    @patch(
        "senita_chem.cli.sys.argv",
        ["senita-chem", "CCO", "--verbose", "--pubchem-method", "api"],
    )
    def test_unexpected_error_verbose_traceback(
        self, mock_enrich: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """main logs full traceback in verbose mode on unexpected exceptions."""
        mock_enrich.side_effect = RuntimeError("boom")
        with caplog.at_level(logging.DEBUG):
            with pytest.raises(SystemExit):
                main()
        assert "Full traceback:" in caplog.text

    @patch("senita_chem.cli.enrich_compounds")
    def test_custom_max_synonyms(self, mock_enrich: MagicMock) -> None:
        """main passes custom --max-synonyms value to enrich_compounds."""
        with patch.object(
            sys,
            "argv",
            ["senita-chem", "CCO", "--max-synonyms", "10", "--pubchem-method", "api"],
        ):
            mock_enrich.return_value = {
                "inchikey1": {"pubchem_cid": 123, "enrichment_source": "pubchem"}
            }
            main()
        assert mock_enrich.call_args.kwargs["max_synonyms"] == 10

    @patch("senita_chem.cli.enrich_compounds")
    def test_api_pubchem_method(self, mock_enrich: MagicMock) -> None:
        """main passes api method to enrich_compounds."""
        with patch.object(
            sys, "argv", ["senita-chem", "CCO", "--pubchem-method", "api"]
        ):
            mock_enrich.return_value = {
                "inchikey1": {"pubchem_cid": 123, "enrichment_source": "pubchem"}
            }
            main()
        assert mock_enrich.call_args.kwargs["pubchem_method"] == "api"
