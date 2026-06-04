"""Unit tests for senita_chem.local_pubchem_batch module."""

from unittest.mock import MagicMock, patch

import pytest

from senita_chem.local_pubchem_batch import batch_lookup_by_inchikeys_sqlite


class TestBatchLookupByInchikeysSqlite:
    """Tests for batch_lookup_by_inchikeys_sqlite function."""

    def test_empty_inchikeys(self) -> None:
        """batch_lookup_by_inchikeys_sqlite returns empty dict for empty input."""
        results = batch_lookup_by_inchikeys_sqlite([], db_path="/fake/db.sqlite")
        assert results == {}

    @patch("senita_chem.local_pubchem_batch.sqlite3.connect")
    def test_single_match_with_synonyms(self, mock_connect: MagicMock) -> None:
        """batch_lookup_by_inchikeys_sqlite returns compound data and synonyms."""
        mock_conn = MagicMock()
        mock_connect.return_value = mock_conn

        # Simulate compound row
        mock_conn.execute.return_value.fetchall.side_effect = [
            [(702, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethanol", "Ethanol")],
            [
                ("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethyl alcohol"),
                ("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "grain alcohol"),
            ],
        ]

        results = batch_lookup_by_inchikeys_sqlite(
            ["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"],
            db_path="/fake/db.sqlite",
        )

        assert len(results) == 1
        record = results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]
        assert record["pubchem_cid"] == 702
        assert record["iupac_name"] == "ethanol"
        assert record["label"] == "Ethanol"
        assert record["raw_synonyms"] == ["ethyl alcohol", "grain alcohol"]

        mock_connect.assert_called_once_with("/fake/db.sqlite")
        mock_conn.close.assert_called_once()

    @patch("senita_chem.local_pubchem_batch.sqlite3.connect")
    def test_single_match_no_synonyms(self, mock_connect: MagicMock) -> None:
        """batch_lookup_by_inchikeys_sqlite returns empty synonyms when none exist."""
        mock_conn = MagicMock()
        mock_connect.return_value = mock_conn

        mock_conn.execute.return_value.fetchall.side_effect = [
            [(702, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethanol", "Ethanol")],
            [],
        ]

        results = batch_lookup_by_inchikeys_sqlite(
            ["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"],
            db_path="/fake/db.sqlite",
        )

        record = results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]
        assert record["raw_synonyms"] == []

    @patch("senita_chem.local_pubchem_batch.sqlite3.connect")
    def test_no_matches(self, mock_connect: MagicMock) -> None:
        """batch_lookup_by_inchikeys_sqlite returns empty dict when no matches."""
        mock_conn = MagicMock()
        mock_connect.return_value = mock_conn

        mock_conn.execute.return_value.fetchall.side_effect = [
            [],
            [],
        ]

        results = batch_lookup_by_inchikeys_sqlite(
            ["UNKNOWN-KEY-12345678-A"],
            db_path="/fake/db.sqlite",
        )

        assert results == {}

    @patch("senita_chem.local_pubchem_batch.sqlite3.connect")
    def test_multiple_inchikeys(self, mock_connect: MagicMock) -> None:
        """batch_lookup_by_inchikeys_sqlite handles multiple InChIKeys."""
        mock_conn = MagicMock()
        mock_connect.return_value = mock_conn

        mock_conn.execute.return_value.fetchall.side_effect = [
            [
                (702, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethanol", "Ethanol"),
                (962, "XLYOFNOQVPJJNP-UHFFFAOYSA-N", "oxidane", "Water"),
            ],
            [
                ("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethyl alcohol"),
                ("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "dihydrogen monoxide"),
            ],
        ]

        results = batch_lookup_by_inchikeys_sqlite(
            [
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            ],
            db_path="/fake/db.sqlite",
        )

        assert len(results) == 2
        assert results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]["raw_synonyms"] == [
            "ethyl alcohol"
        ]
        assert results["XLYOFNOQVPJJNP-UHFFFAOYSA-N"]["raw_synonyms"] == [
            "dihydrogen monoxide"
        ]

    @patch("senita_chem.local_pubchem_batch.sqlite3.connect")
    def test_synonyms_for_unmatched_compound_ignored(
        self, mock_connect: MagicMock
    ) -> None:
        """batch_lookup_by_inchikeys_sqlite ignores synonyms for compounds not in results."""
        mock_conn = MagicMock()
        mock_connect.return_value = mock_conn

        mock_conn.execute.return_value.fetchall.side_effect = [
            [(702, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethanol", "Ethanol")],
            [
                ("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethyl alcohol"),
                ("UNKNOWN-KEY-12345678-A", "orphan synonym"),
            ],
        ]

        results = batch_lookup_by_inchikeys_sqlite(
            ["LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "UNKNOWN-KEY-12345678-A"],
            db_path="/fake/db.sqlite",
        )

        assert len(results) == 1
        assert (
            "orphan synonym"
            not in results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]["raw_synonyms"]
        )

    @patch("senita_chem.local_pubchem_batch.sqlite3.connect")
    def test_connection_closed_on_exception(self, mock_connect: MagicMock) -> None:
        """batch_lookup_by_inchikeys_sqlite closes connection even on exception."""
        mock_conn = MagicMock()
        mock_connect.return_value = mock_conn
        mock_conn.execute.side_effect = RuntimeError("database locked")

        with pytest.raises(RuntimeError, match="database locked"):
            batch_lookup_by_inchikeys_sqlite(
                ["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"],
                db_path="/fake/db.sqlite",
            )

        mock_conn.close.assert_called_once()
