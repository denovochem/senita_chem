"""Unit tests for senita_chem.pubchem_batch module."""

from unittest.mock import MagicMock, patch

import pytest

from senita_chem.pubchem_batch import (
    batch_lookup_by_inchikeys,
    direct_rest_lookup_by_inchikeys,
)


class TestDirectRestLookupByInchikeys:
    """Tests for direct_rest_lookup_by_inchikeys function."""

    def test_empty_inchikeys(self) -> None:
        """direct_rest_lookup_by_inchikeys returns empty dict for empty input."""
        results = direct_rest_lookup_by_inchikeys([])
        assert results == {}

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_single_inchikey_success(self, mock_get: MagicMock) -> None:
        """direct_rest_lookup_by_inchikeys parses properties and synonyms."""
        mock_get.side_effect = [
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(
                    return_value={
                        "PropertyTable": {
                            "Properties": [
                                {
                                    "CID": 702,
                                    "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                                    "IUPACName": "ethanol",
                                    "Title": "Ethanol",
                                    "CanonicalSMILES": "CCO",
                                    "IsomericSMILES": "CCO",
                                    "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                                    "MolecularFormula": "C2H6O",
                                }
                            ]
                        }
                    }
                ),
            ),
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(
                    return_value={
                        "InformationList": {
                            "Information": [
                                {
                                    "CID": 702,
                                    "Synonym": ["ethanol", "ethyl alcohol"],
                                }
                            ]
                        }
                    }
                ),
            ),
        ]

        results = direct_rest_lookup_by_inchikeys(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])

        assert len(results) == 1
        record = results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]
        assert record["pubchem_cid"] == "702"
        assert record["iupac_name"] == "ethanol"
        assert record["preferred_name"] == "Ethanol"
        assert record["raw_synonyms"] == ["ethanol", "ethyl alcohol"]

        assert mock_get.call_count == 2
        props_call, synonyms_call = mock_get.call_args_list
        assert "property" in props_call.args[0]
        assert "synonyms" in synonyms_call.args[0]

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_no_properties_response(self, mock_get: MagicMock) -> None:
        """direct_rest_lookup_by_inchikeys returns empty dict when no properties."""
        mock_get.side_effect = [
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(return_value={"PropertyTable": {"Properties": []}}),
            ),
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(return_value={"InformationList": {"Information": []}}),
            ),
        ]

        results = direct_rest_lookup_by_inchikeys(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])
        assert results == {}

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_no_synonyms_for_match(self, mock_get: MagicMock) -> None:
        """direct_rest_lookup_by_inchikeys ensures raw_synonyms defaults to []."""
        mock_get.side_effect = [
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(
                    return_value={
                        "PropertyTable": {
                            "Properties": [
                                {
                                    "CID": 702,
                                    "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                                    "IUPACName": "ethanol",
                                    "Title": "Ethanol",
                                    "CanonicalSMILES": "CCO",
                                    "IsomericSMILES": "CCO",
                                    "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                                    "MolecularFormula": "C2H6O",
                                }
                            ]
                        }
                    }
                ),
            ),
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(return_value={"InformationList": {"Information": []}}),
            ),
        ]

        results = direct_rest_lookup_by_inchikeys(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])

        assert results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]["raw_synonyms"] == []

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_request_error_continues(
        self, mock_get: MagicMock, caplog: pytest.LogCaptureFixture
    ) -> None:
        """direct_rest_lookup_by_inchikeys logs error and continues on request failure."""
        mock_get.side_effect = Exception("Connection timeout")

        results = direct_rest_lookup_by_inchikeys(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])

        assert results == {}
        assert "Error in REST chunk" in caplog.text

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_chunking_at_50(self, mock_get: MagicMock) -> None:
        """direct_rest_lookup_by_inchikeys chunks requests at 50 keys."""
        mock_get.side_effect = [
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(return_value={"PropertyTable": {"Properties": []}}),
            ),
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(return_value={"InformationList": {"Information": []}}),
            ),
        ] * 3  # 150 keys → 3 chunks of 50 each

        keys = [f"KEY{i:04d}-ABCDEFGHIJ-A" for i in range(150)]
        direct_rest_lookup_by_inchikeys(keys)

        assert mock_get.call_count == 6  # 3 chunks × 2 requests each
        # Verify each chunk has at most 50 keys in the URL
        for call in mock_get.call_args_list:
            url = call.args[0]
            inchikey_part = url.split("/inchikey/")[1].split("/")[0]
            assert len(inchikey_part.split(",")) <= 50

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_duplicate_inchikeys_deduplicated(self, mock_get: MagicMock) -> None:
        """batch_lookup_by_inchikeys deduplicates input InChIKeys."""
        mock_get.side_effect = [
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(return_value={"PropertyTable": {"Properties": []}}),
            ),
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(return_value={"InformationList": {"Information": []}}),
            ),
        ]

        batch_lookup_by_inchikeys(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"] * 3)

        props_call = mock_get.call_args_list[0]
        url = props_call.args[0]
        inchikey_part = url.split("/inchikey/")[1].split("/")[0]
        assert inchikey_part == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_cid_without_matching_inchikey_ignored(self, mock_get: MagicMock) -> None:
        """direct_rest_lookup_by_inchikeys ignores synonyms for unmatched CIDs."""
        mock_get.side_effect = [
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(
                    return_value={
                        "PropertyTable": {
                            "Properties": [
                                {
                                    "CID": 702,
                                    "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                                    "IUPACName": "ethanol",
                                    "Title": "Ethanol",
                                    "CanonicalSMILES": "CCO",
                                    "IsomericSMILES": "CCO",
                                    "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                                    "MolecularFormula": "C2H6O",
                                }
                            ]
                        }
                    }
                ),
            ),
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(
                    return_value={
                        "InformationList": {
                            "Information": [
                                {
                                    "CID": 999,
                                    "Synonym": ["orphan"],
                                }
                            ]
                        }
                    }
                ),
            ),
        ]

        results = direct_rest_lookup_by_inchikeys(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])

        assert results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]["raw_synonyms"] == []

    @patch("senita_chem.pubchem_batch.requests.get")
    def test_multiple_cids_per_chunk(self, mock_get: MagicMock) -> None:
        """direct_rest_lookup_by_inchikeys handles multiple compounds per chunk."""
        mock_get.side_effect = [
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(
                    return_value={
                        "PropertyTable": {
                            "Properties": [
                                {
                                    "CID": 702,
                                    "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                                    "IUPACName": "ethanol",
                                    "Title": "Ethanol",
                                    "CanonicalSMILES": "CCO",
                                    "IsomericSMILES": "CCO",
                                    "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                                    "MolecularFormula": "C2H6O",
                                },
                                {
                                    "CID": 962,
                                    "InChIKey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                                    "IUPACName": "oxidane",
                                    "Title": "Water",
                                    "CanonicalSMILES": "O",
                                    "IsomericSMILES": "O",
                                    "InChI": "InChI=1S/H2O/h1H2",
                                    "MolecularFormula": "H2O",
                                },
                            ]
                        }
                    }
                ),
            ),
            MagicMock(
                raise_for_status=MagicMock(),
                json=MagicMock(
                    return_value={
                        "InformationList": {
                            "Information": [
                                {
                                    "CID": 702,
                                    "Synonym": ["ethanol", "ethyl alcohol"],
                                },
                                {
                                    "CID": 962,
                                    "Synonym": ["water", "dihydrogen monoxide"],
                                },
                            ]
                        }
                    }
                ),
            ),
        ]

        results = direct_rest_lookup_by_inchikeys(
            ["LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"]
        )

        assert len(results) == 2
        assert results["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]["raw_synonyms"] == [
            "ethanol",
            "ethyl alcohol",
        ]
        assert results["XLYOFNOQVPJJNP-UHFFFAOYSA-N"]["raw_synonyms"] == [
            "water",
            "dihydrogen monoxide",
        ]


class TestBatchLookupByInchikeys:
    """Tests for batch_lookup_by_inchikeys wrapper function."""

    @patch("senita_chem.pubchem_batch.direct_rest_lookup_by_inchikeys")
    def test_delegates_to_direct_rest(self, mock_direct: MagicMock) -> None:
        """batch_lookup_by_inchikeys delegates to direct_rest_lookup_by_inchikeys."""
        mock_direct.return_value = {"key": {"pubchem_cid": "1"}}
        results = batch_lookup_by_inchikeys(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])
        assert results == {"key": {"pubchem_cid": "1"}}
        mock_direct.assert_called_once_with(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])

    @patch("senita_chem.pubchem_batch.direct_rest_lookup_by_inchikeys")
    def test_deduplicates_keys(self, mock_direct: MagicMock) -> None:
        """batch_lookup_by_inchikeys removes duplicate InChIKeys before lookup."""
        mock_direct.return_value = {}
        batch_lookup_by_inchikeys(["A", "B", "A", "C", "B"])
        assert mock_direct.call_args.args[0] == ["A", "B", "C"]
