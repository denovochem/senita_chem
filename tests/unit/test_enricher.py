"""Unit tests for senita_chem.enricher module."""

from typing import Dict, Iterator
from unittest.mock import patch

import pytest

from senita_chem.enricher import enrich_compounds


@pytest.fixture(autouse=True)
def clear_caches() -> Iterator[None]:
    """Clear module-level LRU caches between tests."""
    from senita_chem.iupac_naming import cache_clear as iupac_cache_clear
    from senita_chem.rdkit_properties import cache_clear as rdkit_cache_clear

    rdkit_cache_clear()
    iupac_cache_clear()
    yield


@pytest.fixture
def mock_rdkit_props() -> Dict:
    """Return a realistic set of RDKit properties for a single compound."""
    return {
        "canonical_smiles": "CCO",
        "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "is_multi_fragment": False,
        "formula": "C2H6O",
        "mw": 46.0684,
        "mw_exact": 46.0419,
        "logp": -0.0014,
        "tpsa": 20.23,
        "hba": 1,
        "hbd": 1,
        "num_heavy_atoms": 3,
        "num_rotatable_bonds": 0,
        "num_rings": 0,
        "num_aromatic_rings": 0,
        "num_aliphatic_rings": 0,
        "num_heterocycles": 0,
        "frac_csp3": 1.0,
        "num_stereocenters": 0,
        "num_defined_stereocenters": 0,
        "num_undefined_stereocenters": 0,
        "formal_charge": 0,
    }


@pytest.fixture
def mock_pubchem_result() -> Dict:
    """Return a realistic PubChem lookup result for a single compound."""
    return {
        "pubchem_cid": "702",
        "iupac_name": "ethanol",
        "preferred_name": "Ethanol",
        "canonical_smiles_pubchem": "CCO",
        "isomeric_smiles_pubchem": "CCO",
        "inchi_pubchem": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "inchikey_pubchem": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "formula_pubchem": "C2H6O",
        "raw_synonyms": ["ethanol", "ethyl alcohol", "grain alcohol"],
    }


class TestEnrichCompoundsValidation:
    """Tests for enrich_compounds input validation."""

    def test_raises_when_neither_input_given(self) -> None:
        """enrich_compounds raises ValueError when both inputs are None."""
        with pytest.raises(ValueError, match="Provide either compounds or inchikeys"):
            enrich_compounds()

    def test_raises_on_invalid_pubchem_method(self, mock_rdkit_props: Dict) -> None:
        """enrich_compounds raises ValueError for unknown pubchem_method."""
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ):
            with pytest.raises(ValueError, match="Invalid pubchem_method"):
                enrich_compounds(
                    compounds=[{"smiles": "CCO", "name": ""}],
                    pubchem_method="unknown",
                )

    def test_raises_when_local_db_missing_db_path(self, mock_rdkit_props: Dict) -> None:
        """enrich_compounds raises ValueError when local_db is used without db_path."""
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ):
            with pytest.raises(ValueError, match="db_path is required"):
                enrich_compounds(
                    compounds=[{"smiles": "CCO", "name": ""}],
                    pubchem_method="local_db",
                )


class TestEnrichCompoundsWithSmiles:
    """Tests for enrich_compounds using SMILES input."""

    def test_single_smiles_api_success(
        self, mock_rdkit_props: Dict, mock_pubchem_result: Dict
    ) -> None:
        """enrich_compounds merges RDKit and PubChem data via API."""
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ) as mock_rdkit:
            with patch(
                "senita_chem.enricher.batch_lookup_by_inchikeys",
                return_value={mock_rdkit_props["inchikey"]: mock_pubchem_result},
            ) as mock_pubchem:
                with patch(
                    "senita_chem.enricher.clean_synonyms_list",
                    return_value=["ethanol", "ethyl alcohol"],
                ) as mock_clean:
                    with patch(
                        "senita_chem.enricher.get_cas_nos_from_synonyms_list",
                        return_value=["64-17-5"],
                    ) as mock_cas:
                        results = enrich_compounds(
                            compounds=[{"smiles": "CCO", "name": "ethanol"}],
                            pubchem_method="api",
                        )

        assert len(results) == 1
        inchikey = mock_rdkit_props["inchikey"]
        assert inchikey in results
        record = results[inchikey]

        assert record["input_smiles"] == "CCO"
        assert record["input_name"] == "ethanol"
        assert record["pubchem_cid"] == "702"
        assert record["synonyms"] == ["ethanol", "ethyl alcohol"]
        assert record["cas"] == ["64-17-5"]
        assert record["enrichment_source"] == "pubchem"

        mock_rdkit.assert_called_once_with("CCO")
        mock_pubchem.assert_called_once_with([inchikey])
        mock_clean.assert_called_once_with(
            ["ethanol", "ethyl alcohol", "grain alcohol"], 75
        )
        mock_cas.assert_called_once_with(["ethanol", "ethyl alcohol", "grain alcohol"])

    def test_single_smiles_local_db_success(
        self, mock_rdkit_props: Dict, mock_pubchem_result: Dict
    ) -> None:
        """enrich_compounds merges RDKit and PubChem data via local_db."""
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ):
            with patch(
                "senita_chem.enricher.batch_lookup_by_inchikeys_sqlite",
                return_value={mock_rdkit_props["inchikey"]: mock_pubchem_result},
            ) as mock_sqlite:
                with patch(
                    "senita_chem.enricher.clean_synonyms_list",
                    return_value=["ethanol"],
                ):
                    with patch(
                        "senita_chem.enricher.get_cas_nos_from_synonyms_list",
                        return_value=[],
                    ):
                        results = enrich_compounds(
                            compounds=[{"smiles": "CCO", "name": ""}],
                            pubchem_method="local_db",
                            db_path="/fake/db.sqlite",
                        )

        assert len(results) == 1
        mock_sqlite.assert_called_once_with(
            [mock_rdkit_props["inchikey"]], db_path="/fake/db.sqlite"
        )

    def test_invalid_smiles_skipped(self, caplog: pytest.LogCaptureFixture) -> None:
        """enrich_compounds skips compounds with invalid SMILES."""
        with patch("senita_chem.enricher.compute_rdkit_properties", return_value=None):
            results = enrich_compounds(
                compounds=[{"smiles": "INVALID", "name": ""}],
                pubchem_method="api",
            )

        assert results == {}
        assert "Invalid SMILES" in caplog.text

    def test_no_pubchem_data_single_fragment(self, mock_rdkit_props: Dict) -> None:
        """enrich_compounds marks single-fragment compounds as 'failed' when PubChem has no data, using OpenClatura for IUPAC name fallback."""
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ):
            with patch(
                "senita_chem.enricher.batch_lookup_by_inchikeys", return_value={}
            ):
                with patch(
                    "senita_chem.enricher.batch_name_smiles",
                    return_value={"CCO": "ethanol"},
                ) as mock_batch:
                    results = enrich_compounds(
                        compounds=[{"smiles": "CCO", "name": ""}],
                        pubchem_method="api",
                    )

        record = results[mock_rdkit_props["inchikey"]]
        assert record["enrichment_source"] == "failed"
        assert record["synonyms"] == ["ethanol"]
        assert record["cas"] == []
        assert record["pubchem_cid"] is None
        assert record["iupac_name"] == "ethanol"
        assert record["preferred_name"] == "ethanol"
        mock_batch.assert_called_once_with(["CCO"], max_workers=None, chunksize=None)

    def test_no_pubchem_data_multi_fragment(self, mock_rdkit_props: Dict) -> None:
        """enrich_compounds marks multi-fragment compounds as 'rdkit_only' when PubChem has no data."""
        mock_rdkit_props["is_multi_fragment"] = True
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ):
            with patch(
                "senita_chem.enricher.batch_lookup_by_inchikeys", return_value={}
            ):
                with patch(
                    "senita_chem.enricher.batch_name_smiles",
                    return_value={"O.CCO": "ethanol"},
                ):
                    results = enrich_compounds(
                        compounds=[{"smiles": "O.CCO", "name": ""}],
                        pubchem_method="api",
                    )

        record = results[mock_rdkit_props["inchikey"]]
        assert record["enrichment_source"] == "rdkit_only"

    def test_max_synonyms_passthrough(
        self, mock_rdkit_props: Dict, mock_pubchem_result: Dict
    ) -> None:
        """enrich_compounds passes max_synonyms to clean_synonyms_list."""
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ):
            with patch(
                "senita_chem.enricher.batch_lookup_by_inchikeys",
                return_value={mock_rdkit_props["inchikey"]: mock_pubchem_result},
            ):
                with patch(
                    "senita_chem.enricher.clean_synonyms_list",
                    return_value=["ethanol"],
                ) as mock_clean:
                    with patch(
                        "senita_chem.enricher.get_cas_nos_from_synonyms_list",
                        return_value=[],
                    ):
                        enrich_compounds(
                            compounds=[{"smiles": "CCO", "name": ""}],
                            pubchem_method="api",
                            max_synonyms=10,
                        )

        mock_clean.assert_called_once_with(
            ["ethanol", "ethyl alcohol", "grain alcohol"], 10
        )

    def test_multiple_compounds(self, mock_rdkit_props: Dict) -> None:
        """enrich_compounds handles multiple compounds."""
        second_props = {**mock_rdkit_props, "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N"}

        def side_effect(smiles: str) -> Dict:
            if smiles == "CCO":
                return mock_rdkit_props
            return second_props

        with patch(
            "senita_chem.enricher.compute_rdkit_properties", side_effect=side_effect
        ):
            with patch(
                "senita_chem.enricher.batch_lookup_by_inchikeys", return_value={}
            ):
                with patch(
                    "senita_chem.enricher.batch_name_smiles",
                    return_value={"CCO": "ethanol", "O": "ethanol"},
                ):
                    results = enrich_compounds(
                        compounds=[
                            {"smiles": "CCO", "name": "ethanol"},
                            {"smiles": "O", "name": "water"},
                        ],
                        pubchem_method="api",
                    )

        assert len(results) == 2
        assert mock_rdkit_props["inchikey"] in results
        assert second_props["inchikey"] in results


class TestEnrichCompoundsWithInchikeys:
    """Tests for enrich_compounds using InChIKey input."""

    def test_inchikey_only_input(self) -> None:
        """enrich_compounds processes inchikeys-only input via API."""
        inchikey = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        pubchem_result = {
            "pubchem_cid": "702",
            "iupac_name": "ethanol",
            "preferred_name": "Ethanol",
            "raw_synonyms": ["ethanol"],
        }

        with patch(
            "senita_chem.enricher.batch_lookup_by_inchikeys",
            return_value={inchikey: pubchem_result},
        ) as mock_pubchem:
            with patch(
                "senita_chem.enricher.clean_synonyms_list",
                return_value=["ethanol"],
            ):
                with patch(
                    "senita_chem.enricher.get_cas_nos_from_synonyms_list",
                    return_value=[],
                ):
                    results = enrich_compounds(
                        inchikeys=[inchikey],
                        pubchem_method="api",
                    )

        assert len(results) == 1
        record = results[inchikey]
        assert record["inchikey"] == inchikey
        assert record["is_multi_fragment"] is False
        assert record["pubchem_cid"] == "702"
        mock_pubchem.assert_called_once_with([inchikey])

    def test_inchikey_only_local_db(self) -> None:
        """enrich_compounds processes inchikeys-only input via local_db."""
        inchikey = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        with patch(
            "senita_chem.enricher.batch_lookup_by_inchikeys_sqlite",
            return_value={},
        ) as mock_sqlite:
            results = enrich_compounds(
                inchikeys=[inchikey],
                pubchem_method="local_db",
                db_path="/fake/db.sqlite",
            )

        assert len(results) == 1
        record = results[inchikey]
        assert record["enrichment_source"] == "failed"
        mock_sqlite.assert_called_once_with([inchikey], db_path="/fake/db.sqlite")

    def test_empty_inchikeys(self) -> None:
        """enrich_compounds returns empty dict for empty inchikeys list."""
        with patch("senita_chem.enricher.batch_lookup_by_inchikeys", return_value={}):
            results = enrich_compounds(
                inchikeys=[],
                pubchem_method="api",
            )
        assert results == {}


class TestEnrichCompoundsEdgeCases:
    """Tests for edge cases in enrich_compounds."""

    def test_empty_compounds(self) -> None:
        """enrich_compounds returns empty dict for empty compounds list."""
        with patch("senita_chem.enricher.compute_rdkit_properties", return_value=None):
            results = enrich_compounds(
                compounds=[],
                pubchem_method="api",
            )
        assert results == {}

    def test_pubchem_result_without_raw_synonyms(self, mock_rdkit_props: Dict) -> None:
        """enrich_compounds handles PubChem results missing raw_synonyms key."""
        pubchem_result = {
            "pubchem_cid": "702",
            "iupac_name": "ethanol",
            "preferred_name": "Ethanol",
        }
        with patch(
            "senita_chem.enricher.compute_rdkit_properties",
            return_value=mock_rdkit_props,
        ):
            with patch(
                "senita_chem.enricher.batch_lookup_by_inchikeys",
                return_value={mock_rdkit_props["inchikey"]: pubchem_result},
            ):
                with patch(
                    "senita_chem.enricher.clean_synonyms_list",
                    return_value=[],
                ) as mock_clean:
                    with patch(
                        "senita_chem.enricher.get_cas_nos_from_synonyms_list",
                        return_value=[],
                    ):
                        enrich_compounds(
                            compounds=[{"smiles": "CCO", "name": ""}],
                            pubchem_method="api",
                        )

        mock_clean.assert_called_once_with([], 75)
