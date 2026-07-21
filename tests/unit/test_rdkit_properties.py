"""Unit tests for senita_chem.rdkit_properties module."""

import pytest

from senita_chem.rdkit_properties import compute_rdkit_properties


class TestComputeRdkitProperties:
    """Tests for compute_rdkit_properties function."""

    def test_valid_smiles_ethanol(self) -> None:
        """compute_rdkit_properties returns expected keys and values for ethanol."""
        result = compute_rdkit_properties("CCO")

        assert result is not None
        assert result["canonical_smiles"] == "CCO"
        assert result["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert result["is_multi_fragment"] is False
        assert result["formula"] == "C2H6O"
        assert result["mw"] == pytest.approx(46.0684, abs=0.01)
        assert result["num_heavy_atoms"] == 3
        assert result["num_rings"] == 0
        assert result["formal_charge"] == 0
        assert "inchi" in result
        assert result["inchi"] is not None

    def test_valid_smiles_water(self) -> None:
        """compute_rdkit_properties returns expected values for water."""
        result = compute_rdkit_properties("O")

        assert result is not None
        assert result["canonical_smiles"] == "O"
        assert result["formula"] == "H2O"
        assert result["mw"] == pytest.approx(18.015, abs=0.01)
        assert result["num_heavy_atoms"] == 1

    def test_invalid_smiles(self) -> None:
        """compute_rdkit_properties returns None for invalid SMILES."""
        result = compute_rdkit_properties("NOT_A_SMILES")
        assert result is None

    def test_empty_smiles(self) -> None:
        """compute_rdkit_properties returns None for empty string."""
        result = compute_rdkit_properties("")
        assert result is None

    def test_multi_fragment(self) -> None:
        """compute_rdkit_properties detects multi-fragment SMILES."""
        result = compute_rdkit_properties("O.CCO")

        assert result is not None
        assert result["is_multi_fragment"] is True
        assert "." in result["canonical_smiles"]

    def test_single_fragment(self) -> None:
        """compute_rdkit_properties detects single-fragment SMILES."""
        result = compute_rdkit_properties("c1ccccc1")

        assert result is not None
        assert result["is_multi_fragment"] is False

    def test_tautomer_canonicalization(self) -> None:
        """compute_rdkit_properties canonicalizes tautomers to same InChIKey."""
        result1 = compute_rdkit_properties("CC(=O)O")
        result2 = compute_rdkit_properties("C=C(O)O")

        assert result1 is not None
        assert result2 is not None
        assert result1["inchikey"] == result2["inchikey"]

    def test_aromatic_ring_detection(self) -> None:
        """compute_rdkit_properties counts aromatic rings correctly."""
        result = compute_rdkit_properties("c1ccccc1")

        assert result is not None
        assert result["num_rings"] == 1
        assert result["num_aromatic_rings"] == 1
        assert result["num_aliphatic_rings"] == 0

    def test_aliphatic_ring_detection(self) -> None:
        """compute_rdkit_properties counts aliphatic rings correctly."""
        result = compute_rdkit_properties("C1CCCCC1")

        assert result is not None
        assert result["num_rings"] == 1
        assert result["num_aromatic_rings"] == 0
        assert result["num_aliphatic_rings"] == 1

    def test_heterocycle_detection(self) -> None:
        """compute_rdkit_properties counts heterocycles correctly."""
        result = compute_rdkit_properties("c1ccncc1")

        assert result is not None
        assert result["num_heterocycles"] == 1

    def test_stereochemistry(self) -> None:
        """compute_rdkit_properties handles stereochemistry."""
        result = compute_rdkit_properties("C[C@@H](F)Cl")

        assert result is not None
        assert result["num_stereocenters"] == 1

    def test_charged_molecule(self) -> None:
        """compute_rdkit_properties computes formal charge correctly."""
        result = compute_rdkit_properties("[NH4+]")

        assert result is not None
        assert result["formal_charge"] == 1

    def test_returns_all_expected_keys(self) -> None:
        """compute_rdkit_properties returns all expected keys in the result dict."""
        result = compute_rdkit_properties("CCO")

        expected_keys = {
            "canonical_smiles",
            "inchi",
            "inchikey",
            "is_multi_fragment",
            "formula",
            "mw",
            "mw_exact",
            "logp",
            "tpsa",
            "hba",
            "hbd",
            "num_heavy_atoms",
            "num_rotatable_bonds",
            "num_rings",
            "num_aromatic_rings",
            "num_aliphatic_rings",
            "num_heterocycles",
            "frac_csp3",
            "num_stereocenters",
            "num_defined_stereocenters",
            "num_undefined_stereocenters",
            "formal_charge",
            "tautomer_inchikey",
            "tautomer_smiles",
            "original_canonical_smiles",
            "original_inchikey",
        }

        assert result is not None
        assert set(result.keys()) == expected_keys
