from typing import Dict, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]

_taut_opts = rdMolStandardize.CleanupParameters()
_taut_opts.tautomerRemoveSp3Stereo = False
_taut_opts.tautomerRemoveBondStereo = False
TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator(_taut_opts)


def compute_rdkit_properties(smiles: str) -> Optional[Dict]:
    """
    Computes all RDKit physicochemical properties from a SMILES string.

    Returns None if the SMILES is invalid.
    """
    if not smiles or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Stereocenters must be captured before canonicalization strips them
    num_stereocenters = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    num_undefined_stereocenters = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(
        mol
    )

    mol = TAUTOMER_ENUMERATOR.Canonicalize(mol)

    inchi = MolToInchi(mol)
    inchikey = InchiToInchiKey(inchi) if inchi else None
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    is_multi_fragment = "." in canonical_smiles

    return {
        "canonical_smiles": canonical_smiles,
        "inchi": inchi,
        "inchikey": inchikey,
        "is_multi_fragment": is_multi_fragment,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 4),  # type: ignore[attr-defined]
        "mw_exact": round(Descriptors.ExactMolWt(mol), 4),  # type: ignore[attr-defined]
        "logp": round(Descriptors.MolLogP(mol), 4),  # type: ignore[attr-defined]
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 4),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "num_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "num_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "num_heterocycles": rdMolDescriptors.CalcNumHeterocycles(mol),
        "frac_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 4),
        "num_stereocenters": num_stereocenters,
        "num_defined_stereocenters": num_stereocenters,
        "num_undefined_stereocenters": num_undefined_stereocenters,
        "formal_charge": Chem.GetFormalCharge(mol),
    }
