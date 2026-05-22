from typing import Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey


def compute_rdkit_properties(smiles: str) -> Optional[Dict]:
    """
    Computes all RDKit physicochemical properties from a SMILES string.

    Returns None if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    inchi = MolToInchi(mol)
    inchikey = InchiToInchiKey(inchi) if inchi else None
    canonical_smiles = Chem.MolToSmiles(mol)
    is_multi_fragment = "." in canonical_smiles

    return {
        "canonical_smiles": canonical_smiles,
        "inchi": inchi,
        "inchikey": inchikey,
        "is_multi_fragment": is_multi_fragment,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 4),
        "mw_exact": round(Descriptors.ExactMolWt(mol), 4),
        "logp": round(Descriptors.MolLogP(mol), 4),
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
        "num_stereocenters": rdMolDescriptors.CalcNumStereocenters(mol),
        "num_defined_stereocenters": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "num_undefined_stereocenters": rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
    }
