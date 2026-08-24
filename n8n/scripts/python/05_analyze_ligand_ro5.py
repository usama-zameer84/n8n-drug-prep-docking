"""n8n Code node: Analyze Ligand (Ro5).

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    from rdkit import Chem
    from rdkit.Chem import Crippen
    from rdkit.Chem import Descriptors
    from rdkit.Chem import Lipinski
    from rdkit.Chem import rdMolDescriptors

    d = _items[0]["json"]
    mol = Chem.MolFromSmiles(d["canonical_smiles"])
    descriptors = {
        "molecular_weight": round(Descriptors.MolWt(mol), 3),
        "logp": round(Crippen.MolLogP(mol), 3),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "h_bond_donors": Lipinski.NumHDonors(mol),
        "h_bond_acceptors": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
    }
    violations = []
    if descriptors["molecular_weight"] > 500: violations.append("MW > 500")
    if descriptors["logp"] > 5: violations.append("logP > 5")
    if descriptors["h_bond_donors"] > 5: violations.append("HBD > 5")
    if descriptors["h_bond_acceptors"] > 10: violations.append("HBA > 10")
    out = dict(d)
    out["ligand_analysis"] = {
        "descriptors": descriptors,
        "lipinski_rule_of_five_violations": violations,
        "ro5_pass": len(violations) == 0,
        "interpretation_limit": "Rule-of-five descriptors are filters, not efficacy, safety, or binding evidence.",
    }
    return [{"json": out}]
