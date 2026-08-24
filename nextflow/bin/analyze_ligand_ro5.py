#!/usr/bin/env python3
"""Calculate ligand descriptors and Lipinski Rule-of-Five status."""

import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True, type=Path)
    parser.add_argument("--ligand-sdf", required=True, type=Path)
    args = parser.parse_args()

    data = json.loads(args.params.read_text(encoding="utf-8"))
    molecule = next(
        iter(Chem.SDMolSupplier(str(args.ligand_sdf), removeHs=False)), None
    )
    if molecule is None:
        raise SystemExit(f"Could not read {args.ligand_sdf}")

    molecule_without_hydrogens = Chem.RemoveHs(molecule)
    molecular_weight = Descriptors.MolWt(molecule_without_hydrogens)
    alogp = Crippen.MolLogP(molecule_without_hydrogens)
    hydrogen_bond_donors = Lipinski.NumHDonors(molecule_without_hydrogens)
    hydrogen_bond_acceptors = Lipinski.NumHAcceptors(molecule_without_hydrogens)
    violations = sum(
        (
            molecular_weight > 500,
            alogp > 5,
            hydrogen_bond_donors > 5,
            hydrogen_bond_acceptors > 10,
        )
    )
    data["ro5"] = {
        "molecular_weight": round(molecular_weight, 3),
        "alogp": round(alogp, 3),
        "h_bond_donors": hydrogen_bond_donors,
        "h_bond_acceptors": hydrogen_bond_acceptors,
        "tpsa": round(rdMolDescriptors.CalcTPSA(molecule_without_hydrogens), 3),
        "rotatable_bonds": Lipinski.NumRotatableBonds(molecule_without_hydrogens),
        "violations": violations,
        "druglike": violations <= 1,
    }
    args.params.write_text(json.dumps(data, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
