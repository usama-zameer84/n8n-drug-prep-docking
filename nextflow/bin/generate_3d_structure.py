#!/usr/bin/env python3
"""Generate a deterministic three-dimensional ligand conformer."""

import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True, type=Path)
    parser.add_argument("--ligand-sdf", required=True, type=Path)
    parser.add_argument("--ligand-svg", required=True, type=Path)
    args = parser.parse_args()

    data = json.loads(args.params.read_text(encoding="utf-8"))
    molecule = Chem.MolFromSmiles(data["smiles"])
    if molecule is None:
        raise SystemExit("RDKit could not parse the input SMILES")

    molecule = Chem.AddHs(molecule)
    seed = int(data.get("seed", 20260824))
    embedding = AllChem.ETKDGv3()
    embedding.randomSeed = seed
    embedding.useRandomCoords = True
    conformer_id = AllChem.EmbedMolecule(molecule, embedding)
    if conformer_id < 0:
        raise SystemExit("3D embedding failed")

    try:
        AllChem.MMFFOptimizeMolecule(molecule, confId=conformer_id)
        forcefield = "MMFF"
    except Exception:
        AllChem.UFFOptimizeMolecule(molecule, confId=conformer_id)
        forcefield = "UFF"

    molecule_without_hydrogens = Chem.RemoveHs(molecule)
    Chem.MolToMolFile(molecule, str(args.ligand_sdf))
    Draw.MolToFile(molecule_without_hydrogens, str(args.ligand_svg))

    data.update(
        {
            "canonical_smiles": Chem.MolToSmiles(
                molecule_without_hydrogens, isomericSmiles=True
            ),
            "inchi_key": Chem.MolToInchiKey(molecule_without_hydrogens),
            "num_heavy_atoms": molecule.GetNumHeavyAtoms(),
            "num_atoms": molecule.GetNumAtoms(),
            "conformer_generation": {
                "method": "ETKDGv3",
                "forcefield": forcefield,
                "seed": seed,
            },
        }
    )
    args.params.write_text(json.dumps(data, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
