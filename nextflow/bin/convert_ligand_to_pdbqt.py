#!/usr/bin/env python3
"""Prepare the ligand in AutoDock PDBQT format with Meeko."""

import argparse
import json
from pathlib import Path

import meeko
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True, type=Path)
    parser.add_argument("--ligand-sdf", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    data = json.loads(args.params.read_text(encoding="utf-8"))
    molecule = next(
        iter(Chem.SDMolSupplier(str(args.ligand_sdf), removeHs=False)), None
    )
    if molecule is None:
        raise SystemExit(f"Cannot read ligand SDF: {args.ligand_sdf}")

    setups = MoleculePreparation().prepare(molecule)
    if len(setups) != 1:
        raise SystemExit("Meeko must produce exactly one ligand setup")

    result = PDBQTWriterLegacy.write_string(setups[0])
    pdbqt_text = result[0] if isinstance(result, tuple) else result
    if isinstance(result, tuple) and not result[1]:
        raise SystemExit(f"Meeko ligand preparation failed: {result[2]}")
    if "REMARK SMILES" not in pdbqt_text:
        raise SystemExit("Meeko PDBQT is missing the SMILES atom-mapping remarks")

    args.out.write_text(pdbqt_text, encoding="utf-8")
    data["ligand_prep"] = {
        "method": "Meeko MoleculePreparation/PDBQTWriterLegacy",
        "meeko_version": str(getattr(meeko, "__version__", "unknown")),
        "charge_model": "Meeko default Gasteiger",
        "smiles_mapping_present": True,
    }
    args.params.write_text(json.dumps(data, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
