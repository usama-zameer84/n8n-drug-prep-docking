#!/usr/bin/env python3
"""Analyze heavy-atom contacts for the selected Vina pose."""

import argparse
import json
import math


TYPE_TO_ELEMENT = {
    "C": "C",
    "A": "C",
    "N": "N",
    "NA": "N",
    "NS": "N",
    "O": "O",
    "OA": "O",
    "OS": "O",
    "S": "S",
    "SA": "S",
    "P": "P",
    "F": "F",
    "Cl": "Cl",
    "Br": "Br",
    "I": "I",
    "HD": "H",
    "H": "H",
    "Mg": "Mg",
    "Mn": "Mn",
    "Zn": "Zn",
    "Ca": "Ca",
    "Fe": "Fe",
}


def atom_from_line(line):
    parts = line.split()
    atom_type = parts[-1] if parts else ""
    name = line[12:16].strip()
    return {
        "name": name,
        "resname": line[17:20].strip(),
        "chain": line[21:22].strip() or "_",
        "resseq": line[22:27].strip(),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "atom_type": atom_type,
        "element": TYPE_TO_ELEMENT.get(atom_type, line[76:78].strip() or name[:1]),
    }


def read_atoms(path, first_model=False):
    """Read PDBQT atoms, optionally stopping after the first MODEL."""
    atoms = []
    model_started = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                if first_model and model_started:
                    break
                model_started = True
                continue
            if first_model and model_started and line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                atoms.append(atom_from_line(line))
            except (ValueError, IndexError):
                continue
    return atoms


def analyze(receptor, ligand, cutoff):
    by_residue = {}
    ligand_heavy = [atom for atom in ligand if atom["element"] != "H"]
    receptor_heavy = [atom for atom in receptor if atom["element"] != "H"]
    for ligand_atom in ligand_heavy:
        for receptor_atom in receptor_heavy:
            distance = math.sqrt(
                (ligand_atom["x"] - receptor_atom["x"]) ** 2
                + (ligand_atom["y"] - receptor_atom["y"]) ** 2
                + (ligand_atom["z"] - receptor_atom["z"]) ** 2
            )
            if distance > cutoff:
                continue
            if distance < 2.2:
                category = "very_close_contact_check_for_clash"
            elif (
                ligand_atom["element"] in ("N", "O", "S")
                and receptor_atom["element"] in ("N", "O", "S")
                and distance <= 3.5
            ):
                category = "polar_distance_candidate_requires_angle_check"
            else:
                category = "heavy_atom_distance_contact"
            residue_key = ":".join(
                (
                    receptor_atom["resname"],
                    receptor_atom["chain"],
                    receptor_atom["resseq"],
                )
            )
            by_residue.setdefault(residue_key, []).append(
                {
                    "category": category,
                    "distance_A": round(distance, 3),
                    "ligand_atom": ligand_atom["name"],
                    "ligand_element": ligand_atom["element"],
                    "receptor_atom": receptor_atom["name"],
                    "receptor_element": receptor_atom["element"],
                }
            )

    rows = []
    for residue, contacts in by_residue.items():
        contacts.sort(key=lambda item: item["distance_A"])
        rows.append(
            {
                "residue": residue,
                "minimum_distance_A": contacts[0]["distance_A"],
                "contact_count": len(contacts),
                "contacts": contacts[:20],
            }
        )
    rows.sort(key=lambda item: item["minimum_distance_A"])
    return ligand_heavy, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True)
    args = parser.parse_args()
    with open(args.params, encoding="utf-8") as fh:
        data = json.load(fh)

    cutoff = float(data.get("cutoff", 4.5))
    ligand = read_atoms("docked.pdbqt", first_model=True)
    receptor = read_atoms("receptor.pdbqt")
    if not ligand:
        raise SystemExit("Best Vina pose contains no atoms")
    if not receptor:
        raise SystemExit("Docking receptor contains no atoms")
    ligand_heavy, rows = analyze(receptor, ligand, cutoff)
    summary = {
        "status": "CONTACTS_DETECTED_NOT_VALIDATED",
        "method": "Heavy-atom Euclidean distance screening on the selected rigid-receptor Vina pose.",
        "cutoff_A": cutoff,
        "ligand_docked_atom_count": len(ligand_heavy),
        "contact_residue_count": len(rows),
        "residues": rows[:40],
        "limitations": [
            "Polar candidates are not hydrogen-bond assignments because donor state and D-H-A angle are not evaluated.",
            "No water bridges, pi stacking geometry, metal coordination, entropy, dynamics, or free-energy calculation is included.",
            "Very close contacts must be inspected for steric clashes and preparation errors.",
        ],
    }
    data["interactions_summary"] = summary
    with open("interactions.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    with open(args.params, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


if __name__ == "__main__":
    main()
