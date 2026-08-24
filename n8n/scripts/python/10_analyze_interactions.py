"""n8n Code node: Analyze Interactions.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import json
    import math
    import os

    d = _items[0]["json"]
    out_dir = d["output_dir"]
    cutoff = float(d["cutoff"])

    type_to_element = {
        "C": "C", "A": "C", "N": "N", "NA": "N", "NS": "N",
        "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S",
        "P": "P", "F": "F", "Cl": "Cl", "Br": "Br", "I": "I",
        "HD": "H", "H": "H", "Mg": "Mg", "Mn": "Mn", "Zn": "Zn", "Ca": "Ca", "Fe": "Fe",
    }

    def atom_from_line(line):
        parts = line.split()
        ad_type = parts[-1] if parts else ""
        name = line[12:16].strip()
        return {
            "name": name,
            "resname": line[17:20].strip(),
            "chain": line[21:22].strip() or "_",
            "resseq": line[22:27].strip(),
            "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
            "ad_type": ad_type,
            "element": type_to_element.get(ad_type, (line[76:78].strip() or name[:1])),
        }

    def read_atoms(path, first_model=False):
        atoms = []
        in_model = not first_model
        with open(path) as fh:
            for line in fh:
                if first_model and line.startswith("MODEL"):
                    in_model = True
                    continue
                if first_model and line.startswith("ENDMDL"):
                    break
                if in_model and line.startswith(("ATOM", "HETATM")):
                    try:
                        atoms.append(atom_from_line(line))
                    except (ValueError, IndexError):
                        pass
        return atoms

    receptor = read_atoms(os.path.join(out_dir, "receptor.pdbqt"))
    ligand = read_atoms(os.path.join(out_dir, "docked.pdbqt"), first_model=True)
    if not ligand:
        raise ValueError("Best Vina pose contains no atoms")

    by_residue = {}
    for lig_atom in ligand:
        if lig_atom["element"] == "H":
            continue
        for rec_atom in receptor:
            if rec_atom["element"] == "H":
                continue
            distance = math.sqrt(
                (lig_atom["x"] - rec_atom["x"]) ** 2 +
                (lig_atom["y"] - rec_atom["y"]) ** 2 +
                (lig_atom["z"] - rec_atom["z"]) ** 2
            )
            if distance > cutoff:
                continue
            if distance < 2.2:
                category = "very_close_contact_check_for_clash"
            elif lig_atom["element"] in ("N", "O", "S") and rec_atom["element"] in ("N", "O", "S") and distance <= 3.5:
                category = "polar_distance_candidate_requires_angle_check"
            else:
                category = "heavy_atom_distance_contact"
            residue_key = rec_atom["resname"] + ":" + rec_atom["chain"] + ":" + rec_atom["resseq"]
            by_residue.setdefault(residue_key, []).append({
                "category": category,
                "distance_A": round(distance, 3),
                "ligand_atom": lig_atom["name"],
                "ligand_element": lig_atom["element"],
                "receptor_atom": rec_atom["name"],
                "receptor_element": rec_atom["element"],
            })

    rows = []
    for residue, contacts in by_residue.items():
        contacts.sort(key=lambda item: item["distance_A"])
        rows.append({
            "residue": residue,
            "minimum_distance_A": contacts[0]["distance_A"],
            "contact_count": len(contacts),
            "contacts": contacts[:20],
        })
    rows.sort(key=lambda item: item["minimum_distance_A"])
    summary = {
        "status": "CONTACTS_DETECTED_NOT_VALIDATED",
        "method": "Heavy-atom Euclidean distance screening on the selected rigid-receptor Vina pose.",
        "cutoff_A": cutoff,
        "ligand_docked_atom_count": len(ligand),
        "contact_residue_count": len(rows),
        "residues": rows[:40],
        "limitations": [
            "Polar candidates are not hydrogen-bond assignments because donor state and D-H-A angle are not evaluated.",
            "No water bridges, pi stacking geometry, metal coordination, entropy, dynamics, or free-energy calculation is included.",
            "Very close contacts must be inspected for steric clashes and preparation errors.",
        ],
    }
    with open(os.path.join(out_dir, "interactions.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    out = dict(d)
    out["interactions_summary"] = summary
    return [{"json": out}]
