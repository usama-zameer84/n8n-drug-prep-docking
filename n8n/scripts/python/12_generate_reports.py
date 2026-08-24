"""n8n Code node: Generate Reports.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import json
    import os

    d = _items[0]["json"]
    record = {
        "run": {
            "run_id": d["run_id"],
            "run_name": d["run_name"],
            "input_profile": d["input_profile"],
            "source_drive_file_id": d.get("source_file_id"),
            "source_drive_file_name": d.get("source_file_name"),
            "qc_flags": d["qc_flags"],
        },
        "ligand": {
            "input_smiles": d["smiles"],
            "canonical_isomeric_smiles": d["canonical_smiles"],
            "inchi_key": d["inchi_key"],
            "atom_counts": {"all": d["num_atoms"], "heavy": d["num_heavy_atoms"]},
            "conformer_generation": d["conformer_generation"],
            "preparation": d["ligand_prep"],
            "properties": d["ligand_analysis"],
        },
        "receptor": d["protein_stats"],
        "receptor_selection": d.get("receptor_selection", {}),
        "docking_configuration": {
            "vina_version": d["vina_version"],
            "grid_center_A": d["grid_center"],
            "grid_size_A": d["grid_size"],
            "grid_auto_centered": d["grid_auto_centered"],
            "grid_center_source": d.get("grid_center_source"),
            "binding_site_reference": d["protein_stats"].get("binding_site_reference"),
            "receptor_bounds_A": d["receptor_bounds"],
            "exhaustiveness": d["exhaustiveness"],
            "num_modes": d["num_modes"],
            "energy_range_kcal_mol": d["energy_range"],
            "cpu": d["cpu"],
            "replica_seeds": [item["seed"] for item in d["replica_runs"]],
            "timeout_seconds": d["timeout_seconds"],
        },
        "docking_results": d["docking_results"],
        "distance_contacts": d["interactions_summary"],
        "md_handoff": d["md_handoff"],
        "scientific_scope": {
            "completed": "Ligand/receptor preparation, rigid-receptor AutoDock Vina docking, ranked pose export, distance-contact screening, structural MD handoff.",
            "not_completed": "Experimental validation, flexible-receptor/ensemble docking, topology generation, solvation, ions, minimization, NVT/NPT, production MD, convergence, MM/PBSA, alchemical free energy.",
        },
    }
    record_path = os.path.join(d["output_dir"], "report_data.json")
    with open(record_path, "w") as fh:
        json.dump(record, fh, indent=2)
    out = dict(d)
    out["report_data"] = record
    out["report_data_path"] = record_path
    return [{"json": out}]
