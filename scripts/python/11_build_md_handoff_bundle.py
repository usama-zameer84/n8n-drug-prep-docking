"""n8n Code node: Build MD Handoff Bundle.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import base64
    import json
    import os
    import subprocess
    import sys
    from rdkit import Chem

    d = _items[0]["json"]
    out_dir = d["output_dir"]
    docked_pdbqt = os.path.join(out_dir, "docked.pdbqt")
    docked_sdf = os.path.join(out_dir, "docked_poses.sdf")
    export_cmd = [sys.executable, "-m", "meeko.cli.mk_export", docked_pdbqt, "--write_sdf", docked_sdf]
    try:
        proc = subprocess.run(export_cmd, capture_output=True, text=True, timeout=180)
    except subprocess.TimeoutExpired:
        raise ValueError("Official Meeko pose export exceeded 180 seconds")
    export_log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    with open(os.path.join(out_dir, "meeko_export.log"), "w") as fh:
        fh.write(export_log)
    if proc.returncode != 0 or not os.path.exists(docked_sdf):
        raise ValueError("Official Meeko mk_export failed: " + export_log[-1600:])

    poses = [mol for mol in Chem.SDMolSupplier(docked_sdf, removeHs=False) if mol is not None]
    if not poses:
        raise ValueError("Meeko exported no chemically reconstructed ligand poses")
    best = poses[0]
    docked_heavy_atoms = best.GetNumHeavyAtoms()
    if docked_heavy_atoms != int(d["num_heavy_atoms"]):
        raise ValueError("Chemistry-preserving export QC failed: input heavy atoms " + str(d["num_heavy_atoms"]) +
                         " != docked heavy atoms " + str(docked_heavy_atoms))

    atom_mapping = []
    element_counts = {}
    for index, atom in enumerate(best.GetAtoms()):
        symbol = atom.GetSymbol()
        element_counts[symbol] = element_counts.get(symbol, 0) + 1
        atom_name = (symbol + str(element_counts[symbol]))[:4]
        info = Chem.AtomPDBResidueInfo()
        info.SetName(atom_name.rjust(4))
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId("Z")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
        atom_mapping.append({
            "docked_sdf_index_0_based": index,
            "pdb_atom_name": atom_name,
            "element": symbol,
            "atomic_number": atom.GetAtomicNum(),
            "formal_charge": atom.GetFormalCharge(),
        })

    best_sdf = os.path.join(out_dir, "best_pose_ligand.sdf")
    writer = Chem.SDWriter(best_sdf)
    writer.write(best)
    writer.close()
    best_pdb = os.path.join(out_dir, "best_pose_ligand.pdb")
    Chem.MolToPDBFile(best, best_pdb)

    mapping_record = {
        "schema": "docked SDF atom index to unique ligand PDB atom name",
        "atom_mapping": atom_mapping,
        "canonical_isomeric_smiles": d["canonical_smiles"],
        "inchi_key": d["inchi_key"],
        "note": "The SDF is the chemistry-authoritative ligand structure. PDB atom names are a stable coordinate handoff for topology generation.",
    }
    mapping_path = os.path.join(out_dir, "atom_mapping.json")
    with open(mapping_path, "w") as fh:
        json.dump(mapping_record, fh, indent=2)

    prepared_receptor = os.path.join(out_dir, "prepared_receptor.pdb")
    complex_path = os.path.join(out_dir, "complex_best_pose.pdb")
    serial = 1
    receptor_atom_count = 0
    ligand_atom_count = 0
    with open(complex_path, "w") as out_fh:
        with open(prepared_receptor) as rec_fh:
            for line in rec_fh:
                if line.startswith(("ATOM", "HETATM")):
                    out_fh.write(line[:6] + str(serial).rjust(5) + line[11:].rstrip("\n") + "\n")
                    serial += 1
                    receptor_atom_count += 1
        out_fh.write("TER\n")
        with open(best_pdb) as lig_fh:
            for line in lig_fh:
                if line.startswith(("ATOM", "HETATM")):
                    out_fh.write("HETATM" + str(serial).rjust(5) + line[11:].rstrip("\n") + "\n")
                    serial += 1
                    ligand_atom_count += 1
        out_fh.write("END\n")

    if ligand_atom_count != best.GetNumAtoms():
        raise ValueError("Complex assembly lost ligand atoms")
    if serial - 1 != receptor_atom_count + ligand_atom_count:
        raise ValueError("Complex assembly serial/count invariant failed")

    readme = """# MD structural handoff

    Status: not topology-ready.

    This folder contains a chemistry-preserving, coordinate-level handoff from rigid-receptor AutoDock Vina docking. Use docked_poses.sdf or best_pose_ligand.sdf as the ligand chemistry authority, prepared_receptor.pdb as the prepared protein coordinates, and complex_best_pose.pdb only as a combined starting coordinate model.

    Before any real GROMACS MD, independently inspect the binding pose, choose and document protein/ligand force fields, generate ligand parameters and charges, reconcile protonation and atom names with atom_mapping.json, define termini/disulfides/cofactors/metals/waters, build the topology, solvate, add ions, minimize, equilibrate (NVT/NPT), and run replicated production MD with convergence analysis.

    Docking scores are ranking signals. They are not experimental affinities and not MD free energies.
    """
    readme_path = os.path.join(out_dir, "README_MD_HANDOFF.md")
    with open(readme_path, "w") as fh:
        fh.write(readme)

    provenance = {
        "run_id": d["run_id"],
        "source_drive": {
            "file_id": d.get("source_file_id"),
            "file_name": d.get("source_file_name"),
            "mime_type": d.get("source_mime_type"),
            "modified_time": d.get("source_modified_time"),
            "md5": d.get("source_md5"),
        },
        "ligand": {
            "input_smiles": d["smiles"],
            "canonical_isomeric_smiles": d["canonical_smiles"],
            "inchi_key": d["inchi_key"],
            "conformer_generation": d["conformer_generation"],
            "preparation": d["ligand_prep"],
        },
        "receptor": d["protein_stats"],
        "receptor_selection": d.get("receptor_selection", {}),
        "docking": {
            "vina_version": d["vina_version"],
            "grid_center_A": d["grid_center"],
            "grid_size_A": d["grid_size"],
            "grid_auto_centered": d["grid_auto_centered"],
            "grid_center_source": d.get("grid_center_source"),
            "exhaustiveness": d["exhaustiveness"],
            "num_modes": d["num_modes"],
            "energy_range_kcal_mol": d["energy_range"],
            "cpu": d["cpu"],
            "replica_seeds": [item["seed"] for item in d["replica_runs"]],
            "selected_replica": d["selected_replica"],
            "results": d["docking_results"],
        },
        "contacts": d["interactions_summary"],
        "qc_flags": d["qc_flags"],
        "md_status": "STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY",
    }
    provenance_path = os.path.join(out_dir, "provenance.json")
    with open(provenance_path, "w") as fh:
        json.dump(provenance, fh, indent=2)

    files = {
        "receptor_source.cif": d["receptor_source_path"],
        "prepared_receptor.pdb": prepared_receptor,
        "receptor.pdbqt": os.path.join(out_dir, "receptor.pdbqt"),
        "receptor.json": os.path.join(out_dir, "receptor.json"),
        "receptor_prep.log": os.path.join(out_dir, "receptor_prep.log"),
        "ligand_input.sdf": os.path.join(out_dir, "ligand_input.sdf"),
        "ligand.pdbqt": os.path.join(out_dir, "ligand.pdbqt"),
        "ligand_2d.svg": os.path.join(out_dir, "ligand_2d.svg"),
        "docked.pdbqt": docked_pdbqt,
        "docked_poses.sdf": docked_sdf,
        "best_pose_ligand.sdf": best_sdf,
        "best_pose_ligand.pdb": best_pdb,
        "complex_best_pose.pdb": complex_path,
        "atom_mapping.json": mapping_path,
        "docking.log": os.path.join(out_dir, "docking.log"),
        "docking_results.json": os.path.join(out_dir, "docking_results.json"),
        "interactions.json": os.path.join(out_dir, "interactions.json"),
        "protein_stats.json": os.path.join(out_dir, "protein_stats.json"),
        "meeko_export.log": os.path.join(out_dir, "meeko_export.log"),
        "provenance.json": provenance_path,
        "README_MD_HANDOFF.md": readme_path,
    }
    for item in d["replica_runs"]:
        files["docked_replica_" + str(item["replica"]) + ".pdbqt"] = item["pose_path"]
        files["vina_replica_" + str(item["replica"]) + ".log"] = item["log_path"]

    def sha256(path):
        result = subprocess.run(["sha256sum", path], capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            raise ValueError("sha256sum failed for " + path)
        return result.stdout.split()[0]

    manifest = {
        "schema_version": "2.0",
        "run_id": d["run_id"],
        "status": "not topology-ready",
        "chemistry_authority": "best_pose_ligand.sdf",
        "coordinate_complex": "complex_best_pose.pdb",
        "atom_mapping": "atom_mapping.json",
        "input_atom_counts": {"all_atoms": d["num_atoms"], "heavy_atoms": d["num_heavy_atoms"]},
        "exported_atom_counts": {"all_atoms": best.GetNumAtoms(), "heavy_atoms": docked_heavy_atoms},
        "complex_atom_counts": {"receptor": receptor_atom_count, "ligand": ligand_atom_count,
                                "total": receptor_atom_count + ligand_atom_count},
        "files": {name: {"sha256": sha256(path), "bytes": os.path.getsize(path)} for name, path in files.items()},
    }
    manifest_path = os.path.join(out_dir, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    files["manifest.json"] = manifest_path

    md_files_b64 = {}
    for name, path in files.items():
        with open(path, "rb") as fh:
            md_files_b64[name] = base64.b64encode(fh.read()).decode("ascii")

    out = dict(d)
    out["md_handoff"] = {
        "status": "STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY",
        "manifest": manifest,
        "file_count": len(md_files_b64),
        "chemistry_qc_passed": True,
    }
    out["md_files_b64"] = md_files_b64
    return [{"json": out}]
