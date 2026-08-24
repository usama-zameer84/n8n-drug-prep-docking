"""n8n Code node: Build Scientific Report Package.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import base64
    import subprocess
    import time
    import json
    import os
    import shutil

    d = _items[0]["json"]
    r = d["report_data"]
    out_dir = d["output_dir"]

    esc = lambda value: (str(value).replace("&", "&amp;").replace("<", "&lt;")
                         .replace(">", "&gt;").replace('"', "&quot;"))
    json_pre = lambda value: "<pre>" + esc(json.dumps(value, indent=2)) + "</pre>"
    table = lambda rows: "<table>" + "".join("<tr><th>" + esc(key) + "</th><td>" + esc(value) + "</td></tr>" for key, value in rows) + "</table>"

    style = """
    :root{color-scheme:light;font-family:Inter,system-ui,sans-serif;color:#172033;background:#f4f7fb}
    body{max-width:1120px;margin:0 auto;padding:28px}header{background:#102a43;color:white;padding:24px;border-radius:14px}
    nav{margin:18px 0;padding:14px;background:white;border-radius:10px}nav a{margin-right:14px}
    section,.card{background:white;margin:16px 0;padding:20px;border-radius:12px;box-shadow:0 2px 10px #183b5620}
    table{border-collapse:collapse;width:100%}th,td{padding:9px;border-bottom:1px solid #d9e2ec;text-align:left;vertical-align:top}
    th{width:34%;color:#334e68}code,pre{background:#eef2f7;padding:3px 6px;border-radius:5px}pre{white-space:pre-wrap;overflow-wrap:anywhere;padding:14px}
    .warn{border-left:5px solid #d97706;background:#fff7ed;padding:12px}.bad{border-left:5px solid #b91c1c;background:#fef2f2;padding:12px}
    .ok{border-left:5px solid #15803d;background:#f0fdf4;padding:12px}.muted{color:#627d98}.score{font-size:2rem;font-weight:700}
    ul{line-height:1.55}svg{max-width:100%;height:auto}
    """

    names = [
        ("index.html", "Index"),
        ("01_input_provenance.html", "Input & provenance"),
        ("02_ligand_preparation.html", "Ligand preparation"),
        ("03_receptor_preparation.html", "Receptor preparation"),
        ("03a_receptor_selection.html", "Receptor selection findings"),
        ("04_docking_results.html", "Docking results"),
        ("05_distance_contacts.html", "Distance contacts"),
        ("06_methods_qc.html", "Methods & QC"),
        ("07_md_handoff_readiness.html", "MD handoff readiness"),
        ("08_raw_data.html", "Raw data"),
        ("09_visualization.html", "3D visualization"),
    ]
    nav = "<nav>" + "".join("<a href='" + file_name + "'>" + label + "</a>" for file_name, label in names) + "</nav>"

    page = lambda title, body: "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width'><title>" + esc(title) + "</title><style>" + style + "</style></head><body><header><h1>" + esc(title) + "</h1><div>Run " + esc(r["run"]["run_id"]) + "</div></header>" + nav + body + "</body></html>"

    top = r["docking_results"]["top_affinity_kcal_mol"]
    flags = r["run"]["qc_flags"]
    index_body = "<section><h2>Outcome</h2><div class='score'>" + esc(top) + " kcal/mol</div><p>Best AutoDock Vina score for the selected replica. This is a ranking signal, not measured affinity and not MD free energy.</p></section>"
    index_body += "<section><h2>Identity</h2>" + table([
        ("Ligand", r["ligand"]["canonical_isomeric_smiles"]),
        ("InChIKey", r["ligand"]["inchi_key"]),
        ("Receptor", r["receptor"]["source"]["pdb_id"]),
        ("Receptor selection mode", (r.get("receptor_selection") or {}).get("mode")),
        ("Selected target", ((r.get("receptor_selection") or {}).get("selected") or {}).get("target_name")),
        ("Selected chains", ", ".join(r["receptor"]["selected_chains"])),
        ("Run ID", r["run"]["run_id"]),
    ]) + "</section>"
    index_body += "<section><h2>QC flags</h2>" + ("<div class='ok'>No workflow QC flags.</div>" if not flags else "<div class='warn'><ul>" + "".join("<li>" + esc(flag) + "</li>" for flag in flags) + "</ul></div>") + "</section>"
    index_body += "<section><h2>What this run establishes</h2><p>" + esc(r["scientific_scope"]["completed"]) + "</p><h3>What it does not establish</h3><p>" + esc(r["scientific_scope"]["not_completed"]) + "</p></section>"

    input_body = "<section><h2>Google Drive source</h2>" + table([
        ("Source file ID", r["run"]["source_drive_file_id"]),
        ("Source file name", r["run"]["source_drive_file_name"]),
        ("Input profile", r["run"]["input_profile"]),
        ("Run ID", r["run"]["run_id"]),
    ]) + "</section><section><h2>Ligand/receptor request</h2>" + table([
        ("Input SMILES", r["ligand"]["input_smiles"]),
        ("Canonical isomeric SMILES", r["ligand"]["canonical_isomeric_smiles"]),
        ("Receptor selection mode", (r.get("receptor_selection") or {}).get("mode")),
        ("PDB ID", r["receptor"]["source"]["pdb_id"]),
        ("RCSB URL", r["receptor"]["source"]["url"]),
        ("Source mmCIF SHA-256", r["receptor"]["source"]["sha256"]),
    ]) + "</section>"

    lig_body = "<section><h2>Deterministic 3D preparation</h2>" + table([
        ("Method", r["ligand"]["conformer_generation"]["method"]),
        ("Random seed", r["ligand"]["conformer_generation"]["random_seed"]),
        ("Force field", r["ligand"]["conformer_generation"]["force_field"]),
        ("Optimization converged", r["ligand"]["conformer_generation"]["optimization_converged"]),
        ("PDBQT method", r["ligand"]["preparation"]["method"]),
        ("Meeko version", r["ligand"]["preparation"]["meeko_version"]),
        ("SMILES mapping retained", r["ligand"]["preparation"]["smiles_mapping_present"]),
    ]) + "</section><section><h2>Descriptors</h2>" + table(list(r["ligand"]["properties"]["descriptors"].items())) + "<p class='muted'>" + esc(r["ligand"]["properties"]["interpretation_limit"]) + "</p></section>"
    try:
        with open(os.path.join(out_dir, "ligand_2d.svg")) as fh:
            lig_body += "<section><h2>2D structure</h2>" + fh.read() + "</section>"
    except OSError:
        pass

    rec = r["receptor"]
    rec_body = "<section><h2>Source and selection</h2>" + table([
        ("Source format", rec["source"]["format"]),
        ("Available chains", ", ".join(rec["available_chains"])),
        ("Selected chains", ", ".join(rec["selected_chains"])),
        ("Heterogen policy", rec["heterogen_policy"]),
        ("Waters detected", rec["waters_detected"]),
        ("Non-water heterogens detected", len(rec["non_water_heterogens_detected"])),
    ]) + "</section><section><h2>Repair and preparation QC</h2>" + table([
        ("Missing residue segments detected", len(rec["missing_residues_detected"])),
        ("Missing residues built", rec["missing_residues_built"]),
        ("Missing atoms added", rec["missing_atoms_added"]),
        ("Missing terminal atoms added", rec["missing_terminal_atoms_added"]),
        ("Prepared residues", rec["prepared_residue_count"]),
        ("Prepared atoms", rec["prepared_atom_count"]),
        ("Hydrogen pH", rec["ph"]),
        ("PDBQT method", rec["receptor_pdbqt_method"]),
        ("Charge model", rec["charge_model"]),
    ]) + "</section><section><h2>Binding-site reference</h2>" + json_pre(rec.get("binding_site_reference")) + "</section><section><h2>Detected non-water heterogens</h2>" + json_pre(rec["non_water_heterogens_detected"]) + "</section>"

    selection = r.get("receptor_selection") or {}
    selected_receptor = selection.get("selected") or {}
    candidate_targets = selection.get("candidate_targets") or []
    target_rows = ""
    for index, candidate in enumerate(candidate_targets, 1):
        best_structure = candidate.get("best_structure") or {}
        target_rows += "<tr><td>" + esc(index) + "</td><td>" + esc(candidate.get("target_name")) + "</td><td>" + esc(candidate.get("target_chembl_id")) + "</td><td>" + esc(candidate.get("uniprot_accession")) + "</td><td>" + esc(candidate.get("max_similarity_percent")) + "</td><td>" + esc(candidate.get("evidence_records")) + "</td><td>" + esc(candidate.get("max_pchembl")) + "</td><td>" + esc(best_structure.get("pdb_id")) + "</td><td>" + esc(best_structure.get("resolution_A")) + "</td><td>" + esc(candidate.get("evidence_score")) + "</td></tr>"
    if not target_rows:
        target_rows = "<tr><td colspan='10'>No inferred candidates: a user-supplied receptor was used.</td></tr>"
    selection_body = "<section><h2>Selection outcome</h2>" + table([
        ("Mode", selection.get("mode")),
        ("Status", selection.get("status")),
        ("Method", selection.get("method")),
        ("Selected target", selected_receptor.get("target_name")),
        ("Selected target ChEMBL ID", selected_receptor.get("target_chembl_id")),
        ("Selected UniProt accession", selected_receptor.get("uniprot_accession")),
        ("Selected organism", selected_receptor.get("organism")),
        ("Selected PDB ID", selected_receptor.get("pdb_id")),
        ("Selected chain(s)", ", ".join(selected_receptor.get("chain_ids") or [])),
        ("Experimental method", selected_receptor.get("experimental_method")),
        ("Resolution (A)", selected_receptor.get("resolution_A")),
        ("Combined selection score", selected_receptor.get("selection_score")),
    ]) + "</section>"
    selection_body += "<section><h2>Ranked target findings</h2><table><tr><th>Rank</th><th>Target</th><th>ChEMBL</th><th>UniProt</th><th>Max similarity %</th><th>Evidence records</th><th>Max pChEMBL</th><th>Best PDB</th><th>Resolution A</th><th>Evidence score</th></tr>" + target_rows + "</table></section>"
    selection_body += "<section><h2>Why this receptor was used</h2><ol>" + "".join("<li>" + esc(item) + "</li>" for item in selection.get("rationale", [])) + "</ol><h3>How it entered docking</h3><p>The selected PDB ID and UniProt/SIFTS-mapped chain were passed into PDBFixer and official Meeko receptor preparation. When no grid centre was supplied, an eligible co-crystallized component was ranked by RDKit fingerprint similarity to the query ligand and its coordinate centroid was used as the Vina box centre.</p></section>"
    selection_body += "<section><h2>Preparation compatibility preflight</h2>" + json_pre(selection.get("preparation_preflight", [])) + "</section><section><h2>Evidence query</h2>" + json_pre(selection.get("query", {})) + "<h3>Similar ChEMBL molecules</h3>" + json_pre(selection.get("candidate_molecules", [])) + "<h3>Data sources</h3>" + json_pre(selection.get("data_sources", {})) + "</section>"
    selection_body += "<section><h2>Mandatory limitations</h2><div class='warn'><ul>" + "".join("<li>" + esc(item) + "</li>" for item in selection.get("limitations", [])) + "</ul></div></section>"

    cfg = r["docking_configuration"]
    poses = r["docking_results"]["poses"]
    pose_rows = "".join("<tr><td>" + esc(p["rank"]) + "</td><td>" + esc(p["affinity_kcal_mol"]) + "</td><td>" + esc(p["rmsd_lb_A"]) + "</td><td>" + esc(p["rmsd_ub_A"]) + "</td></tr>" for p in poses)
    dock_body = "<section><h2>Configuration</h2>" + table([
        ("Vina version", cfg["vina_version"]),
        ("Grid center (A)", cfg["grid_center_A"]),
        ("Grid size (A)", cfg["grid_size_A"]),
        ("Auto-centered", cfg["grid_auto_centered"]),
        ("Grid center source", cfg.get("grid_center_source")),
        ("Binding-site reference", cfg.get("binding_site_reference")),
        ("Exhaustiveness", cfg["exhaustiveness"]),
        ("Requested modes", cfg["num_modes"]),
        ("Energy range (kcal/mol)", cfg["energy_range_kcal_mol"]),
        ("CPU", cfg["cpu"]),
        ("Replica seeds", cfg["replica_seeds"]),
        ("Timeout per replica (s)", cfg["timeout_seconds"]),
    ]) + "</section><section><h2>Ranked poses: selected replica</h2><table><tr><th>Rank</th><th>Vina score (kcal/mol)</th><th>RMSD lower (A)</th><th>RMSD upper (A)</th></tr>" + pose_rows + "</table><p class='warn'>" + esc(r["docking_results"]["score_semantics"]) + "</p></section><section><h2>Replica consistency</h2>" + table([
        ("Top scores", r["docking_results"]["replica_top_scores_kcal_mol"]),
        ("Mean", r["docking_results"]["replica_top_score_mean"]),
        ("Std dev", r["docking_results"]["replica_top_score_stddev"]),
        ("Range", r["docking_results"]["replica_top_score_range"]),
    ]) + "<p>" + esc(r["docking_results"]["convergence_note"]) + "</p></section>"

    contacts = r["distance_contacts"]
    contact_rows = "".join("<tr><td>" + esc(item["residue"]) + "</td><td>" + esc(item["minimum_distance_A"]) + "</td><td>" + esc(item["contact_count"]) + "</td><td>" + esc(item["contacts"][0]["category"]) + "</td></tr>" for item in contacts["residues"])
    contact_body = "<section><h2>Transparent geometric screen</h2>" + table([
        ("Status", contacts["status"]),
        ("Method", contacts["method"]),
        ("Cutoff (A)", contacts["cutoff_A"]),
        ("Contact residues", contacts["contact_residue_count"]),
    ]) + "</section><section><table><tr><th>Residue</th><th>Minimum distance (A)</th><th>Contacts</th><th>Closest category</th></tr>" + contact_rows + "</table></section><section><h2>Limitations</h2><ul>" + "".join("<li>" + esc(item) + "</li>" for item in contacts["limitations"]) + "</ul></section>"

    qc_body = "<section><h2>Recorded software</h2>" + table(list(rec["versions"].items()) + [("vina", cfg["vina_version"])]) + "</section><section><h2>Reproducibility parameters</h2>" + json_pre(cfg) + "</section><section><h2>QC flags</h2>" + json_pre(flags) + "</section><section><h2>Interpretation rules</h2><ul><li>Review the grid in 3D; receptor-centroid auto-placement is only a smoke-test default.</li><li>Review removed heterogens, metals and waters against the biological question.</li><li>Inspect protonation, tautomers, stereochemistry, missing residues and close contacts.</li><li>Use independent replicas/ensembles and experimental evidence for decisions.</li></ul></section>"

    manifest = r["md_handoff"]["manifest"]
    md_body = "<section><h2>Status</h2><div class='bad'>STRUCTURAL HANDOFF ONLY — NOT TOPOLOGY-READY</div>" + table([
        ("Chemistry QC passed", r["md_handoff"]["chemistry_qc_passed"]),
        ("Chemistry authority", manifest["chemistry_authority"]),
        ("Coordinate complex", manifest["coordinate_complex"]),
        ("Atom mapping", manifest["atom_mapping"]),
        ("Input heavy atoms", manifest["input_atom_counts"]["heavy_atoms"]),
        ("Exported heavy atoms", manifest["exported_atom_counts"]["heavy_atoms"]),
        ("Complex atoms", manifest["complex_atom_counts"]["total"]),
    ]) + "</section><section><h2>Required before real MD</h2><ol><li>Inspect pose and receptor chemistry.</li><li>Select protein/water/ion and ligand force fields.</li><li>Generate and validate ligand topology/charges.</li><li>Reconcile atom names using atom_mapping.json.</li><li>Define cofactors, metals, waters, termini and disulfides.</li><li>Build topology, periodic box, solvent and ions.</li><li>Minimize, equilibrate NVT/NPT, then run replicated production MD.</li><li>Analyze RMSD/RMSF, contacts, convergence and uncertainty.</li></ol></section>"

    raw_body = "<section><h2>Machine-readable report</h2>" + json_pre(r) + "</section>"

    with open(os.path.join(out_dir, "prepared_receptor.pdb")) as fh:
        receptor_text = fh.read()
    with open(os.path.join(out_dir, "docked.pdbqt")) as fh:
        ligand_text = fh.read()
    visual_body = """<section><h2>Selected pose viewer</h2><div class='warn'>Requires internet access to load the 3Dmol.js viewer. Always inspect the SDF/PDB files in a validated molecular viewer.</div><div id='viewer' style='height:680px;width:100%;position:relative'></div></section>
    <script src='https://3Dmol.org/build/3Dmol-min.js'></script><script>
    const receptor = """ + json.dumps(receptor_text) + """;
    const ligand = """ + json.dumps(ligand_text) + """;
    const viewer = $3Dmol.createViewer('viewer',{backgroundColor:'white'});
    viewer.addModel(receptor,'pdb'); viewer.setStyle({model:0},{cartoon:{color:'spectrum'}});
    viewer.addModel(ligand,'pdbqt'); viewer.setStyle({model:1},{stick:{colorscheme:'greenCarbon'}});
    viewer.zoomTo(); viewer.render();
    </script>"""

    overlay_profile = {
        "pdb_id": rec["source"]["pdb_id"],
        "source_sha256": rec["source"]["sha256"],
        "selected_chains": rec["selected_chains"],
        "heterogen_policy": rec["heterogen_policy"],
        "ph": rec["ph"],
    }
    overlay_profile_text = json.dumps(overlay_profile, sort_keys=True, separators=(",", ":"))
    overlay_hash = subprocess.run(["sha256sum"], input=overlay_profile_text, text=True, capture_output=True, check=True).stdout.split()[0][:16]
    chain_label = "-".join(rec["selected_chains"]) or "all"
    overlay_group_id = rec["source"]["pdb_id"] + "_" + chain_label + "_" + overlay_hash
    overlay_root = os.path.join("/md_project/data/receptor_overlays", overlay_group_id)
    pose_root = os.path.join(overlay_root, "poses")
    os.makedirs(pose_root, exist_ok=True)
    lock_path = os.path.join(overlay_root, ".registry.lock")
    registry_path = os.path.join(overlay_root, "overlay_manifest.json")
    reference_path = os.path.join(overlay_root, "receptor_reference.pdb")
    current_pose_source = os.path.join(out_dir, "best_pose_ligand.sdf")
    current_pose_name = d["run_id"] + ".sdf"
    current_pose_path = os.path.join(pose_root, current_pose_name)
    with open(lock_path, "a+") as lock_fh:
        lock_dir = lock_path + ".d"
        lock_acquired = False
        for _ in range(300):
            try:
                os.mkdir(lock_dir)
                lock_acquired = True
                break
            except FileExistsError:
                try:
                    if time.time() - os.path.getmtime(lock_dir) > 600:
                        shutil.rmtree(lock_dir, ignore_errors=True)
                        continue
                except OSError:
                    pass
                time.sleep(0.1)
        if not lock_acquired:
            raise RuntimeError("Timed out acquiring receptor-overlay registry lock")
        try:
            try:
                with open(registry_path) as fh:
                    overlay_registry = json.load(fh)
            except (OSError, ValueError):
                overlay_registry = {"schema_version": "1.0", "group_id": overlay_group_id, "receptor_profile": overlay_profile, "entries": []}
            if not os.path.exists(reference_path):
                shutil.copyfile(os.path.join(out_dir, "prepared_receptor.pdb"), reference_path)
            shutil.copyfile(current_pose_source, current_pose_path)
            entries = [entry for entry in overlay_registry.get("entries", []) if entry.get("run_id") != d["run_id"] and os.path.exists(os.path.join(pose_root, entry.get("pose_file", "")))]
            selected_target = (r.get("receptor_selection") or {}).get("selected") or {}
            entries.append({
                "run_id": d["run_id"],
                "run_name": d["run_name"],
                "ligand_name": d["ligand_name"],
                "canonical_smiles": r["ligand"]["canonical_isomeric_smiles"],
                "inchi_key": r["ligand"]["inchi_key"],
                "top_vina_score_kcal_mol": r["docking_results"]["top_affinity_kcal_mol"],
                "grid_center_A": r["docking_configuration"]["grid_center_A"],
                "grid_size_A": r["docking_configuration"]["grid_size_A"],
                "grid_center_source": r["docking_configuration"].get("grid_center_source"),
                "target_name": selected_target.get("target_name"),
                "target_chembl_id": selected_target.get("target_chembl_id"),
                "pose_file": current_pose_name,
            })
            entries.sort(key=lambda entry: (entry.get("ligand_name") or "", entry.get("run_id") or ""))
            overlay_registry["entries"] = entries
            overlay_registry["ligand_count"] = len(entries)
            overlay_registry["coordinate_policy"] = "All poses share the exact prepared receptor profile identified by source SHA-256, chains, heterogen policy and pH; no structural transform is required."
            temp_registry = registry_path + ".tmp"
            with open(temp_registry, "w") as fh:
                json.dump(overlay_registry, fh, indent=2)
            os.replace(temp_registry, registry_path)
        finally:
            shutil.rmtree(lock_dir, ignore_errors=True)

    with open(reference_path) as fh:
        overlay_receptor_text = fh.read()
    overlay_records = []
    combined_sdf_parts = []
    public_entries = []
    palette = ["#e11d48", "#2563eb", "#16a34a", "#9333ea", "#ea580c", "#0891b2", "#ca8a04", "#db2777", "#4f46e5", "#059669", "#7c3aed", "#dc2626"]
    for index, entry in enumerate(overlay_registry["entries"]):
        pose_path = os.path.join(pose_root, entry["pose_file"])
        try:
            with open(pose_path) as fh:
                sdf_text = fh.read()
        except OSError:
            continue
        if "$$$$" not in sdf_text:
            sdf_text = sdf_text.rstrip() + "\n$$$$\n"
        combined_sdf_parts.append(sdf_text.rstrip() + "\n")
        public_entry = dict(entry)
        public_entry["color"] = palette[index % len(palette)]
        public_entries.append(public_entry)
        overlay_records.append({"metadata": public_entry, "sdf": sdf_text})
    combined_sdf = "\n".join(combined_sdf_parts)
    overlay_public = {
        "schema_version": overlay_registry["schema_version"],
        "group_id": overlay_group_id,
        "receptor_profile": overlay_profile,
        "ligand_count": len(public_entries),
        "coordinate_policy": overlay_registry["coordinate_policy"],
        "entries": public_entries,
    }
    control_rows = "".join(
        "<tr><td><input type='checkbox' checked onchange='togglePose(" + str(index) + ",this.checked)'></td><td><span style='display:inline-block;width:13px;height:13px;border-radius:50%;background:" + esc(entry["color"]) + "'></span></td><td>" + esc(entry.get("ligand_name")) + "</td><td>" + esc(entry.get("run_id")) + "</td><td>" + esc(entry.get("top_vina_score_kcal_mol")) + "</td><td>" + esc(entry.get("grid_center_source")) + "</td></tr>"
        for index, entry in enumerate(public_entries)
    )
    visual_body = "<section><h2>Superimposed best poses</h2><div class='warn'>This is a receptor-profile-matched overlay of predicted docking poses, not proof of binding or covalent attachment. It requires internet access to load 3Dmol.js.</div>" + table([
        ("Overlay group", overlay_group_id),
        ("Receptor", overlay_profile["pdb_id"]),
        ("Chains", ", ".join(overlay_profile["selected_chains"])),
        ("Ligands in this snapshot", len(public_entries)),
        ("Coordinate policy", overlay_registry["coordinate_policy"]),
    ]) + "<div style='margin:12px 0'><button onclick='setAll(true)'>Show all</button> <button onclick='setAll(false)'>Hide all</button> <button onclick='viewer.zoomTo();viewer.render()'>Fit view</button></div><div id='viewer' style='height:720px;width:100%;position:relative'></div></section>"
    visual_body += "<section><h2>Ligand controls and scores</h2><table><tr><th>Show</th><th>Color</th><th>Ligand</th><th>Run ID</th><th>Top Vina score</th><th>Grid source</th></tr>" + control_rows + "</table><p class='muted'>The overlay is a snapshot generated during this run. Later runs with the same receptor profile will appear in later snapshots.</p></section>"
    visual_body += "<script src='https://3Dmol.org/build/3Dmol-min.js'></script><script>const receptor=" + json.dumps(overlay_receptor_text) + ";const poseData=" + json.dumps(overlay_records) + ";const viewer=$3Dmol.createViewer('viewer',{backgroundColor:'white'});viewer.addModel(receptor,'pdb');viewer.setStyle({model:0},{cartoon:{color:'spectrum'}});const ligandModels=[];const ligandStyles=[];poseData.forEach((pose,i)=>{const model=viewer.addModel(pose.sdf,'sdf');const style={stick:{color:pose.metadata.color,radius:0.22},sphere:{color:pose.metadata.color,scale:0.22}};model.setStyle({},style);ligandModels.push(model);ligandStyles.push(style);});function togglePose(i,on){ligandModels[i].setStyle({},on?ligandStyles[i]:{});viewer.render();}function setAll(on){document.querySelectorAll(\"input[type='checkbox']\").forEach((box,i)=>{box.checked=on;togglePose(i,on);});}viewer.zoomTo();viewer.render();</script>"

    documents = {
        "index.html": page("Docking run report", index_body),
        "01_input_provenance.html": page("Input and provenance", input_body),
        "02_ligand_preparation.html": page("Ligand preparation", lig_body),
        "03_receptor_preparation.html": page("Receptor preparation", rec_body),
        "03a_receptor_selection.html": page("Automatic receptor-selection findings", selection_body),
        "04_docking_results.html": page("Docking configuration and ranked poses", dock_body),
        "05_distance_contacts.html": page("Distance-contact screen", contact_body),
        "06_methods_qc.html": page("Methods, reproducibility and QC", qc_body),
        "07_md_handoff_readiness.html": page("MD handoff readiness", md_body),
        "08_raw_data.html": page("Raw machine-readable data", raw_body),
        "09_visualization.html": page("3D visualization", visual_body),
    }
    documents["receptor_overlay_reference.pdb"] = overlay_receptor_text
    documents["all_best_poses.sdf"] = combined_sdf
    documents["ligand_overlay.json"] = json.dumps(overlay_public, indent=2)
    run_summary = {
        "run_id": r["run"]["run_id"],
        "run_name": d["run_name"],
        "status": "DOCKING_COMPLETED_STRUCTURAL_MD_HANDOFF_CREATED",
        "top_vina_score_kcal_mol": top,
        "selected_replica": r["docking_results"]["selected_replica"],
        "qc_flags": flags,
        "md_status": r["md_handoff"]["status"],
        "receptor_selection": {"mode": selection.get("mode"), "status": selection.get("status"), "selected": selected_receptor},
        "overlay_group": {"group_id": overlay_group_id, "ligand_count": len(public_entries), "artifacts": ["09_visualization.html", "receptor_overlay_reference.pdb", "all_best_poses.sdf", "ligand_overlay.json"]},
        "report_files": sorted([name for name in documents] + ["run_summary.json"]),
    }
    documents["run_summary.json"] = json.dumps(run_summary, indent=2)

    reports_b64 = {name: base64.b64encode(content.encode("utf-8")).decode("ascii") for name, content in documents.items()}
    out = dict(d)
    out["reports_b64"] = reports_b64
    out["report_summary"] = run_summary
    return [{"json": out}]
