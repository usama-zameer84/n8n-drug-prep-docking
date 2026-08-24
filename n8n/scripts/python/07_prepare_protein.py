"""n8n Code node: Prepare Protein.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import json
    import os
    import subprocess
    import sys
    import urllib.request
    import urllib.parse
    import meeko
    import rdkit
    import openmm
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from openmm import unit
    from rdkit import Chem, DataStructs

    d = _items[0]["json"]
    input_dir = d["input_dir"]
    out_dir = d["output_dir"]
    source_cif = os.path.join(input_dir, d["pdb_id"] + ".cif")
    url = "https://files.rcsb.org/download/" + d["pdb_id"] + ".cif"
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            source_bytes = response.read()
    except Exception as exc:
        raise ValueError("RCSB mmCIF download failed for " + d["pdb_id"] + ": " + str(exc))
    if len(source_bytes) < 500:
        raise ValueError("RCSB returned an unexpectedly small mmCIF file")
    with open(source_cif, "wb") as fh:
        fh.write(source_bytes)

    fixer = PDBFixer(filename=source_cif)
    chains = list(fixer.topology.chains())
    available_chains = [chain.id or str(index) for index, chain in enumerate(chains)]
    protein_names = set(("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASX GLX SEC PYL DA DC DG DT A C G U".split()))
    water_names = set(("HOH WAT H2O DOD".split()))
    heterogens = []
    waters = []
    for residue in fixer.topology.residues():
        item = {"name": residue.name, "chain": residue.chain.id, "id": residue.id}
        if residue.name in water_names:
            waters.append(item)
        elif residue.name not in protein_names:
            heterogens.append(item)

    binding_site_candidates = []
    query_mol = Chem.MolFromSmiles(d["smiles"])
    query_fp = Chem.RDKFingerprint(query_mol) if query_mol is not None else None
    excluded_components = set(("HOH WAT H2O DOD NA CL K CA MG MN ZN FE CU CO NI CD HG SO4 PO4 NO3 NH4 GOL EDO PEG PG4 DMS ACT ACE".split()))
    for residue in fixer.topology.residues():
        if residue.name in protein_names or residue.name in water_names or residue.name in excluded_components:
            continue
        heavy_atoms = [atom for atom in residue.atoms() if atom.element is not None and atom.element.symbol != "H"]
        if len(heavy_atoms) < 5:
            continue
        coordinates_A = []
        for atom in heavy_atoms:
            point = fixer.positions[atom.index].value_in_unit(unit.angstrom)
            coordinates_A.append([float(point.x), float(point.y), float(point.z)])
        similarity = None
        component_smiles = None
        component_name = residue.name
        try:
            component_url = "https://data.rcsb.org/rest/v1/core/chemcomp/" + urllib.parse.quote(residue.name)
            req = urllib.request.Request(component_url, headers={"Accept": "application/json", "User-Agent": "ligand-docking-workbench/3.0 pocket-selection"})
            with urllib.request.urlopen(req, timeout=25) as response:
                component_data = json.load(response)
            descriptor = component_data.get("rcsb_chem_comp_descriptor") or {}
            component_smiles = descriptor.get("SMILES_stereo") or descriptor.get("SMILES")
            component_name = (component_data.get("chem_comp") or {}).get("name") or residue.name
            component_mol = Chem.MolFromSmiles(component_smiles) if component_smiles else None
            if query_fp is not None and component_mol is not None:
                similarity = float(DataStructs.TanimotoSimilarity(query_fp, Chem.RDKFingerprint(component_mol)))
        except Exception:
            pass
        center_A = [sum(row[index] for row in coordinates_A) / len(coordinates_A) for index in range(3)]
        binding_site_candidates.append({
            "component_id": residue.name,
            "component_name": component_name,
            "chain": residue.chain.id,
            "residue_id": residue.id,
            "heavy_atom_count": len(heavy_atoms),
            "query_similarity": round(similarity, 4) if similarity is not None else None,
            "center_A": [round(value, 3) for value in center_A],
            "component_smiles": component_smiles,
        })
    binding_site_candidates.sort(key=lambda item: (-(item["query_similarity"] if item["query_similarity"] is not None else -1.0), -item["heavy_atom_count"], item["component_id"]))
    binding_site_reference = None
    if binding_site_candidates:
        binding_site_reference = dict(binding_site_candidates[0])
        binding_site_reference["selection_method"] = "Highest RDKit fingerprint similarity between the query ligand and eligible co-crystallized non-polymer components; heavy-atom count breaks ties."
        binding_site_reference["candidate_count"] = len(binding_site_candidates)

    if d["chain_ids"]:
        missing_chains = [chain for chain in d["chain_ids"] if chain not in available_chains]
        if missing_chains:
            raise ValueError("Requested chains not present: " + ",".join(missing_chains) + "; available: " + ",".join(available_chains))
        remove_indices = [index for index, chain in enumerate(chains) if (chain.id or str(index)) not in d["chain_ids"]]
        if remove_indices:
            fixer.removeChains(remove_indices)

    fixer.findMissingResidues()
    missing_residues_detected = {str(key): list(value) for key, value in fixer.missingResidues.items()}
    if not d["add_missing_residues"]:
        fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    nonstandard = [{"name": residue.name, "chain": residue.chain.id, "id": residue.id, "replacement": replacement}
                   for residue, replacement in (fixer.nonstandardResidues or [])]
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=d["heterogen_policy"] == "keep_water")
    fixer.findMissingAtoms()
    missing_atoms_count = sum(len(atoms) for atoms in fixer.missingAtoms.values())
    missing_terminals_count = sum(len(atoms) for atoms in fixer.missingTerminals.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(float(d["ph"]))

    prepared_pdb = os.path.join(out_dir, "prepared_receptor.pdb")
    with open(prepared_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    receptor_pdbqt = os.path.join(out_dir, "receptor.pdbqt")
    receptor_json = os.path.join(out_dir, "receptor.json")
    prep_cmd = [
        sys.executable, "-m", "meeko.cli.mk_prepare_receptor",
        "--read_pdb", prepared_pdb,
        "--write_pdbqt", receptor_pdbqt,
        "--write_json", receptor_json,
        "--charge_model", "gasteiger",
        "--allow_bad_res",
    ]
    try:
        proc = subprocess.run(prep_cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        raise ValueError("Meeko receptor preparation exceeded 300 seconds")
    prep_log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    with open(os.path.join(out_dir, "receptor_prep.log"), "w") as fh:
        fh.write(prep_log)
    if proc.returncode != 0 or not os.path.exists(receptor_pdbqt) or os.path.getsize(receptor_pdbqt) == 0:
        raise ValueError("Official Meeko receptor preparation failed: " + prep_log[-1600:])

    def sha256(path):
        result = subprocess.run(["sha256sum", path], capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            raise ValueError("sha256sum failed for " + path)
        return result.stdout.split()[0]

    atom_count = 0
    residue_keys = set()
    with open(prepared_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atom_count += 1
                residue_keys.add((line[21:22], line[22:27], line[17:20]))

    stats = {
        "source": {
            "pdb_id": d["pdb_id"],
            "format": "PDBx/mmCIF",
            "url": url,
            "sha256": sha256(source_cif),
            "drive_input_file_id": d.get("source_file_id"),
            "drive_input_file_name": d.get("source_file_name"),
        },
        "available_chains": available_chains,
        "selected_chains": d["chain_ids"] or available_chains,
        "heterogen_policy": d["heterogen_policy"],
        "non_water_heterogens_detected": heterogens,
        "waters_detected": len(waters),
        "missing_residues_detected": missing_residues_detected,
        "missing_residues_built": d["add_missing_residues"],
        "nonstandard_replacements": nonstandard,
        "missing_atoms_added": missing_atoms_count,
        "missing_terminal_atoms_added": missing_terminals_count,
        "prepared_atom_count": atom_count,
        "prepared_residue_count": len(residue_keys),
        "ph": d["ph"],
        "receptor_pdbqt_method": "Meeko mk_prepare_receptor --read_pdb",
        "charge_model": "gasteiger",
        "meeko_receptor_options": ["--allow_bad_res"],
        "binding_site_reference": binding_site_reference,
        "binding_site_candidates": binding_site_candidates,
        "versions": {
            "meeko": str(getattr(meeko, "__version__", "unknown")),
            "rdkit": str(getattr(rdkit, "__version__", "unknown")),
            "openmm": str(getattr(openmm, "__version__", "unknown")),
        },
    }
    with open(os.path.join(out_dir, "protein_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=2)

    qc = [flag for flag in d["qc_flags"] if flag != "GRID_CENTER_TO_BE_SELECTED_OR_INFERRED"]
    qc.append("MEEKO_ALLOW_BAD_RES_USED_REQUIRES_REVIEW")
    if d["cx"] is None:
        if binding_site_reference:
            qc.append("GRID_CENTER_FROM_COCRYSTALLIZED_LIGAND_REQUIRES_REVIEW")
        else:
            qc.append("AUTO_CENTER_NOT_SCIENTIFICALLY_VALIDATED")
    if heterogens:
        qc.append("NON_WATER_HETEROGENS_REMOVED_EXPLICITLY")
    if waters and d["heterogen_policy"] == "remove_all":
        qc.append("CRYSTAL_WATERS_REMOVED_EXPLICITLY")
    if missing_residues_detected and not d["add_missing_residues"]:
        qc.append("MISSING_RESIDUES_NOT_MODELLED")
    if d["add_missing_residues"] and missing_residues_detected:
        qc.append("MODELLED_MISSING_RESIDUES_REQUIRE_REVIEW")
    out = dict(d)
    out["qc_flags"] = sorted(set(qc))
    out["protein_stats"] = stats
    out["receptor_source_path"] = source_cif
    out["prepared_receptor_path"] = prepared_pdb
    out["auto_grid_center"] = binding_site_reference["center_A"] if d["cx"] is None and binding_site_reference else None
    out["grid_center_source"] = ("user_supplied" if d["cx"] is not None else ("co_crystallized_ligand:" + binding_site_reference["component_id"] if binding_site_reference else "receptor_centroid_fallback"))
    return [{"json": out}]
