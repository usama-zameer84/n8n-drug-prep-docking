"""n8n Code node: Auto Select Receptor.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import json
    import os
    import shutil
    import subprocess
    import sys
    import tempfile
    import time
    import urllib.parse
    import urllib.request
    from rdkit import Chem
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    d = dict(_items[0]["json"])
    mode = d.get("receptor_selection_mode", "auto")

    def get_json(url, payload=None, attempts=3):
        headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "n8n-drug-prep-docking/3.0 receptor-selection"}
        body = json.dumps(payload).encode("utf-8") if payload is not None else None
        last_error = None
        for attempt in range(attempts):
            try:
                req = urllib.request.Request(url, data=body, headers=headers)
                with urllib.request.urlopen(req, timeout=35) as response:
                    return json.load(response)
            except Exception as exc:
                last_error = exc
                if attempt + 1 < attempts:
                    time.sleep(0.6 * (attempt + 1))
        raise ValueError("External evidence query failed for " + url.split("?")[0] + ": " + str(last_error))

    if mode == "provided":
        d["receptor_selection"] = {
            "mode": "provided",
            "status": "USER_SPECIFIED_RECEPTOR_USED",
            "query": {"input_smiles": d["smiles"], "target_organism": d["target_organism"], "provided_pdb_id": d["pdb_id"], "provided_chain_ids": d["chain_ids"]},
            "method": "No ligand-to-target inference was performed because the input explicitly selected a receptor.",
            "candidate_molecules": [],
            "candidate_targets": [],
            "selected": {"pdb_id": d["pdb_id"], "chain_ids": d["chain_ids"], "target_name": None, "target_chembl_id": None, "uniprot_accession": None, "organism": None, "resolution_A": None, "experimental_method": None, "selection_score": None},
            "rationale": ["The input explicitly supplied pdb_id=" + d["pdb_id"] + ".", "The supplied chain and grid settings were passed to receptor preparation and docking."],
            "limitations": ["A supplied PDB ID is not evidence that the ligand binds this receptor.", "The receptor identity, biological assembly, chain, construct and binding site remain the user's scientific responsibility."],
            "data_sources": {"rcsb": "https://www.rcsb.org/", "chembl": "not queried in provided mode"},
        }
        return [{"json": d}]

    query_mol = Chem.MolFromSmiles(d["smiles"])
    if query_mol is None:
        raise ValueError("RDKit could not parse the ligand for automatic receptor selection")
    canonical = Chem.MolToSmiles(query_mol, isomericSmiles=True)
    threshold = float(d.get("target_similarity_threshold", 70))
    limit = int(d.get("target_candidate_limit", 5))
    organism = d.get("target_organism", "Homo sapiens")
    similarity_url = "https://www.ebi.ac.uk/chembl/api/data/similarity/" + urllib.parse.quote(canonical, safe="") + "/" + str(int(threshold)) + ".json?limit=12"
    similarity_data = get_json(similarity_url)
    molecule_hits = []
    seen_molecules = set()
    for item in similarity_data.get("molecules", []):
        molecule_id = item.get("molecule_chembl_id")
        if not molecule_id or molecule_id in seen_molecules:
            continue
        try:
            similarity = float(item.get("similarity") or 0)
        except (TypeError, ValueError):
            similarity = 0
        if similarity < threshold:
            continue
        structures = item.get("molecule_structures") or {}
        molecule_hits.append({"molecule_chembl_id": molecule_id, "similarity_percent": round(similarity, 3), "canonical_smiles": structures.get("canonical_smiles")})
        seen_molecules.add(molecule_id)
        if len(molecule_hits) >= 8:
            break
    if not molecule_hits:
        raise ValueError("Automatic receptor selection found no ChEMBL molecule above the configured similarity threshold")

    def load_activities(hit):
        url = "https://www.ebi.ac.uk/chembl/api/data/activity.json?" + urllib.parse.urlencode({"molecule_chembl_id": hit["molecule_chembl_id"], "limit": 300})
        return {"hit": hit, "activities": get_json(url).get("activities", [])}
    activity_batches = [load_activities(hit) for hit in molecule_hits]

    aggregates = {}
    for activity_batch in activity_batches:
        hit = activity_batch["hit"]
        activities = activity_batch["activities"]
        for activity in activities:
            target_id = activity.get("target_chembl_id")
            target_name = activity.get("target_pref_name")
            if not target_id or target_name in (None, "No relevant target") or activity.get("target_organism") != organism:
                continue
            if activity.get("assay_type") not in ("B", "F") or activity.get("data_validity_comment"):
                continue
            record = aggregates.setdefault(target_id, {"target_chembl_id": target_id, "target_name": target_name, "organism": organism, "evidence_records": 0, "similar_molecules": set(), "max_similarity_percent": 0.0, "pchembl_values": [], "assay_types": set(), "standard_types": set(), "exact_query_evidence_records": 0})
            record["evidence_records"] += 1
            record["similar_molecules"].add(hit["molecule_chembl_id"])
            record["max_similarity_percent"] = max(record["max_similarity_percent"], hit["similarity_percent"])
            record["assay_types"].add(str(activity.get("assay_type")))
            if activity.get("standard_type"):
                record["standard_types"].add(str(activity.get("standard_type")))
            try:
                if activity.get("pchembl_value") is not None:
                    record["pchembl_values"].append(float(activity["pchembl_value"]))
            except (TypeError, ValueError):
                pass
            if hit["similarity_percent"] >= 99.99:
                record["exact_query_evidence_records"] += 1

    preliminary_score = lambda record: (
        0.42 * (record["max_similarity_percent"] / 100.0)
        + 0.20 * (min(record["evidence_records"], 20) / 20.0)
        + 0.13 * (min(len(record["similar_molecules"]), 5) / 5.0)
        + 0.10 * (min(max((max(record["pchembl_values"]) if record["pchembl_values"] else 0.0) - 3.0, 0.0), 6.0) / 6.0)
        + 0.15 * (1.0 if record["exact_query_evidence_records"] else 0.0)
    )

    ranked_raw = sorted(aggregates.values(), key=lambda item: (-preliminary_score(item), item["target_chembl_id"]))[:12]
    if not ranked_raw:
        raise ValueError("Automatic receptor selection found no organism-matched binding or functional target evidence in ChEMBL")

    def load_target(record):
        return {"record": record, "detail": get_json("https://www.ebi.ac.uk/chembl/api/data/target/" + urllib.parse.quote(record["target_chembl_id"]) + ".json")}
    target_candidates = []
    for target_result in [load_target(record) for record in ranked_raw]:
        record = target_result["record"]
        detail = target_result["detail"]
        if detail.get("target_type") != "SINGLE PROTEIN":
            continue
        accessions = [component.get("accession") for component in (detail.get("target_components") or []) if component.get("accession")]
        if not accessions:
            continue
        values = record["pchembl_values"]
        target_candidates.append({"target_chembl_id": record["target_chembl_id"], "target_name": detail.get("pref_name") or record["target_name"], "organism": detail.get("organism") or record["organism"], "target_type": detail.get("target_type"), "uniprot_accession": accessions[0], "evidence_records": record["evidence_records"], "similar_molecule_count": len(record["similar_molecules"]), "max_similarity_percent": round(record["max_similarity_percent"], 3), "exact_query_evidence_records": record["exact_query_evidence_records"], "max_pchembl": round(max(values), 3) if values else None, "assay_types": sorted(record["assay_types"]), "standard_types": sorted(record["standard_types"]), "evidence_score": round(preliminary_score(record), 6)})
    target_candidates.sort(key=lambda item: (-item["evidence_score"], item["target_chembl_id"]))
    target_candidates = target_candidates[:limit]
    if not target_candidates:
        raise ValueError("ChEMBL evidence did not resolve to an organism-matched SINGLE PROTEIN target with a UniProt accession")

    def search_structures(target):
        accession = target["uniprot_accession"]
        query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "exact_match", "value": accession}},
            {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.nonpolymer_entity_count", "operator": "greater", "value": 0}}
        ]}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 8}, "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]}}
        response = get_json("https://search.rcsb.org/rcsbsearch/v2/query", query)
        return {"target": target, "pdb_ids": [item.get("identifier") for item in response.get("result_set", []) if item.get("identifier")]}
    structure_requests = [search_structures(target) for target in target_candidates]

    def load_entry(target, pdb_id):
        detail = get_json("https://data.rcsb.org/rest/v1/core/entry/" + pdb_id)
        info = detail.get("rcsb_entry_info") or {}
        resolutions = info.get("resolution_combined") or []
        resolution = min(resolutions) if resolutions else None
        methods = [item.get("method") for item in detail.get("exptl", []) if item.get("method")]
        nonpolymer_count = int(info.get("nonpolymer_entity_count") or 0)
        resolution_quality = max(0.0, min(1.0, (4.5 - float(resolution)) / 3.0)) if resolution is not None else 0.0
        combined = 0.72 * target["evidence_score"] + 0.20 * resolution_quality + (0.08 if nonpolymer_count > 0 else 0.0)
        return {"target_chembl_id": target["target_chembl_id"], "target_name": target["target_name"], "organism": target["organism"], "uniprot_accession": target["uniprot_accession"], "pdb_id": pdb_id, "resolution_A": round(float(resolution), 3) if resolution is not None else None, "experimental_method": "; ".join(methods) if methods else "unknown", "nonpolymer_entity_count": nonpolymer_count, "selection_score": round(combined, 6), "polymer_entity_ids": (detail.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids") or []}
    structures = []
    for request in structure_requests:
        for pdb_id in request["pdb_ids"][:6]:
            try:
                structures.append(load_entry(request["target"], pdb_id))
            except Exception:
                pass
    if not structures:
        raise ValueError("No ligand-containing experimental RCSB structure was found for the ranked target candidates")
    structures.sort(key=lambda item: (-item["selection_score"], item["resolution_A"] if item["resolution_A"] is not None else 99.0, item["pdb_id"]))

    def mapped_target_chains(structure):
        chains = []
        for entity_id in structure.get("polymer_entity_ids", []):
            try:
                entity = get_json("https://data.rcsb.org/rest/v1/core/polymer_entity/" + structure["pdb_id"] + "/" + str(entity_id))
            except Exception:
                continue
            identifiers = entity.get("rcsb_polymer_entity_container_identifiers") or {}
            refs = identifiers.get("reference_sequence_identifiers") or []
            if any(ref.get("database_accession") == structure["uniprot_accession"] for ref in refs):
                chains.extend(identifiers.get("auth_asym_ids") or [])
        return {"chain_ids": sorted(set(chain for chain in chains if chain))}

    def receptor_preflight(structure, chain_ids):
        work = tempfile.mkdtemp(prefix="receptor_candidate_", dir="/md_project/data")
        result = {"pdb_id": structure["pdb_id"], "chain_ids": chain_ids, "passed": False}
        try:
            source_cif = os.path.join(work, structure["pdb_id"] + ".cif")
            with urllib.request.urlopen("https://files.rcsb.org/download/" + structure["pdb_id"] + ".cif", timeout=120) as response:
                with open(source_cif, "wb") as fh:
                    fh.write(response.read())
            fixer = PDBFixer(filename=source_cif)
            chains = list(fixer.topology.chains())
            available = [chain.id or str(index) for index, chain in enumerate(chains)]
            missing = [chain for chain in chain_ids if chain not in available]
            result["available_chains"] = available
            if missing:
                result["reason"] = "UniProt/SIFTS author chain(s) not exposed by PDBFixer: " + ",".join(missing)
                return result
            remove_indices = [index for index, chain in enumerate(chains) if (chain.id or str(index)) not in chain_ids]
            if remove_indices:
                fixer.removeChains(remove_indices)
            fixer.findMissingResidues()
            fixer.missingResidues = {}
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.removeHeterogens(False)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(float(d["ph"]))
            prepared_pdb = os.path.join(work, "prepared.pdb")
            with open(prepared_pdb, "w") as fh:
                PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
            receptor_pdbqt = os.path.join(work, "receptor.pdbqt")
            receptor_json = os.path.join(work, "receptor.json")
            cmd = [sys.executable, "-m", "meeko.cli.mk_prepare_receptor", "--read_pdb", prepared_pdb, "--write_pdbqt", receptor_pdbqt, "--write_json", receptor_json, "--charge_model", "gasteiger", "--allow_bad_res"]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
            result["returncode"] = proc.returncode
            result["log_tail"] = log_text[-600:]
            result["passed"] = proc.returncode == 0 and os.path.exists(receptor_pdbqt) and os.path.getsize(receptor_pdbqt) > 0
            if not result["passed"]:
                result["reason"] = "PDBFixer-to-Meeko preparation preflight failed"
            return result
        except Exception as exc:
            result["reason"] = str(exc)
            return result
        finally:
            shutil.rmtree(work, ignore_errors=True)

    preparation_preflight = []
    selected_structure = None
    selected_chains = []
    for structure in structures[:12]:
        chain_result = mapped_target_chains(structure)
        candidate_chains = chain_result["chain_ids"]
        if not candidate_chains:
            preparation_preflight.append({"pdb_id": structure["pdb_id"], "chain_ids": [], "passed": False, "reason": "No UniProt/SIFTS author-chain mapping"})
            continue
        check = receptor_preflight(structure, candidate_chains)
        preparation_preflight.append(check)
        if check["passed"]:
            selected_structure = dict(structure)
            selected_chains = candidate_chains
            selected_structure["chain_ids"] = candidate_chains
            selected_structure["preparation_preflight"] = "PASSED_PDBFIXER_MEEKO"
            break
    if selected_structure is None:
        raise ValueError("No ranked receptor structure passed the bounded PDBFixer-to-Meeko preparation preflight")
    selected_target = next(item for item in target_candidates if item["target_chembl_id"] == selected_structure["target_chembl_id"])
    for target in target_candidates:
        matches = [item for item in structures if item["target_chembl_id"] == target["target_chembl_id"]]
        target["best_structure"] = ({key: value for key, value in matches[0].items() if key != "polymer_entity_ids"} if matches else None)
    selected_public = {key: value for key, value in selected_structure.items() if key != "polymer_entity_ids"}
    selected_public["chain_ids"] = selected_chains
    d["pdb_id"] = selected_structure["pdb_id"]
    d["chain_ids"] = selected_chains
    d["run_name"] = d["pdb_id"] + "_" + d["ligand_name"] + "_" + d["run_id"]
    d["receptor_selection"] = {
        "mode": "auto",
        "status": "AUTOMATIC_HYPOTHESIS_SELECTED_FOR_DOCKING",
        "query": {"input_smiles": d["smiles"], "canonical_smiles": canonical, "target_organism": organism, "similarity_threshold_percent": threshold, "target_candidate_limit": limit},
        "method": "ChEMBL ligand-similarity target fishing using organism-matched binding/functional activity evidence, followed by RCSB experimental structure ranking and bounded PDBFixer-to-Meeko preparation preflight.",
        "preparation_preflight": preparation_preflight,
        "candidate_molecules": molecule_hits,
        "candidate_targets": target_candidates,
        "selected": selected_public,
        "rationale": ["The selected target had the highest combined ligand-similarity, curated activity-evidence and experimental-structure quality score.", "The selected PDB entry contains non-polymer component(s), has recorded experimental resolution, and maps to the target UniProt accession.", "Target chain selection used the RCSB UniProt/SIFTS polymer-entity mapping.", "The highest-ranked structure that passed an actual PDBFixer-to-Meeko preparation preflight was selected; rejected structures and reasons are retained.", "The chosen receptor and chain were passed automatically to PDBFixer, official Meeko receptor preparation and AutoDock Vina."],
        "limitations": ["Ligand-based target fishing generates hypotheses and cannot establish the true biological receptor.", "ChEMBL evidence may be sparse, biased toward well-studied chemotypes, or measured in heterogeneous assays.", "A high-quality PDB structure may contain engineered mutations, fusion partners, stabilizing ligands or a non-native conformational state.", "The selected receptor, pocket and pose require expert review and experimental validation."],
        "data_sources": {"chembl_similarity": similarity_url, "chembl_api": "https://www.ebi.ac.uk/chembl/api/data/", "rcsb_search_api": "https://search.rcsb.org/rcsbsearch/v2/query", "rcsb_data_api": "https://data.rcsb.org/rest/v1/core/"}
    }
    d["qc_flags"] = sorted(set(list(d.get("qc_flags", [])) + ["AUTOMATED_RECEPTOR_SELECTION_REQUIRES_REVIEW"]))
    return [{"json": d}]
