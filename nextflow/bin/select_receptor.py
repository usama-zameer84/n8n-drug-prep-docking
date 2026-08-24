#!/usr/bin/env python3
"""Select a receptor and record its preparation preflight."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request


def get_json(url, payload=None, attempts=3):
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": "ligand-docking-workbench/1.0 receptor-selection",
    }
    body = json.dumps(payload).encode("utf-8") if payload is not None else None
    last = None
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, data=body, headers=headers)
            with urllib.request.urlopen(req, timeout=35) as r:
                raw = r.read()
            # 204 No Content / empty body = definitive "no results" for this query; do not retry.
            if not raw or not raw.strip():
                raise ValueError(f"empty response (HTTP {r.status})")
            return json.loads(raw)
        except ValueError:
            raise
        except Exception as exc:
            last = exc
            if attempt + 1 < attempts:
                time.sleep(0.6 * (attempt + 1))
    raise ValueError(f"External evidence query failed for {url.split('?')[0]}: {last}")


def receptor_preflight(pdb_id, chain_ids, ph, work):
    result = {"pdb_id": pdb_id, "chain_ids": chain_ids, "passed": False}
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        result["reason"] = f"pdbfixer/openmm not available: {exc}"
        return result
    source_cif = os.path.join(work, pdb_id + ".cif")
    with urllib.request.urlopen(
        f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=120
    ) as r:
        open(source_cif, "wb").write(r.read())
    fixer = PDBFixer(filename=source_cif)
    chains = list(fixer.topology.chains())
    available = [c.id or str(i) for i, c in enumerate(chains)]
    missing = [c for c in chain_ids if c not in available]
    result["available_chains"] = available
    if missing:
        result["reason"] = (
            "UniProt/SIFTS chain(s) not exposed by PDBFixer: " + ",".join(missing)
        )
        return result
    remove = [i for i, c in enumerate(chains) if (c.id or str(i)) not in chain_ids]
    if remove:
        fixer.removeChains(remove)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(float(ph))
    prepared = os.path.join(work, "prepared.pdb")
    with open(prepared, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    pdbqt = os.path.join(work, "receptor.pdbqt")
    rjson = os.path.join(work, "receptor.json")
    cmd = [
        sys.executable,
        "-m",
        "meeko.cli.mk_prepare_receptor",
        "--read_pdb",
        prepared,
        "--write_pdbqt",
        pdbqt,
        "--write_json",
        rjson,
        "--charge_model",
        "gasteiger",
        "--allow_bad_res",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    result["returncode"] = proc.returncode
    result["log_tail"] = log[-600:]
    result["passed"] = (
        proc.returncode == 0 and os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 0
    )
    if not result["passed"]:
        result["reason"] = "PDBFixer-to-Meeko preflight failed"
    return result


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    args = ap.parse_args()
    d = json.load(open(args.params))
    mode = d.get("receptor_selection_mode", "auto")

    if mode == "provided":
        d["receptor_selection"] = {
            "mode": "provided",
            "status": "USER_SPECIFIED_RECEPTOR_USED",
            "method": "No ligand-to-target inference; the input explicitly selected a receptor.",
            "selected": {"pdb_id": d["pdb_id"], "chain_ids": d["chain_ids"]},
            "rationale": [f"The input explicitly supplied pdb_id={d['pdb_id']}."],
            "limitations": [
                "A supplied PDB ID is not evidence that the ligand binds this receptor."
            ],
        }
        json.dump(d, open(args.params, "w"), indent=2)
        return

    from rdkit import Chem

    query_mol = Chem.MolFromSmiles(d["smiles"])
    if query_mol is None:
        raise SystemExit(
            "RDKit could not parse the ligand for automatic receptor selection"
        )
    canonical = Chem.MolToSmiles(query_mol, isomericSmiles=True)
    threshold = float(d.get("target_similarity_threshold", 70))
    limit = int(d.get("target_candidate_limit", 5))
    organism = d.get("target_organism", "Homo sapiens")

    sim_url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{urllib.parse.quote(canonical, safe='')}/{int(threshold)}.json?limit=12"
    sim_data = get_json(sim_url)
    molecule_hits, seen = [], set()
    for item in sim_data.get("molecules", []):
        mid = item.get("molecule_chembl_id")
        if not mid or mid in seen:
            continue
        try:
            sim = float(item.get("similarity") or 0)
        except (TypeError, ValueError):
            sim = 0
        if sim < threshold:
            continue
        st = item.get("molecule_structures") or {}
        molecule_hits.append(
            {
                "molecule_chembl_id": mid,
                "similarity_percent": round(sim, 3),
                "canonical_smiles": st.get("canonical_smiles"),
            }
        )
        seen.add(mid)
        if len(molecule_hits) >= 8:
            break
    if not molecule_hits:
        raise SystemExit("No ChEMBL molecule above the configured similarity threshold")

    aggregates = {}
    for hit in molecule_hits:
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/activity.json?"
            + urllib.parse.urlencode(
                {"molecule_chembl_id": hit["molecule_chembl_id"], "limit": 300}
            )
        )
        for activity in get_json(url).get("activities", []):
            tid = activity.get("target_chembl_id")
            tname = activity.get("target_pref_name")
            if (
                not tid
                or tname in (None, "No relevant target")
                or activity.get("target_organism") != organism
            ):
                continue
            if activity.get("assay_type") not in ("B", "F") or activity.get(
                "data_validity_comment"
            ):
                continue
            rec = aggregates.setdefault(
                tid,
                {
                    "target_chembl_id": tid,
                    "target_name": tname,
                    "organism": organism,
                    "evidence_records": 0,
                    "similar_molecules": set(),
                    "max_similarity_percent": 0.0,
                    "pchembl_values": [],
                    "exact_query_evidence_records": 0,
                },
            )
            rec["evidence_records"] += 1
            rec["similar_molecules"].add(hit["molecule_chembl_id"])
            rec["max_similarity_percent"] = max(
                rec["max_similarity_percent"], hit["similarity_percent"]
            )
            if activity.get("pchembl_value") is not None:
                try:
                    rec["pchembl_values"].append(float(activity["pchembl_value"]))
                except (TypeError, ValueError):
                    pass
            if hit["similarity_percent"] >= 99.99:
                rec["exact_query_evidence_records"] += 1

    def score(r):
        return (
            0.42 * (r["max_similarity_percent"] / 100.0)
            + 0.20 * (min(r["evidence_records"], 20) / 20.0)
            + 0.13 * (min(len(r["similar_molecules"]), 5) / 5.0)
            + 0.10
            * (
                min(
                    max(
                        (max(r["pchembl_values"]) if r["pchembl_values"] else 0.0)
                        - 3.0,
                        0.0,
                    ),
                    6.0,
                )
                / 6.0
            )
            + 0.15 * (1.0 if r["exact_query_evidence_records"] else 0.0)
        )

    ranked = sorted(
        aggregates.values(), key=lambda r: (-score(r), r["target_chembl_id"])
    )[:12]
    if not ranked:
        raise SystemExit(
            "No organism-matched binding/functional target evidence in ChEMBL"
        )

    target_candidates = []
    for r in ranked:
        detail = get_json(
            f"https://www.ebi.ac.uk/chembl/api/data/target/{urllib.parse.quote(r['target_chembl_id'])}.json"
        )
        if detail.get("target_type") != "SINGLE PROTEIN":
            continue
        accessions = [
            c.get("accession")
            for c in (detail.get("target_components") or [])
            if c.get("accession")
        ]
        if not accessions:
            continue
        target_candidates.append(
            {
                "target_chembl_id": r["target_chembl_id"],
                "target_name": detail.get("pref_name") or r["target_name"],
                "organism": detail.get("organism") or r["organism"],
                "uniprot_accession": accessions[0],
                "evidence_records": r["evidence_records"],
                "similar_molecule_count": len(r["similar_molecules"]),
                "max_similarity_percent": round(r["max_similarity_percent"], 3),
                "evidence_score": round(score(r), 6),
            }
        )
    target_candidates.sort(key=lambda t: (-t["evidence_score"], t["target_chembl_id"]))
    target_candidates = target_candidates[:limit]
    if not target_candidates:
        raise SystemExit(
            "ChEMBL evidence did not resolve to an organism-matched SINGLE PROTEIN target"
        )

    structures = []
    for target in target_candidates:
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                            "operator": "exact_match",
                            "value": target["uniprot_accession"],
                        },
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                            "operator": "greater",
                            "value": 0,
                        },
                    },
                ],
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 8},
                "sort": [
                    {
                        "sort_by": "rcsb_entry_info.resolution_combined",
                        "direction": "asc",
                    }
                ],
            },
        }
        # RCSB returns HTTP 204 (No Content) when a target has no ligand-containing structures.
        # Skip that target and try the next ranked candidate instead of failing the whole run.
        try:
            response = get_json("https://search.rcsb.org/rcsbsearch/v2/query", query)
        except ValueError:
            continue
        for pdb_id in [
            i.get("identifier")
            for i in response.get("result_set", [])
            if i.get("identifier")
        ][:6]:
            try:
                detail = get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
                info = detail.get("rcsb_entry_info") or {}
                resolutions = info.get("resolution_combined") or []
                resolution = min(resolutions) if resolutions else None
                methods = [
                    m.get("method") for m in detail.get("exptl", []) if m.get("method")
                ]
                nonpoly = int(info.get("nonpolymer_entity_count") or 0)
                rq = (
                    max(0.0, min(1.0, (4.5 - float(resolution)) / 3.0))
                    if resolution is not None
                    else 0.0
                )
                combined = (
                    0.72 * target["evidence_score"]
                    + 0.20 * rq
                    + (0.08 if nonpoly > 0 else 0.0)
                )
                structures.append(
                    {
                        "target_chembl_id": target["target_chembl_id"],
                        "target_name": target["target_name"],
                        "organism": target["organism"],
                        "uniprot_accession": target["uniprot_accession"],
                        "pdb_id": pdb_id,
                        "resolution_A": round(float(resolution), 3)
                        if resolution is not None
                        else None,
                        "experimental_method": "; ".join(methods)
                        if methods
                        else "unknown",
                        "selection_score": round(combined, 6),
                        "polymer_entity_ids": (
                            detail.get("rcsb_entry_container_identifiers") or {}
                        ).get("polymer_entity_ids")
                        or [],
                    }
                )
            except Exception:
                pass
    if not structures:
        raise SystemExit("No ligand-containing experimental RCSB structure found")
    structures.sort(
        key=lambda s: (
            -s["selection_score"],
            s["resolution_A"] if s["resolution_A"] is not None else 99.0,
            s["pdb_id"],
        )
    )

    preflight = []
    selected = None
    selected_chains = []
    for structure in structures[:12]:
        chains = []
        for entity_id in structure.get("polymer_entity_ids", []):
            try:
                entity = get_json(
                    f"https://data.rcsb.org/rest/v1/core/polymer_entity/{structure['pdb_id']}/{entity_id}"
                )
                ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                refs = ids.get("reference_sequence_identifiers") or []
                if any(
                    ref.get("database_accession") == structure["uniprot_accession"]
                    for ref in refs
                ):
                    chains.extend(ids.get("auth_asym_ids") or [])
            except Exception:
                continue
        chains = sorted(set(c for c in chains if c))
        if not chains:
            preflight.append(
                {
                    "pdb_id": structure["pdb_id"],
                    "chain_ids": [],
                    "passed": False,
                    "reason": "No SIFTS mapping",
                }
            )
            continue
        work = tempfile.mkdtemp(prefix="receptor_candidate_")
        check = receptor_preflight(structure["pdb_id"], chains, d["ph"], work)
        shutil.rmtree(work, ignore_errors=True)
        preflight.append(check)
        if check["passed"]:
            selected = dict(structure)
            selected_chains = chains
            selected["chain_ids"] = chains
            selected["preparation_preflight"] = "PASSED_PDBFIXER_MEEKO"
            break
    if selected is None:
        raise SystemExit(
            "No ranked receptor structure passed the bounded PDBFixer-to-Meeko preflight"
        )

    d["pdb_id"] = selected["pdb_id"]
    d["chain_ids"] = selected_chains
    d["run_name"] = d["pdb_id"] + "_" + d["ligand_name"]
    d["receptor_selection"] = {
        "mode": "auto",
        "status": "AUTOMATIC_HYPOTHESIS_SELECTED_FOR_DOCKING",
        "query": {
            "input_smiles": d["smiles"],
            "canonical_smiles": canonical,
            "target_organism": organism,
            "similarity_threshold_percent": threshold,
            "target_candidate_limit": limit,
        },
        "method": "ChEMBL ligand-similarity target fishing + RCSB ranking + bounded PDBFixer-to-Meeko preflight.",
        "preparation_preflight": preflight,
        "candidate_molecules": molecule_hits,
        "candidate_targets": target_candidates,
        "selected": {k: v for k, v in selected.items() if k != "polymer_entity_ids"},
        "data_sources": {
            "chembl_similarity": sim_url,
            "rcsb_search_api": "https://search.rcsb.org/rcsbsearch/v2/query",
        },
    }
    d["qc_flags"] = sorted(
        set(d.get("qc_flags", []) + ["AUTOMATED_RECEPTOR_SELECTION_REQUIRES_REVIEW"])
    )
    json.dump(d, open(args.params, "w"), indent=2)


if __name__ == "__main__":
    main()
