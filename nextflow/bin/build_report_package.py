#!/usr/bin/env python3
"""Validate and assemble the published report package."""

import argparse
import hashlib
import json
import os
import shutil
from html.parser import HTMLParser


class ReportHTMLParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tags = set()
        self.local_links = []
        self.script_sources = []

    def handle_starttag(self, tag, attrs):
        self.tags.add(tag)
        attributes = dict(attrs)
        if tag == "a" and attributes.get("href"):
            self.local_links.append(attributes["href"])
        if tag == "script" and attributes.get("src"):
            self.script_sources.append(attributes["src"])


def validate_reports(out_dir, html_names):
    required_tags = {"html", "head", "title", "body"}
    template_values = (
        "$title",
        "$payload",
        "$ligand_name",
        "$run_id",
        "$table_rows",
        "$qc_flags",
        "$pdb_id",
        "$chains",
        "$selection_score",
        "$candidate_rows",
        "$preflight_rows",
        "$method",
        "$pose_rows",
        "$receptor_atoms",
        "$ligand_atoms",
    )
    for name in html_names:
        source = open(os.path.join(out_dir, name), encoding="utf-8").read()
        if not source.lower().lstrip().startswith("<!doctype html>"):
            raise SystemExit(f"Invalid HTML report {name}: missing HTML5 doctype")
        if any(value in source for value in template_values):
            raise SystemExit(f"Invalid HTML report {name}: unresolved template value")
        parser = ReportHTMLParser()
        parser.feed(source)
        parser.close()
        missing_tags = required_tags - parser.tags
        if missing_tags:
            raise SystemExit(
                f"Invalid HTML report {name}: missing {', '.join(sorted(missing_tags))}"
            )
        for target in parser.local_links:
            if "://" in target or target.startswith(("#", "mailto:")):
                continue
            target_path = target.split("#", 1)[0].split("?", 1)[0]
            if target_path and not os.path.exists(os.path.join(out_dir, target_path)):
                raise SystemExit(f"Broken report link in {name}: {target}")

    visualization = open(
        os.path.join(out_dir, "09_visualization.html"), encoding="utf-8"
    ).read()
    if "fetch(" in visualization or "<script src=" in visualization:
        raise SystemExit("Visualization must not depend on network or file fetches")
    if "receptorAtoms = [" not in visualization or "ligandAtoms = [" not in visualization:
        raise SystemExit("Visualization does not contain embedded structure data")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--report-dir")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    d = json.load(open(args.params))
    if args.report_dir:
        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        shutil.copytree(args.report_dir, args.out_dir)
    else:
        os.makedirs(args.out_dir, exist_ok=True)

    md_src = os.path.join(args.out_dir, "MD_Handoff")
    required = [
        "index.html",
        "01_input_provenance.html",
        "02_ligand_preparation.html",
        "03_receptor_preparation.html",
        "03a_receptor_selection.html",
        "04_docking_results.html",
        "05_distance_contacts.html",
        "06_methods_qc.html",
        "07_md_handoff_readiness.html",
        "08_raw_data.html",
        "09_visualization.html",
        "run_summary.json",
        "MD_Handoff/best_pose_ligand.sdf",
        "MD_Handoff/prepared_receptor.pdb",
        "MD_Handoff/manifest.json",
    ]
    missing = [
        name
        for name in required
        if not os.path.exists(os.path.join(args.out_dir, name))
    ]
    if missing:
        raise SystemExit("Report package inputs are missing: " + ", ".join(missing))

    html_names = [name for name in required if name.endswith(".html")]
    validate_reports(args.out_dir, html_names)
    summary = json.load(open(os.path.join(args.out_dir, "run_summary.json")))
    summary_fields = {
        "run_id",
        "ligand_id",
        "smiles",
        "pdb_id",
        "best_affinity_kcal_mol",
        "grid_center_source",
        "qc_flags",
        "md_status",
    }
    if summary_fields - summary.keys():
        raise SystemExit("run_summary.json is missing required fields")

    # all_best_poses.sdf: the chemistry-authoritative best pose
    best_sdf = os.path.join(md_src, "best_pose_ligand.sdf")
    shutil.copyfile(best_sdf, os.path.join(args.out_dir, "all_best_poses.sdf"))

    # receptor_overlay_reference.pdb: the prepared receptor for overlay comparisons
    receptor_path = os.path.join(md_src, "prepared_receptor.pdb")
    shutil.copyfile(
        receptor_path, os.path.join(args.out_dir, "receptor_overlay_reference.pdb")
    )

    # ligand_overlay.json: the overlay-group fingerprint (per docs/OUTPUTS.md)
    profile = {
        "pdb_id": d.get("pdb_id"),
        "chains": d.get("chain_ids"),
        "heterogen_policy": d.get("heterogen_policy"),
        "ph": d.get("ph"),
        "receptor_source_hash": (d.get("protein_stats") or {})
        .get("source", {})
        .get("sha256"),
    }
    profile_key = json.dumps(profile, sort_keys=True, separators=(",", ":")).encode()
    overlay = {
        "schema_version": "1.0",
        "overlay_group_id": hashlib.sha256(profile_key).hexdigest()[:16],
        "receptor_profile": profile,
        "entries": [
            {
                "ligand_id": d.get("ligand_id"),
                "ligand_name": d.get("ligand_name"),
                "run_id": d.get("run_id"),
                "score_kcal_mol": (d.get("docking_results") or {}).get(
                    "top_affinity_kcal_mol"
                ),
                "pose_file": "all_best_poses.sdf",
                "color": "#2f80ed",
            }
        ],
    }
    json.dump(
        overlay, open(os.path.join(args.out_dir, "ligand_overlay.json"), "w"), indent=2
    )


if __name__ == "__main__":
    main()
