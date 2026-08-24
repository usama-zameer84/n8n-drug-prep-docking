#!/usr/bin/env python3
"""Build the report data and write the per-run report set."""

import argparse
import json
import shutil
from pathlib import Path

from report_renderer import ReportRenderer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--templates", required=True, type=Path)
    parser.add_argument("--structure-dir", type=Path)
    parser.add_argument("--handoff-dir", type=Path)
    args = parser.parse_args()

    data = json.loads(args.params.read_text())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.handoff_dir:
        structure_dir = args.out_dir / "MD_Handoff"
        shutil.copytree(args.handoff_dir, structure_dir, dirs_exist_ok=True)
    elif args.structure_dir:
        structure_dir = args.structure_dir
    else:
        raise SystemExit("Use --handoff-dir or --structure-dir")
    pages = ReportRenderer(args.templates).render_all(data, structure_dir)
    for name, content in pages.items():
        (args.out_dir / name).write_text(content, encoding="utf-8")

    docking = data.get("docking_results", {}) or {}
    summary = {
        "run_id": data.get("run_id"),
        "ligand_id": data.get("ligand_id"),
        "smiles": data.get("smiles"),
        "pdb_id": data.get("pdb_id"),
        "best_affinity_kcal_mol": docking.get(
            "top_affinity_kcal_mol", docking.get("best_affinity_kcal_mol")
        ),
        "grid_center_source": data.get("grid_center_source"),
        "qc_flags": data.get("qc_flags", []),
        "md_status": "STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY",
    }
    summary_name = "run_summary.json"
    (args.out_dir / summary_name).write_text(json.dumps(summary, indent=2))
    data["report_data"] = {"files": [*pages, summary_name]}
    args.params.write_text(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
