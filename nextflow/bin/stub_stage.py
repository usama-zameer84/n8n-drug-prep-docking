#!/usr/bin/env python3
"""Create deterministic placeholder outputs for Nextflow stub runs."""

import argparse
import base64
import json
import shutil
from pathlib import Path


REPORT_FILES = (
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
)


def load_params(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_params(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def touch_files(*names: str) -> None:
    for name in names:
        Path(name).touch()


def normalize(args: argparse.Namespace) -> None:
    try:
        raw = base64.b64decode(args.input_base64, validate=True)
        data = json.loads(raw)
    except (ValueError, json.JSONDecodeError) as error:
        raise SystemExit(f"Invalid stub input: {error}") from error
    data["qc_flags"] = ["STUB"]
    save_params(args.params, data)


def select_receptor(args: argparse.Namespace) -> None:
    data = load_params(args.params)
    data.setdefault(
        "receptor_selection",
        {"mode": data.get("receptor_selection_mode", "auto"), "status": "STUB"},
    )
    save_params(args.params, data)


def generate_3d(args: argparse.Namespace) -> None:
    touch_files("ligand_input.sdf", "ligand_2d.svg")
    data = load_params(args.params)
    data.update(
        {
            "canonical_smiles": data["smiles"],
            "inchi_key": "STUB",
            "num_heavy_atoms": 0,
            "num_atoms": 0,
            "conformer_generation": {"method": "stub"},
        }
    )
    save_params(args.params, data)


def analyze_ro5(args: argparse.Namespace) -> None:
    data = load_params(args.params)
    data["ro5"] = {"druglike": True, "violations": 0, "STUB": True}
    save_params(args.params, data)


def convert_pdbqt(args: argparse.Namespace) -> None:
    touch_files("ligand.pdbqt")
    data = load_params(args.params)
    data["ligand_prep"] = {"method": "stub"}
    save_params(args.params, data)


def prepare_protein(args: argparse.Namespace) -> None:
    touch_files(
        "prepared_receptor.pdb",
        "receptor.pdbqt",
        "receptor.json",
        "protein_stats.json",
        "receptor_source.cif",
        "receptor_prep.log",
    )
    data = load_params(args.params)
    data["protein_stats"] = {"STUB": True}
    data["auto_grid_center"] = None
    data["grid_center_source"] = "receptor_centroid_fallback"
    save_params(args.params, data)


def run_docking(args: argparse.Namespace) -> None:
    touch_files(
        "docked.pdbqt",
        "docking.log",
        "docked_replica_1.pdbqt",
        "vina_replica_1.log",
    )
    data = load_params(args.params)
    data["docking_results"] = {"STUB": True}
    data["grid_center"] = [0, 0, 0]
    save_params(args.params, data)


def parse_results(args: argparse.Namespace) -> None:
    Path("docking_results.json").write_text(
        json.dumps({"STUB": True}), encoding="utf-8"
    )
    data = load_params(args.params)
    data["docking_results"] = {"STUB": True}
    save_params(args.params, data)


def analyze_interactions(args: argparse.Namespace) -> None:
    Path("interactions.json").write_text(
        json.dumps({"STUB": True}), encoding="utf-8"
    )
    data = load_params(args.params)
    data["interactions_summary"] = {"STUB": True}
    save_params(args.params, data)


def build_md_handoff(args: argparse.Namespace) -> None:
    output = Path("MD_Handoff")
    output.mkdir(parents=True, exist_ok=True)
    (output / "README_MD_HANDOFF.md").touch()
    status = "STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY"
    (output / "manifest.json").write_text(
        json.dumps({"status": status, "STUB": True}), encoding="utf-8"
    )
    data = load_params(args.params)
    data["md_handoff"] = {"status": status, "STUB": True}
    save_params(args.params, data)


def generate_reports(args: argparse.Namespace) -> None:
    output = Path("report_files")
    output.mkdir(parents=True, exist_ok=True)
    if args.input_dir:
        shutil.copytree(args.input_dir, output / "MD_Handoff", dirs_exist_ok=True)
    for name in REPORT_FILES:
        (output / name).write_text(f"stub {name}\n", encoding="utf-8")
    data = load_params(args.params)
    data["report_data"] = {"STUB": True}
    save_params(args.params, data)


def build_package(args: argparse.Namespace) -> None:
    output = Path("package")
    if output.exists():
        shutil.rmtree(output)
    (output / "MD_Handoff").mkdir(parents=True)
    (output / "index.html").write_text("stub\n", encoding="utf-8")
    (output / "run_summary.json").write_text("stub\n", encoding="utf-8")
    pose = f"stub_{args.ligand_id}\n$$$$\n"
    (output / "all_best_poses.sdf").write_text(pose, encoding="utf-8")
    (output / "MD_Handoff" / "best_pose_ligand.sdf").write_text(
        pose, encoding="utf-8"
    )


HANDLERS = {
    "normalize": normalize,
    "select-receptor": select_receptor,
    "generate-3d": generate_3d,
    "analyze-ro5": analyze_ro5,
    "convert-pdbqt": convert_pdbqt,
    "prepare-protein": prepare_protein,
    "run-docking": run_docking,
    "parse-results": parse_results,
    "analyze-interactions": analyze_interactions,
    "build-md-handoff": build_md_handoff,
    "generate-reports": generate_reports,
    "build-package": build_package,
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=HANDLERS)
    parser.add_argument("--params", type=Path, default=Path("params.json"))
    parser.add_argument("--input-base64")
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--ligand-id")
    args = parser.parse_args()

    if args.stage == "normalize" and not args.input_base64:
        parser.error("normalize requires --input-base64")
    if args.stage == "build-package" and not args.ligand_id:
        parser.error("build-package requires --ligand-id")
    HANDLERS[args.stage](args)


if __name__ == "__main__":
    main()
