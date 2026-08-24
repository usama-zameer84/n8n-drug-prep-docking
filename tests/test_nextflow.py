"""Behavioral and structural tests for the Nextflow workflow."""

import glob
import base64
import json
import py_compile
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NF = ROOT / "nextflow"
REPORT_PAGES = [
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
]


def run_helper(name, *args, cwd):
    return subprocess.run(
        [sys.executable, str(NF / "bin" / name), *map(str, args)],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def pdbqt_atom(serial, name, residue, chain, resid, x, y, z, atom_type):
    return (
        f"ATOM  {serial:5d} {name:<4} {residue:>3} {chain}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00    0.000 {atom_type:>2}\n"
    )


class NextflowTreeTests(unittest.TestCase):
    def test_process_files_present(self):
        names = sorted(path.name for path in (NF / "processes").glob("*.nf"))
        self.assertEqual(len(names), 12)
        prefixes = {name.split("_")[0] for name in names}
        self.assertEqual(
            prefixes,
            {
                "02",
                "03",
                "04",
                "05",
                "06",
                "07",
                "08",
                "09",
                "10",
                "11",
                "12",
                "13",
            },
        )

    def test_main_and_config_safety_contracts(self):
        main = (NF / "main.nf").read_text()
        config = (NF / "nextflow.config").read_text()
        normalize = (NF / "processes" / "02_normalize_input.nf").read_text()
        self.assertIn("Input must contain exactly one ligand record", main)
        self.assertIn("conda.enabled = true", config)
        self.assertIn("--input-base64", normalize)
        self.assertNotIn('--smiles "$smiles"', normalize)

    def test_processes_only_invoke_external_python_helpers(self):
        helper_mapping = {
            "02_normalize_input.nf": "normalize_input.py",
            "03_auto_select_receptor.nf": "select_receptor.py",
            "04_generate_3d_structure.nf": "generate_3d_structure.py",
            "05_analyze_ligand_ro5.nf": "analyze_ligand_ro5.py",
            "06_convert_to_pdbqt.nf": "convert_ligand_to_pdbqt.py",
            "07_prepare_protein.nf": "prepare_protein.py",
            "08_run_docking.nf": "run_docking.py",
            "09_parse_results.nf": "parse_results.py",
            "10_analyze_interactions.nf": "analyze_interactions.py",
            "11_build_md_handoff.nf": "build_md_handoff.py",
            "12_generate_reports.nf": "generate_reports.py",
            "13_build_report_package.nf": "build_report_package.py",
        }
        for name, helper in helper_mapping.items():
            source = (NF / "processes" / name).read_text()
            self.assertIn(f"bin/{helper}", source)
            self.assertNotIn("python -c", source)
            self.assertNotIn("<<", source)
            self.assertNotIn("import json", source)
            self.assertNotIn("json.load(", source)
            self.assertNotIn("json.dump(", source)

    def test_bin_helpers_compile(self):
        for helper in glob.glob(str(NF / "bin" / "*.py")):
            py_compile.compile(helper, doraise=True)

    def test_examples_present(self):
        self.assertTrue((NF / "examples" / "input.smi").is_file())
        self.assertTrue((NF / "examples" / "input.test.smi").is_file())

    def test_normalizer_enforces_bounds_and_json_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            bad = run_helper(
                "normalize_input.py",
                "--ligand-id",
                "x",
                "--smiles",
                "CC",
                "--ph",
                "3",
                "--out",
                tmp_path / "bad.json",
                cwd=tmp,
            )
            self.assertNotEqual(bad.returncode, 0)
            payload = {
                "ligand_id": "safe_name",
                "smiles": "CC",
                "target_organism": "$(touch sentinel) `touch other` Homo sapiens",
                "allow_multicomponent": False,
            }
            (tmp_path / "input.json").write_text(json.dumps(payload))
            good = run_helper(
                "normalize_input.py",
                "--input-json",
                tmp_path / "input.json",
                "--out",
                tmp_path / "params.json",
                cwd=tmp,
            )
            self.assertEqual(good.returncode, 0, good.stderr)
            result = json.loads((tmp_path / "params.json").read_text())
            self.assertEqual(result["target_organism"], payload["target_organism"])
            self.assertFalse((tmp_path / "sentinel").exists())
            self.assertFalse((tmp_path / "other").exists())

            encoded = base64.b64encode(json.dumps(payload).encode()).decode()
            base64_result = run_helper(
                "normalize_input.py",
                "--input-base64",
                encoded,
                "--out",
                tmp_path / "params-base64.json",
                cwd=tmp,
            )
            self.assertEqual(base64_result.returncode, 0, base64_result.stderr)
            decoded = json.loads((tmp_path / "params-base64.json").read_text())
            self.assertEqual(decoded["target_organism"], payload["target_organism"])

    def test_parse_results_preserves_replica_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            replicas = [
                {
                    "replica": 1,
                    "seed": 10,
                    "poses": [{"rank": 1, "affinity_kcal_mol": -6.0}],
                    "top_affinity_kcal_mol": -6.0,
                },
                {
                    "replica": 2,
                    "seed": 11,
                    "poses": [{"rank": 1, "affinity_kcal_mol": -7.0}],
                    "top_affinity_kcal_mol": -7.0,
                },
                {
                    "replica": 3,
                    "seed": 12,
                    "poses": [{"rank": 1, "affinity_kcal_mol": -5.0}],
                    "top_affinity_kcal_mol": -5.0,
                },
            ]
            params = {"replica_runs": replicas, "selected_replica": 2}
            (tmp_path / "params.json").write_text(json.dumps(params))
            result = run_helper("parse_results.py", "--params", "params.json", cwd=tmp)
            self.assertEqual(result.returncode, 0, result.stderr)
            parsed = json.loads((tmp_path / "params.json").read_text())
            docking = parsed["docking_results"]
            self.assertEqual(docking["selected_replica"], 2)
            self.assertEqual(docking["top_affinity_kcal_mol"], -7.0)
            self.assertEqual(docking["replicas"], replicas)
            self.assertEqual(docking["replica_top_score_mean"], -6.0)
            self.assertEqual(docking["replica_top_score_range"], 2.0)

    def test_interactions_use_first_model_and_heavy_atoms(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            receptor = pdbqt_atom(1, "C1", "ALA", "A", 1, 0, 0, 0, "C")
            receptor += pdbqt_atom(2, "H1", "ALA", "A", 1, 0, 0, 0, "HD")
            ligand = "MODEL 1\n"
            ligand += pdbqt_atom(1, "C1", "LIG", "Z", 1, 1, 0, 0, "C")
            ligand += pdbqt_atom(2, "H1", "LIG", "Z", 1, 1, 0, 0, "HD")
            ligand += "ENDMDL\nMODEL 2\n"
            ligand += pdbqt_atom(3, "C2", "LIG", "Z", 1, 1, 0, 0, "C")
            ligand += "ENDMDL\n"
            (tmp_path / "receptor.pdbqt").write_text(receptor)
            (tmp_path / "docked.pdbqt").write_text(ligand)
            (tmp_path / "params.json").write_text(json.dumps({"cutoff": 4.5}))
            result = run_helper(
                "analyze_interactions.py", "--params", "params.json", cwd=tmp
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((tmp_path / "interactions.json").read_text())
            self.assertEqual(summary["ligand_docked_atom_count"], 1)
            self.assertEqual(summary["contact_residue_count"], 1)
            self.assertEqual(summary["residues"][0]["contact_count"], 1)

    def test_reports_and_package_are_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            params = {
                "run_id": "run-a",
                "ligand_id": "a",
                "ligand_name": "a",
                "smiles": "CC",
                "pdb_id": "1ABC",
                "chain_ids": ["A"],
                "heterogen_policy": "remove_all",
                "ph": 7.4,
                "docking_results": {
                    "selected_replica": 1,
                    "top_affinity_kcal_mol": -6.0,
                    "replica_top_scores_kcal_mol": [-6.0],
                    "poses": [{"rank": 1, "affinity_kcal_mol": -6.0}],
                },
                "protein_stats": {"source": {"sha256": "abc"}},
            }
            (tmp_path / "params.json").write_text(json.dumps(params))
            md = tmp_path / "handoff"
            md.mkdir()
            (md / "best_pose_ligand.sdf").write_text("pose\n$$$$\n")
            (md / "best_pose_ligand.pdb").write_text(
                "HETATM    1  C1  LIG Z   1      10.000  11.000  12.000  1.00 20.00           C\n"
            )
            (md / "prepared_receptor.pdb").write_text(
                "ATOM      1  CA  ALA A   1       9.000  10.000  11.000  1.00 20.00           C\n"
            )
            (md / "manifest.json").write_text("{}")
            result = run_helper(
                "generate_reports.py",
                "--params",
                "params.json",
                "--out-dir",
                tmp_path / "report_files",
                "--templates",
                NF / "templates",
                "--handoff-dir",
                md,
                cwd=tmp,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            packaged = run_helper(
                "build_report_package.py",
                "--params",
                "params.json",
                "--report-dir",
                tmp_path / "report_files",
                "--out-dir",
                tmp_path / "package",
                cwd=tmp,
            )
            self.assertEqual(packaged.returncode, 0, packaged.stderr)
            report_dir = tmp_path / "package"
            for name in [
                *REPORT_PAGES,
                "all_best_poses.sdf",
                "receptor_overlay_reference.pdb",
                "ligand_overlay.json",
            ]:
                self.assertTrue((report_dir / name).is_file(), name)
            viewer = (report_dir / "09_visualization.html").read_text()
            self.assertIn("receptorAtoms = [", viewer)
            self.assertIn("ligandAtoms = [", viewer)
            self.assertIn("viewerReady", viewer)
            self.assertNotIn("fetch(", viewer)
            self.assertNotIn("<script src=", viewer)
            self.assertIn("<tbody><tr>", viewer)

if __name__ == "__main__":
    unittest.main()
