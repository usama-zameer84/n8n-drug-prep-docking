from __future__ import annotations

import ast
import json
import shutil
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
N8N = ROOT / "n8n"
WORKFLOW_PATH = N8N / "workflow" / "drug-prep-docking.workflow.json"


class RepositoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workflow = json.loads(WORKFLOW_PATH.read_text(encoding="utf-8"))
        cls.nodes = cls.workflow["nodes"]
        cls.by_name = {node["name"]: node for node in cls.nodes}

    def test_workflow_structure(self):
        self.assertFalse(self.workflow["active"])
        self.assertEqual(len(self.nodes), 29)
        self.assertEqual(
            self.workflow["connections"]["Normalize Input"]["main"][0][0]["node"],
            "Auto Select Receptor",
        )
        self.assertEqual(
            self.workflow["connections"]["Auto Select Receptor"]["main"][0][0]["node"],
            "Generate 3D Structure",
        )

    def test_public_drive_configuration(self):
        self.assertTrue(all("credentials" not in node for node in self.nodes))
        trigger = self.by_name["SMILES Drop (Drive Trigger)"]
        report_folder = self.by_name["Create Run Folder (Drive)"]
        self.assertEqual(
            trigger["parameters"]["folderToWatch"]["value"],
            "YOUR_GOOGLE_DRIVE_INPUT_FOLDER_ID",
        )
        self.assertEqual(
            report_folder["parameters"]["folderId"]["value"],
            "YOUR_GOOGLE_DRIVE_REPORTS_PARENT_FOLDER_ID",
        )

    def test_code_node_exports_match(self):
        result = subprocess.run(
            [sys.executable, "tools/extract_code_nodes.py", "--check"],
            cwd=N8N,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        manifest = json.loads((N8N / "scripts" / "manifest.json").read_text())
        self.assertEqual(manifest["code_node_count"], 17)

    def test_python_exports_parse(self):
        for path in sorted((N8N / "scripts" / "python").glob("*.py")):
            with self.subTest(path=path.name):
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    @unittest.skipUnless(shutil.which("node"), "Node.js is not installed")
    def test_javascript_exports_parse(self):
        for path in sorted((N8N / "scripts" / "javascript").glob("*.js")):
            with self.subTest(path=path.name):
                result = subprocess.run(
                    ["node", "--check", str(path)],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_report_and_handoff_artifacts(self):
        report_source = self.by_name["Build Scientific Report Package"]["parameters"][
            "pythonCode"
        ]
        for name in (
            "03a_receptor_selection.html",
            "09_visualization.html",
            "receptor_overlay_reference.pdb",
            "all_best_poses.sdf",
            "ligand_overlay.json",
        ):
            self.assertIn(name, report_source)

        handoff_source = self.by_name["Build MD Handoff Bundle"]["parameters"][
            "pythonCode"
        ]
        for name in (
            "best_pose_ligand.sdf",
            "complex_best_pose.pdb",
            "atom_mapping.json",
            "provenance.json",
            "manifest.json",
        ):
            self.assertIn(name, handoff_source)

    def test_repository_hygiene(self):
        blocked_names = {".DS_Store", ".env", "credentials.json"}
        blocked_directories = {
            ".nextflow",
            ".ruff_cache",
            "__pycache__",
            "node_modules",
            "results",
            "work",
        }
        result = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
            cwd=ROOT,
            capture_output=True,
            check=True,
        )
        publication_candidates = [
            Path(raw.decode()) for raw in result.stdout.split(b"\0") if raw
        ]
        for path in publication_candidates:
            self.assertNotIn(path.name, blocked_names)
            self.assertFalse(
                blocked_directories.intersection(path.parts),
                f"generated path is publishable: {path}",
            )

        for ignored in (
            ".DS_Store",
            ".env",
            ".idea/workspace.xml",
            ".nextflow/history",
            ".nextflow.log.1",
            "report.html",
            "node_modules/package/index.js",
            "work/task/file",
            "results/run/file",
            ".ruff_cache/cache",
        ):
            check = subprocess.run(
                ["git", "check-ignore", "--quiet", ignored], cwd=ROOT
            )
            self.assertEqual(check.returncode, 0, f"expected Git to ignore {ignored}")

    def test_manifest_paths_resolve_under_n8n(self):
        manifest = json.loads((N8N / "scripts" / "manifest.json").read_text())
        self.assertEqual(manifest["code_node_count"], len(manifest["nodes"]))
        for node in manifest["nodes"]:
            with self.subTest(node=node["node"]):
                self.assertTrue((N8N / node["file"]).is_file())


if __name__ == "__main__":
    unittest.main()
