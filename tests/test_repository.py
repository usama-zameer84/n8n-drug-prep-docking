from __future__ import annotations

import ast
import json
import shutil
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_PATH = ROOT / "workflow" / "drug-prep-docking.workflow.json"


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
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        manifest = json.loads((ROOT / "scripts" / "manifest.json").read_text())
        self.assertEqual(manifest["code_node_count"], 17)

    def test_python_exports_parse(self):
        for path in sorted((ROOT / "scripts" / "python").glob("*.py")):
            with self.subTest(path=path.name):
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    @unittest.skipUnless(shutil.which("node"), "Node.js is not installed")
    def test_javascript_exports_parse(self):
        for path in sorted((ROOT / "scripts" / "javascript").glob("*.js")):
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
        blocked_names = {".mcp.json", ".env", "credentials.json"}
        for path in ROOT.rglob("*"):
            if ".git" in path.parts or path.is_dir():
                continue
            self.assertNotIn(path.name, blocked_names)


if __name__ == "__main__":
    unittest.main()
