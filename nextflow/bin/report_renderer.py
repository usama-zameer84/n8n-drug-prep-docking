"""Render report pages from files in the report template directory."""

import html
import json
from pathlib import Path
from string import Template


def escape(value):
    return html.escape(str(value) if value is not None else "—")


class ReportRenderer:
    def __init__(self, template_dir):
        self.template_dir = Path(template_dir)

    def _render(self, name, **values):
        source = (self.template_dir / name).read_text(encoding="utf-8")
        return Template(source).substitute(values)

    def _data_page(self, title, value):
        return self._render(
            "data.html",
            title=escape(title),
            payload=escape(json.dumps(value, indent=2, sort_keys=True)),
        )

    def _index(self, data):
        docking = data.get("docking_results", {}) or {}
        ro5 = data.get("ro5", {}) or {}
        replica_scores = docking.get("replica_top_scores_kcal_mol", [])
        rows = [
            ("SMILES", data.get("smiles")),
            ("Ligand", data.get("ligand_name")),
            ("Receptor (PDB)", data.get("pdb_id")),
            ("Selection mode", data.get("receptor_selection_mode")),
            ("Best affinity (kcal/mol)", docking.get("best_affinity_kcal_mol")),
            ("Replicas run", len(replica_scores) or len(docking.get("replicas", []))),
            ("Grid center source", data.get("grid_center_source")),
            ("Ro5 druglike", ro5.get("druglike")),
            ("Ro5 violations", ro5.get("violations")),
        ]
        table_rows = "".join(
            f"<tr><th>{escape(label)}</th><td>{escape(value)}</td></tr>"
            for label, value in rows
        )
        qc_flags = ", ".join(data.get("qc_flags", [])) or "none"
        return self._render(
            "index.html",
            ligand_name=escape(data.get("ligand_name")),
            run_id=escape(data.get("run_id")),
            table_rows=table_rows,
            qc_flags=escape(qc_flags),
        )

    def _receptor_selection(self, data):
        selection = data.get("receptor_selection", {}) or {}
        selected = selection.get("selected", {}) or {}
        candidates = selection.get("candidate_targets", []) or []
        preflight = selection.get("preparation_preflight", []) or []
        candidate_rows = "".join(
            "<tr>"
            f"<td>{escape(item.get('target_chembl_id'))}</td>"
            f"<td>{escape(item.get('target_name'))}</td>"
            f"<td>{escape(item.get('uniprot_accession'))}</td>"
            f"<td>{escape(item.get('evidence_score'))}</td>"
            "</tr>"
            for item in candidates
        )
        preflight_rows = "".join(
            "<tr>"
            f"<td>{escape(item.get('pdb_id'))}</td>"
            f"<td>{escape(','.join(item.get('chain_ids') or []))}</td>"
            f"<td>{escape(item.get('passed'))}</td>"
            f"<td>{escape(item.get('reason', ''))}</td>"
            "</tr>"
            for item in preflight
        )
        return self._render(
            "receptor_selection.html",
            ligand_name=escape(data.get("ligand_name")),
            pdb_id=escape(selected.get("pdb_id")),
            chains=escape(",".join(selected.get("chain_ids") or [])),
            selection_score=escape(selected.get("selection_score")),
            candidate_rows=candidate_rows,
            preflight_rows=preflight_rows,
            method=escape(selection.get("method")),
        )

    @staticmethod
    def _read_pdb_atoms(source, alpha_carbons_only=False):
        atoms = []
        for line in Path(source).read_text(encoding="utf-8").splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            if alpha_carbons_only and atom_name != "CA":
                continue
            try:
                x = round(float(line[30:38]), 3)
                y = round(float(line[38:46]), 3)
                z = round(float(line[46:54]), 3)
            except ValueError:
                continue
            element = line[76:78].strip() or atom_name[:1]
            atoms.append(
                {
                    "x": x,
                    "y": y,
                    "z": z,
                    "element": element.upper(),
                    "chain": line[21:22].strip(),
                    "residue": line[22:26].strip(),
                }
            )
        if not atoms:
            raise ValueError(f"No readable atoms in {source}")
        return atoms

    def _visualization(self, data, structure_dir):
        docking = data.get("docking_results", {}) or {}
        selected_replica = docking.get("selected_replica")
        rows = "".join(
            "<tr>"
            f"<td>{escape(selected_replica)}</td>"
            f"<td>{escape(pose.get('rank'))}</td>"
            f"<td>{escape(pose.get('affinity_kcal_mol'))}</td>"
            f"<td>{escape(pose.get('rmsd_lb_A'))}</td>"
            f"<td>{escape(pose.get('rmsd_ub_A'))}</td>"
            "</tr>"
            for pose in docking.get("poses", [])[:20]
        )
        return self._render(
            "visualization.html",
            ligand_name=escape(data.get("ligand_name")),
            pose_rows=rows,
            receptor_atoms=json.dumps(
                self._read_pdb_atoms(
                    Path(structure_dir) / "prepared_receptor.pdb",
                    alpha_carbons_only=True,
                ),
                separators=(",", ":"),
            ),
            ligand_atoms=json.dumps(
                self._read_pdb_atoms(Path(structure_dir) / "best_pose_ligand.pdb"),
                separators=(",", ":"),
            ),
        )

    def render_all(self, data, structure_dir):
        return {
            "index.html": self._index(data),
            "01_input_provenance.html": self._data_page(
                "Input provenance",
                {
                    "run_id": data.get("run_id"),
                    "ligand_id": data.get("ligand_id"),
                    "input_smiles": data.get("smiles"),
                    "canonical_smiles": data.get("canonical_smiles"),
                    "inchi_key": data.get("inchi_key"),
                    "receptor_source": (data.get("protein_stats") or {}).get("source"),
                },
            ),
            "02_ligand_preparation.html": self._data_page(
                "Ligand preparation",
                {
                    "conformer_generation": data.get("conformer_generation"),
                    "ligand_prep": data.get("ligand_prep"),
                    "rule_of_five": data.get("ro5"),
                },
            ),
            "03_receptor_preparation.html": self._data_page(
                "Receptor preparation", data.get("protein_stats")
            ),
            "03a_receptor_selection.html": self._receptor_selection(data),
            "04_docking_results.html": self._data_page(
                "Docking results",
                {
                    "vina_version": data.get("vina_version"),
                    "grid_center": data.get("grid_center"),
                    "grid_size": data.get("grid_size"),
                    "grid_center_source": data.get("grid_center_source"),
                    "docking_results": data.get("docking_results"),
                },
            ),
            "05_distance_contacts.html": self._data_page(
                "Distance contacts", data.get("interactions_summary")
            ),
            "06_methods_qc.html": self._data_page(
                "Methods and QC",
                {
                    "qc_flags": data.get("qc_flags"),
                    "protein_versions": (data.get("protein_stats") or {}).get("versions"),
                    "reproducibility": {
                        "seed": data.get("seed"),
                        "replicas": data.get("replicas"),
                        "ph": data.get("ph"),
                    },
                },
            ),
            "07_md_handoff_readiness.html": self._data_page(
                "MD handoff readiness", data.get("md_handoff")
            ),
            "08_raw_data.html": self._data_page("Raw run data", data),
            "09_visualization.html": self._visualization(data, structure_dir),
        }
