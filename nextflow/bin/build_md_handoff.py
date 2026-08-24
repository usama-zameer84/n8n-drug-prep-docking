#!/usr/bin/env python3
"""Build the structural MD handoff and its checksum manifest."""

import argparse
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    d = json.load(open(args.params))
    os.makedirs(args.out_dir, exist_ok=True)
    from rdkit import Chem

    required_inputs = [
        "receptor_source.cif",
        "prepared_receptor.pdb",
        "receptor.pdbqt",
        "receptor.json",
        "receptor_prep.log",
        "protein_stats.json",
        "ligand_input.sdf",
        "ligand.pdbqt",
        "ligand_2d.svg",
        "docked.pdbqt",
        "docking.log",
        "docking_results.json",
        "interactions.json",
    ]
    missing = [name for name in required_inputs if not os.path.isfile(name)]
    if missing:
        raise SystemExit("MD handoff inputs are missing: " + ", ".join(missing))
    for name in required_inputs:
        shutil.copyfile(name, os.path.join(args.out_dir, name))
    for pattern in ("docked_replica_*.pdbqt", "vina_replica_*.log"):
        for source in sorted(glob.glob(pattern)):
            shutil.copyfile(
                source, os.path.join(args.out_dir, os.path.basename(source))
            )

    # Meeko mk_export: docked.pdbqt -> docked_poses.sdf (chemistry-preserving)
    docked_sdf = os.path.join(args.out_dir, "docked_poses.sdf")
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "meeko.cli.mk_export",
            "docked.pdbqt",
            "--write_sdf",
            docked_sdf,
        ],
        capture_output=True,
        text=True,
        timeout=180,
    )
    open(os.path.join(args.out_dir, "meeko_export.log"), "w").write(
        (proc.stdout or "") + "\n" + (proc.stderr or "")
    )
    if proc.returncode != 0 or not os.path.exists(docked_sdf):
        raise SystemExit("Meeko mk_export failed")
    poses = [m for m in Chem.SDMolSupplier(docked_sdf, removeHs=False) if m is not None]
    if not poses:
        raise SystemExit("Meeko exported no ligand poses")
    best = poses[0]
    docked_heavy = best.GetNumHeavyAtoms()
    expected_heavy = int(d.get("num_heavy_atoms", docked_heavy))
    if docked_heavy != expected_heavy:
        raise SystemExit(
            f"Heavy-atom QC failed: input {expected_heavy} != docked {docked_heavy}"
        )

    # atom mapping + PDB residue info
    element_counts = {}
    atom_mapping = []
    for index, atom in enumerate(best.GetAtoms()):
        symbol = atom.GetSymbol()
        element_counts[symbol] = element_counts.get(symbol, 0) + 1
        atom_name = (symbol + str(element_counts[symbol]))[:4]
        info = Chem.AtomPDBResidueInfo()
        info.SetName(atom_name.rjust(4))
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId("Z")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
        atom_mapping.append(
            {
                "docked_sdf_index_0_based": index,
                "pdb_atom_name": atom_name,
                "element": symbol,
                "atomic_number": atom.GetAtomicNum(),
                "formal_charge": atom.GetFormalCharge(),
            }
        )
    best_sdf = os.path.join(args.out_dir, "best_pose_ligand.sdf")
    w = Chem.SDWriter(best_sdf)
    w.write(best)
    w.close()
    best_pdb = os.path.join(args.out_dir, "best_pose_ligand.pdb")
    Chem.MolToPDBFile(best, best_pdb)
    json.dump(
        {
            "schema": "docked SDF atom index -> ligand PDB atom name",
            "atom_mapping": atom_mapping,
            "canonical_isomeric_smiles": d.get("canonical_smiles"),
            "inchi_key": d.get("inchi_key"),
        },
        open(os.path.join(args.out_dir, "atom_mapping.json"), "w"),
        indent=2,
    )

    # complex_best_pose.pdb: prepared receptor + best ligand pose
    complex_path = os.path.join(args.out_dir, "complex_best_pose.pdb")
    serial, rec_n, lig_n = 1, 0, 0
    with open(complex_path, "w") as out:
        for line in open("prepared_receptor.pdb"):
            if line.startswith(("ATOM", "HETATM")):
                out.write(
                    line[:6] + str(serial).rjust(5) + line[11:].rstrip("\n") + "\n"
                )
                serial += 1
                rec_n += 1
        out.write("TER\n")
        for line in open(best_pdb):
            if line.startswith(("ATOM", "HETATM")):
                out.write(
                    "HETATM" + str(serial).rjust(5) + line[11:].rstrip("\n") + "\n"
                )
                serial += 1
                lig_n += 1
        out.write("END\n")

    readme = """# MD structural handoff

Status: not topology-ready.

This folder contains a chemistry-preserving, coordinate-level handoff from rigid-receptor
AutoDock Vina docking. Use docked_poses.sdf or best_pose_ligand.sdf as the ligand chemistry
authority, prepared_receptor.pdb as the prepared protein coordinates, and complex_best_pose.pdb
only as a combined starting coordinate model.

Before any real MD: inspect the pose, choose protein/ligand force fields, generate ligand
parameters/charges, reconcile protonation/atom names with atom_mapping.json, define
termini/disulfides/cofactors/metals/waters, build topology, solvate, ionize, minimize,
equilibrate (NVT/NPT), run replicated production MD with convergence analysis.

Docking scores are ranking signals, not experimental affinities and not MD free energies.
"""
    open(os.path.join(args.out_dir, "README_MD_HANDOFF.md"), "w").write(readme)

    provenance = {
        "run_id": d.get("run_id"),
        "ligand": {
            "input_smiles": d["smiles"],
            "canonical_isomeric_smiles": d.get("canonical_smiles"),
            "inchi_key": d.get("inchi_key"),
            "preparation": d.get("ligand_prep"),
        },
        "receptor": d.get("protein_stats"),
        "receptor_selection": d.get("receptor_selection", {}),
        "docking": {
            "vina_version": d.get("vina_version"),
            "grid_center_A": d.get("grid_center"),
            "grid_size_A": d.get("grid_size"),
            "results": d.get("docking_results", {}),
        },
        "contacts": d.get("interactions_summary"),
        "qc_flags": d.get("qc_flags", []),
        "md_status": "STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY",
    }
    json.dump(
        provenance, open(os.path.join(args.out_dir, "provenance.json"), "w"), indent=2
    )

    files = {
        os.path.basename(path): path
        for path in glob.glob(os.path.join(args.out_dir, "*"))
        if os.path.isfile(path) and os.path.basename(path) != "manifest.json"
    }
    manifest = {
        "schema_version": "2.0",
        "status": "not topology-ready",
        "chemistry_authority": "best_pose_ligand.sdf",
        "coordinate_complex": "complex_best_pose.pdb",
        "input_atom_counts": {
            "all_atoms": d.get("num_atoms"),
            "heavy_atoms": d.get("num_heavy_atoms"),
        },
        "exported_atom_counts": {
            "all_atoms": best.GetNumAtoms(),
            "heavy_atoms": docked_heavy,
        },
        "complex_atom_counts": {
            "receptor": rec_n,
            "ligand": lig_n,
            "total": rec_n + lig_n,
        },
        "files": {
            n: {"sha256": sha256(p), "bytes": os.path.getsize(p)}
            for n, p in files.items()
        },
    }
    json.dump(
        manifest, open(os.path.join(args.out_dir, "manifest.json"), "w"), indent=2
    )

    d["md_handoff"] = {
        "status": "STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY",
        "manifest": manifest,
        "chemistry_qc_passed": True,
    }
    json.dump(d, open(args.params, "w"), indent=2)


if __name__ == "__main__":
    main()
