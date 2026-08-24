#!/usr/bin/env python3
"""Run AutoDock Vina with the configured replicas and search box."""

import argparse
import json
import os
import re
import shutil
import subprocess


SCORE_RE = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$"
)


def coords(path):
    values = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    values.append(
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    )
                except ValueError:
                    pass
    if not values:
        raise ValueError("Prepared receptor has no readable coordinates")
    return values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    args = ap.parse_args()
    d = json.load(open(args.params))

    receptor_path = "receptor.pdbqt"
    ligand_path = "ligand.pdbqt"
    c = coords(receptor_path)
    mins = [min(row[i] for row in c) for i in range(3)]
    maxs = [max(row[i] for row in c) for i in range(3)]
    auto = d["cx"] is None
    if not auto:
        center = [d["cx"], d["cy"], d["cz"]]
        source = "user_supplied"
    elif d.get("auto_grid_center"):
        center = [float(v) for v in d["auto_grid_center"]]
        source = d.get("grid_center_source") or "co_crystallized_ligand"
    else:
        center = [sum(row[i] for row in c) / len(c) for i in range(3)]
        source = "receptor_centroid_fallback"
    size = [d["sx"], d["sy"], d["sz"]]
    overlaps = [
        center[i] + size[i] / 2 >= mins[i] and center[i] - size[i] / 2 <= maxs[i]
        for i in range(3)
    ]
    if not all(overlaps):
        raise SystemExit("Docking box does not overlap the receptor on all three axes")

    vina_path = os.environ.get("VINA_PATH", "vina")
    if not shutil.which(vina_path) and not os.path.exists(vina_path):
        raise SystemExit(
            "AutoDock Vina is not installed (expected `vina` on PATH or $VINA_PATH)"
        )
    vproc = subprocess.run(
        [vina_path, "--version"], capture_output=True, text=True, timeout=30
    )
    vina_version = ((vproc.stdout or "") + " " + (vproc.stderr or "")).strip()

    replica_runs = []
    for replica in range(d["replicas"]):
        seed = int(d["seed"]) + replica
        pose = f"docked_replica_{replica + 1}.pdbqt"
        log = f"vina_replica_{replica + 1}.log"
        cmd = [
            vina_path,
            "--receptor",
            receptor_path,
            "--ligand",
            ligand_path,
            "--center_x",
            str(center[0]),
            "--center_y",
            str(center[1]),
            "--center_z",
            str(center[2]),
            "--size_x",
            str(size[0]),
            "--size_y",
            str(size[1]),
            "--size_z",
            str(size[2]),
            "--exhaustiveness",
            str(d["exhaustiveness"]),
            "--num_modes",
            str(d["num_modes"]),
            "--energy_range",
            str(d["energy_range"]),
            "--seed",
            str(seed),
            "--cpu",
            str(d["cpu"]),
            "--out",
            pose,
        ]
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=d["timeout_seconds"]
            )
        except subprocess.TimeoutExpired:
            raise SystemExit(f"Vina replica {replica + 1} exceeded timeout_seconds")
        log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
        open(log, "w").write(log_text)
        if proc.returncode != 0:
            raise SystemExit(f"Vina replica {replica + 1} failed: " + log_text[-1600:])
        scores = []
        for line in log_text.splitlines():
            m = SCORE_RE.match(line)
            if m:
                scores.append(
                    {
                        "rank": int(m.group(1)),
                        "affinity_kcal_mol": float(m.group(2)),
                        "rmsd_lb_A": float(m.group(3)),
                        "rmsd_ub_A": float(m.group(4)),
                    }
                )
        if not scores:
            raise SystemExit(f"No Vina score table found for replica {replica + 1}")
        replica_runs.append(
            {
                "replica": replica + 1,
                "seed": seed,
                "pose_path": pose,
                "log_path": log,
                "poses": scores,
                "top_affinity_kcal_mol": scores[0]["affinity_kcal_mol"],
            }
        )

    best = min(replica_runs, key=lambda r: r["top_affinity_kcal_mol"])
    shutil.copyfile(best["pose_path"], "docked.pdbqt")
    combined = "\n\n".join(
        f"=== REPLICA {r['replica']} SEED {r['seed']} ===\n"
        + open(r["log_path"]).read()
        for r in replica_runs
    )
    open("docking.log", "w").write(combined)

    d["grid_center"] = [round(v, 3) for v in center]
    d["grid_size"] = size
    d["grid_auto_centered"] = auto
    d["grid_center_source"] = source
    d["receptor_bounds"] = {
        "min": [round(v, 3) for v in mins],
        "max": [round(v, 3) for v in maxs],
    }
    d["vina_version"] = vina_version
    d["replica_runs"] = replica_runs
    d["selected_replica"] = best["replica"]
    d["docking_results"] = {
        "best_affinity_kcal_mol": best["top_affinity_kcal_mol"],
        "replicas": replica_runs,
        "grid_center": d["grid_center"],
        "grid_size": size,
        "grid_center_source": source,
    }
    if source == "receptor_centroid_fallback":
        d["qc_flags"] = sorted(
            set(d.get("qc_flags", []) + ["AUTO_CENTER_NOT_SCIENTIFICALLY_VALIDATED"])
        )
    json.dump(d, open(args.params, "w"), indent=2)


if __name__ == "__main__":
    main()
