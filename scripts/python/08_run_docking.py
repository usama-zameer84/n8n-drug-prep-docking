"""n8n Code node: Run Docking.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import json
    import math
    import os
    import re
    import shutil
    import subprocess

    d = _items[0]["json"]
    out_dir = d["output_dir"]
    receptor_path = os.path.join(out_dir, "receptor.pdbqt")
    ligand_path = os.path.join(out_dir, "ligand.pdbqt")

    def coordinates(path):
        values = []
        with open(path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        values.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    except ValueError:
                        pass
        if not values:
            raise ValueError("Prepared receptor has no readable coordinates")
        return values

    coords = coordinates(receptor_path)
    mins = [min(row[i] for row in coords) for i in range(3)]
    maxs = [max(row[i] for row in coords) for i in range(3)]
    auto = d["cx"] is None
    if not auto:
        center = [d["cx"], d["cy"], d["cz"]]
        center_source = "user_supplied"
    elif d.get("auto_grid_center"):
        center = [float(value) for value in d["auto_grid_center"]]
        center_source = d.get("grid_center_source") or "co_crystallized_ligand"
    else:
        center = [sum(row[i] for row in coords) / len(coords) for i in range(3)]
        center_source = "receptor_centroid_fallback"
    size = [d["sx"], d["sy"], d["sz"]]
    overlaps = [
        center[i] + size[i] / 2 >= mins[i] and center[i] - size[i] / 2 <= maxs[i]
        for i in range(3)
    ]
    if not all(overlaps):
        raise ValueError("Docking box does not overlap the receptor on all three axes; verify center and size")

    vina_path = "/usr/local/bin/vina"
    if not os.path.exists(vina_path):
        raise ValueError("AutoDock Vina is not installed at /usr/local/bin/vina")
    version_proc = subprocess.run([vina_path, "--version"], capture_output=True, text=True, timeout=30)
    vina_version = ((version_proc.stdout or "") + " " + (version_proc.stderr or "")).strip()

    score_pattern = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$")
    replica_runs = []
    for replica in range(d["replicas"]):
        replica_seed = int(d["seed"]) + replica
        pose_path = os.path.join(out_dir, "docked_replica_" + str(replica + 1) + ".pdbqt")
        log_path = os.path.join(out_dir, "vina_replica_" + str(replica + 1) + ".log")
        cmd = [
            vina_path,
            "--receptor", receptor_path,
            "--ligand", ligand_path,
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
            "--exhaustiveness", str(d["exhaustiveness"]),
            "--num_modes", str(d["num_modes"]),
            "--energy_range", str(d["energy_range"]),
            "--seed", str(replica_seed),
            "--cpu", str(d["cpu"]),
            "--out", pose_path,
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=d["timeout_seconds"])
        except subprocess.TimeoutExpired:
            raise ValueError("Vina replica " + str(replica + 1) + " exceeded timeout_seconds")
        log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
        with open(log_path, "w") as fh:
            fh.write(log_text)
        if proc.returncode != 0:
            raise ValueError("Vina replica " + str(replica + 1) + " failed: " + log_text[-1600:])
        scores = []
        for line in log_text.splitlines():
            match = score_pattern.match(line)
            if match:
                scores.append({
                    "rank": int(match.group(1)),
                    "affinity_kcal_mol": float(match.group(2)),
                    "rmsd_lb_A": float(match.group(3)),
                    "rmsd_ub_A": float(match.group(4)),
                })
        if not scores:
            raise ValueError("No Vina score table found for replica " + str(replica + 1))
        replica_runs.append({
            "replica": replica + 1,
            "seed": replica_seed,
            "pose_path": pose_path,
            "log_path": log_path,
            "poses": scores,
            "top_affinity_kcal_mol": scores[0]["affinity_kcal_mol"],
        })

    best = min(replica_runs, key=lambda item: item["top_affinity_kcal_mol"])
    shutil.copyfile(best["pose_path"], os.path.join(out_dir, "docked.pdbqt"))
    combined_log = "\n\n".join("=== REPLICA " + str(item["replica"]) + " SEED " + str(item["seed"]) + " ===\n" + open(item["log_path"]).read() for item in replica_runs)
    with open(os.path.join(out_dir, "docking.log"), "w") as fh:
        fh.write(combined_log)

    out = dict(d)
    out["grid_center"] = [round(value, 3) for value in center]
    out["grid_size"] = size
    out["grid_auto_centered"] = auto
    out["grid_center_source"] = center_source
    out["receptor_bounds"] = {"min": [round(v, 3) for v in mins], "max": [round(v, 3) for v in maxs]}
    out["vina_version"] = vina_version
    out["replica_runs"] = replica_runs
    out["selected_replica"] = best["replica"]
    if center_source == "receptor_centroid_fallback":
        out["qc_flags"] = sorted(set(list(out["qc_flags"]) + ["AUTO_CENTER_NOT_SCIENTIFICALLY_VALIDATED"]))
    return [{"json": out}]
