"""n8n Code node: Normalize Input.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import json
    import os
    import re
    import shutil
    import tempfile
    import time

    j = _items[0]["json"]
    raw = str(j.get("data", "")).strip()
    if not raw:
        raise ValueError("Drive input is empty")
    name = str(j.get("source_file_name") or "drive-input.txt")
    ext = os.path.splitext(name.lower())[1]
    if ext not in (".txt", ".smi", ".smiles", ".json"):
        raise ValueError("Unsupported Drive input extension " + ext + "; use .txt, .smi, .smiles, or .json")
    is_json = raw.lstrip().startswith("{")
    if is_json:
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError("Invalid JSON input: " + str(exc))
        if not isinstance(data, dict):
            raise ValueError("JSON input must be one object")
    else:
        data = {}
        for line in raw.splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                data["smiles"] = line.split()[0]
                data["ligand_name"] = " ".join(line.split()[1:]) or "test_ligand"
                break
    smiles = str(data.get("smiles", "")).strip()
    if not smiles:
        raise ValueError("Missing required field: smiles")
    if len(smiles) > 1000:
        raise ValueError("SMILES exceeds the 1000-character safety limit")
    if "." in smiles and not bool(data.get("allow_multicomponent", False)):
        raise ValueError("Multi-component SMILES are rejected; provide one ligand or explicitly set allow_multicomponent=true")
    pdb_value = str(data.get("pdb_id") or "").strip().upper()
    if pdb_value == "AUTO":
        pdb_value = ""
    selection_mode = str(data.get("receptor_selection_mode") or ("provided" if pdb_value else "auto")).strip().lower()
    if selection_mode not in ("auto", "provided"):
        raise ValueError("receptor_selection_mode must be auto or provided")
    if selection_mode == "provided":
        if not re.match(r"^[0-9][A-Z0-9]{3}$", pdb_value):
            raise ValueError("Provided-receptor mode requires a four-character pdb_id such as 1AKI")
    elif pdb_value:
        raise ValueError("Set receptor_selection_mode=provided to force a pdb_id, or omit pdb_id for automatic selection")
    pdb_id = pdb_value if selection_mode == "provided" else ""
    def number(key, default, low, high):
        try:
            value = float(data.get(key, default))
        except (TypeError, ValueError):
            raise ValueError(key + " must be numeric")
        if value < low or value > high:
            raise ValueError(key + " must be between " + str(low) + " and " + str(high))
        return value
    def integer(key, default, low, high):
        value = number(key, default, low, high)
        if int(value) != value:
            raise ValueError(key + " must be an integer")
        return int(value)
    center_values = [data.get("center_x"), data.get("center_y"), data.get("center_z")]
    has_any_center = any(value is not None for value in center_values)
    if has_any_center and not all(value is not None for value in center_values):
        raise ValueError("Provide all of center_x, center_y, and center_z or omit all three")
    center = [number(key, 0, -10000, 10000) for key in ("center_x", "center_y", "center_z")] if has_any_center else [None, None, None]
    sizes = [number(key, 20, 8, 30) for key in ("size_x", "size_y", "size_z")]
    if sizes[0] * sizes[1] * sizes[2] > 27000:
        raise ValueError("Docking box volume exceeds the 27000 A^3 safety limit")
    chains_raw = data.get("chain_ids", [])
    if isinstance(chains_raw, str):
        chain_ids = [part.strip() for part in chains_raw.split(",") if part.strip()]
    elif isinstance(chains_raw, list):
        chain_ids = [str(part).strip() for part in chains_raw if str(part).strip()]
    else:
        raise ValueError("chain_ids must be a comma-separated string or a JSON array")
    heterogen_policy = str(data.get("heterogen_policy") or "remove_all").strip().lower()
    if heterogen_policy not in ("remove_all", "keep_water"):
        raise ValueError("heterogen_policy must be remove_all or keep_water")
    target_organism = str(data.get("target_organism") or "Homo sapiens").strip()
    if not target_organism or len(target_organism) > 100:
        raise ValueError("target_organism must be a non-empty organism name under 100 characters")
    run_root = "/md_project/data/runs"
    os.makedirs(run_root, exist_ok=True)
    now = time.time()
    for entry in os.listdir(run_root):
        path = os.path.join(run_root, entry)
        if entry.startswith("dock_") and os.path.isdir(path):
            try:
                if now - os.stat(path).st_mtime > 172800:
                    shutil.rmtree(path)
            except OSError:
                pass
    run_dir = tempfile.mkdtemp(prefix="dock_", dir=run_root)
    input_dir = os.path.join(run_dir, "input")
    output_dir = os.path.join(run_dir, "output")
    reports_dir = os.path.join(run_dir, "reports")
    for path in (input_dir, output_dir, reports_dir):
        os.makedirs(path, exist_ok=True)
    run_id = os.path.basename(run_dir)
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(data.get("ligand_name") or "ligand"))[:40].strip("_") or "ligand"
    result = {
        "run_id": run_id,
        "run_name": (pdb_id or "AUTO") + "_" + safe_name + "_" + run_id,
        "run_dir": run_dir,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "reports_dir": reports_dir,
        "input_profile": "json-production-config" if is_json else "plain-smiles-auto-selection",
        "smiles": smiles,
        "ligand_name": safe_name,
        "pdb_id": pdb_id,
        "chain_ids": chain_ids,
        "receptor_selection_mode": selection_mode,
        "target_organism": target_organism,
        "target_similarity_threshold": number("target_similarity_threshold", 70, 40, 100),
        "target_candidate_limit": integer("target_candidate_limit", 5, 1, 10),
        "heterogen_policy": heterogen_policy,
        "add_missing_residues": bool(data.get("add_missing_residues", False)),
        "cx": center[0], "cy": center[1], "cz": center[2],
        "sx": sizes[0], "sy": sizes[1], "sz": sizes[2],
        "ph": number("ph", 7.4, 4, 10),
        "exhaustiveness": integer("exhaustiveness", 8, 8, 64),
        "num_modes": integer("num_modes", 9, 1, 20),
        "energy_range": number("energy_range", 3, 1, 10),
        "seed": integer("seed", 20260824, 1, 2147483647),
        "cpu": integer("cpu", 1, 1, 8),
        "replicas": integer("replicas", 1, 1, 3),
        "timeout_seconds": integer("timeout_seconds", 900, 60, 1800),
        "cutoff": number("cutoff", 4.5, 2.5, 6),
        "source_file_id": j.get("source_file_id"),
        "source_file_name": name,
        "source_mime_type": j.get("source_mime_type"),
        "source_modified_time": j.get("source_modified_time"),
        "source_md5": j.get("source_md5"),
        "qc_flags": [],
    }
    if not has_any_center:
        result["qc_flags"].append("GRID_CENTER_TO_BE_SELECTED_OR_INFERRED")
    if selection_mode == "auto":
        result["qc_flags"].append("AUTOMATED_RECEPTOR_SELECTION_REQUIRES_REVIEW")
    if not is_json:
        result["qc_flags"].append("PLAIN_SMILES_DEFAULTS_USED")
    return [{"json": result}]
