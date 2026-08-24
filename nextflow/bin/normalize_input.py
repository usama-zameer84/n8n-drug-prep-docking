#!/usr/bin/env python3
"""Validate ligand parameters and write params.json."""

import argparse
import base64
import json
import re


def number(value, default, low, high, key):
    try:
        v = float(value if value not in (None, "", "null") else default)
    except (TypeError, ValueError):
        raise ValueError(f"{key} must be numeric")
    if v < low or v > high:
        raise ValueError(f"{key} must be between {low} and {high}")
    return v


def integer(value, default, low, high, key):
    v = number(value, default, low, high, key)
    if int(v) != v:
        raise ValueError(f"{key} must be an integer")
    return int(v)


def boolean(value):
    if isinstance(value, bool):
        return value
    return str(value).lower() in ("true", "1", "yes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-json")
    ap.add_argument("--input-base64")
    ap.add_argument("--ligand-id")
    ap.add_argument("--smiles")
    ap.add_argument("--receptor-selection-mode", default="auto")
    ap.add_argument("--pdb-id", default="")
    ap.add_argument("--chain-ids", default="")
    ap.add_argument("--center-x", default=None)
    ap.add_argument("--center-y", default=None)
    ap.add_argument("--center-z", default=None)
    ap.add_argument("--size-x", default=20)
    ap.add_argument("--size-y", default=20)
    ap.add_argument("--size-z", default=20)
    ap.add_argument("--ph", default=7.4)
    ap.add_argument("--exhaustiveness", default=8)
    ap.add_argument("--num-modes", default=9)
    ap.add_argument("--energy-range", default=3)
    ap.add_argument("--seed", default=20260824)
    ap.add_argument("--cpu", default=1)
    ap.add_argument("--replicas", default=1)
    ap.add_argument("--timeout-seconds", default=900)
    ap.add_argument("--cutoff", default=4.5)
    ap.add_argument("--target-organism", default="Homo sapiens")
    ap.add_argument("--target-similarity-threshold", default=70)
    ap.add_argument("--target-candidate-limit", default=5)
    ap.add_argument("--heterogen-policy", default="remove_all")
    ap.add_argument("--add-missing-residues", default="false")
    ap.add_argument("--allow-multicomponent", default="false")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.input_json and args.input_base64:
        raise SystemExit("Use only one of --input-json or --input-base64")
    if args.input_json:
        with open(args.input_json, encoding="utf-8") as fh:
            payload = json.load(fh)
    elif args.input_base64:
        try:
            payload = json.loads(
                base64.b64decode(args.input_base64, validate=True).decode("utf-8")
            )
        except (ValueError, UnicodeDecodeError, json.JSONDecodeError) as error:
            raise SystemExit(f"Invalid base64 input JSON: {error}") from error
    else:
        payload = None

    if payload is not None:
        if not isinstance(payload, dict):
            raise SystemExit("input JSON must contain an object")
        allowed = {action.dest for action in ap._actions}
        unknown = sorted(set(payload) - allowed)
        if unknown:
            raise SystemExit("Unsupported input field(s): " + ", ".join(unknown))
        for key, value in payload.items():
            setattr(args, key, value)

    if args.ligand_id is None or args.smiles is None:
        raise SystemExit("ligand_id and smiles are required")

    smiles = args.smiles.strip()
    if not smiles:
        raise SystemExit("Missing required field: smiles")
    if len(smiles) > 1000:
        raise SystemExit("SMILES exceeds the 1000-character safety limit")
    if "." in smiles and not boolean(args.allow_multicomponent):
        raise SystemExit(
            "Multi-component SMILES are rejected; set allow_multicomponent=true"
        )

    pdb_value = (args.pdb_id or "").strip().upper()
    if pdb_value == "AUTO":
        pdb_value = ""
    mode = (
        (args.receptor_selection_mode or ("provided" if pdb_value else "auto"))
        .strip()
        .lower()
    )
    if mode not in ("auto", "provided"):
        raise SystemExit("receptor_selection_mode must be auto or provided")
    if mode == "provided":
        if not re.match(r"^[0-9][A-Z0-9]{3}$", pdb_value):
            raise SystemExit(
                "Provided-receptor mode requires a four-character pdb_id such as 1AKI"
            )
    elif pdb_value:
        raise SystemExit(
            "Set receptor_selection_mode=provided to force a pdb_id, or omit pdb_id for automatic selection"
        )
    pdb_id = pdb_value if mode == "provided" else ""

    # Nextflow interpolates `params.center_x = null` as the literal string "null"; treat it as missing.
    centers = [
        None if c in (None, "", "null") else c
        for c in (args.center_x, args.center_y, args.center_z)
    ]
    has_any = any(c is not None for c in centers)
    if has_any and not all(c is not None for c in centers):
        raise SystemExit(
            "Provide all of center_x, center_y, center_z or omit all three"
        )
    center = (
        [
            number(c, 0, -10000, 10000, k)
            for c, k in zip(centers, ("center_x", "center_y", "center_z"))
        ]
        if has_any
        else [None, None, None]
    )
    sizes = [
        number(getattr(args, k), 20, 8, 30, k) for k in ("size_x", "size_y", "size_z")
    ]
    if sizes[0] * sizes[1] * sizes[2] > 27000:
        raise SystemExit("Docking box volume exceeds the 27000 A^3 safety limit")

    chains_raw = args.chain_ids
    if isinstance(chains_raw, str):
        chain_ids = (
            [p.strip() for p in chains_raw.split(",") if p.strip()]
            if chains_raw
            else []
        )
    else:
        chain_ids = list(chains_raw)
    heterogen_policy = (args.heterogen_policy or "remove_all").strip().lower()
    if heterogen_policy not in ("remove_all", "keep_water"):
        raise SystemExit("heterogen_policy must be remove_all or keep_water")
    organism = (args.target_organism or "Homo sapiens").strip()
    if not organism or len(organism) > 100:
        raise SystemExit(
            "target_organism must be a non-empty organism name under 100 characters"
        )

    qc_flags = []
    if not has_any:
        qc_flags.append("GRID_CENTER_TO_BE_SELECTED_OR_INFERRED")
    if mode == "auto":
        qc_flags.append("AUTOMATED_RECEPTOR_SELECTION_REQUIRES_REVIEW")

    safe_name = (
        re.sub(r"[^A-Za-z0-9_.-]+", "_", args.ligand_id)[:40].strip("_") or "ligand"
    )
    d = {
        "run_id": f"nextflow_{safe_name}",
        "run_name": (pdb_id or "AUTO") + "_" + safe_name,
        "ligand_id": args.ligand_id,
        "smiles": smiles,
        "ligand_name": safe_name,
        "pdb_id": pdb_id,
        "chain_ids": chain_ids,
        "receptor_selection_mode": mode,
        "target_organism": organism,
        "target_similarity_threshold": number(
            args.target_similarity_threshold, 70, 40, 100, "target_similarity_threshold"
        ),
        "target_candidate_limit": integer(
            args.target_candidate_limit, 5, 1, 10, "target_candidate_limit"
        ),
        "heterogen_policy": heterogen_policy,
        "add_missing_residues": boolean(args.add_missing_residues),
        "cx": center[0],
        "cy": center[1],
        "cz": center[2],
        "sx": sizes[0],
        "sy": sizes[1],
        "sz": sizes[2],
        "ph": number(args.ph, 7.4, 4, 10, "ph"),
        "exhaustiveness": integer(args.exhaustiveness, 8, 8, 64, "exhaustiveness"),
        "num_modes": integer(args.num_modes, 9, 1, 20, "num_modes"),
        "energy_range": number(args.energy_range, 3, 1, 10, "energy_range"),
        "seed": integer(args.seed, 20260824, 1, 2147483647, "seed"),
        "cpu": integer(args.cpu, 1, 1, 8, "cpu"),
        "replicas": integer(args.replicas, 1, 1, 3, "replicas"),
        "timeout_seconds": integer(
            args.timeout_seconds, 900, 60, 1800, "timeout_seconds"
        ),
        "cutoff": number(args.cutoff, 4.5, 2.5, 6, "cutoff"),
        "qc_flags": qc_flags,
    }
    with open(args.out, "w") as fh:
        json.dump(d, fh, indent=2)


if __name__ == "__main__":
    main()
