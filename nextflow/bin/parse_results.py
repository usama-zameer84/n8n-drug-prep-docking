#!/usr/bin/env python3
"""Create the stable docking-results contract from Vina replica metadata."""

import argparse
import json
import math


def build_results(data):
    replica_runs = data.get("replica_runs") or []
    if not replica_runs:
        raise ValueError("No replica metadata from docking stage")
    selected_replica = data.get("selected_replica")
    selected = next(
        (run for run in replica_runs if run.get("replica") == selected_replica), None
    )
    if selected is None or not selected.get("poses"):
        raise ValueError("Selected docking replica is missing or has no poses")
    top_scores = [float(run["top_affinity_kcal_mol"]) for run in replica_runs]
    mean = sum(top_scores) / len(top_scores)
    variance = sum((value - mean) ** 2 for value in top_scores) / len(top_scores)
    return {
        "status": "COMPLETED",
        "score_semantics": "AutoDock Vina scoring values in kcal/mol; use for pose ranking, not as experimental affinity or MD free energy.",
        "selected_replica": selected_replica,
        "top_affinity_kcal_mol": selected["top_affinity_kcal_mol"],
        "best_affinity_kcal_mol": selected["top_affinity_kcal_mol"],
        "poses": selected["poses"],
        "num_poses": len(selected["poses"]),
        "replicas": replica_runs,
        "replica_top_scores_kcal_mol": top_scores,
        "replica_top_score_mean": round(mean, 3),
        "replica_top_score_stddev": round(math.sqrt(variance), 3),
        "replica_top_score_range": round(max(top_scores) - min(top_scores), 3),
        "convergence_note": "Replica score spread is descriptive only; structural clustering and experimental validation are still required.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True)
    parser.add_argument("--out", default="docking_results.json")
    args = parser.parse_args()
    with open(args.params, encoding="utf-8") as fh:
        data = json.load(fh)
    data["docking_results"] = build_results(data)
    with open(args.params, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(data["docking_results"], fh, indent=2)


if __name__ == "__main__":
    main()
