"""n8n Code node: Parse Results.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import json
    import math
    import os

    d = _items[0]["json"]
    out_dir = d["output_dir"]
    selected = [item for item in d["replica_runs"] if item["replica"] == d["selected_replica"]][0]
    top_scores = [item["top_affinity_kcal_mol"] for item in d["replica_runs"]]
    mean = sum(top_scores) / len(top_scores)
    variance = sum((value - mean) ** 2 for value in top_scores) / len(top_scores)
    results = {
        "status": "COMPLETED",
        "score_semantics": "AutoDock Vina scoring values in kcal/mol; use for pose ranking, not as experimental affinity or MD free energy.",
        "selected_replica": d["selected_replica"],
        "top_affinity_kcal_mol": selected["top_affinity_kcal_mol"],
        "poses": selected["poses"],
        "num_poses": len(selected["poses"]),
        "replica_top_scores_kcal_mol": top_scores,
        "replica_top_score_mean": round(mean, 3),
        "replica_top_score_stddev": round(math.sqrt(variance), 3),
        "replica_top_score_range": round(max(top_scores) - min(top_scores), 3),
        "convergence_note": "Replica score spread is descriptive only; structural clustering and experimental validation are still required.",
    }
    with open(os.path.join(out_dir, "docking_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    out = dict(d)
    out["docking_results"] = results
    return [{"json": out}]
