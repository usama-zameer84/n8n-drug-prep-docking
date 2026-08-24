"""n8n Code node: Cleanup Run Workspace.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import os
    import shutil

    d = _items[0]["json"]
    run_dir = str(d.get("run_dir", ""))
    root = "/md_project/data/runs"
    cleaned = False
    if run_dir and os.path.dirname(run_dir) == root and os.path.basename(run_dir).startswith("dock_") and os.path.isdir(run_dir):
        shutil.rmtree(run_dir)
        cleaned = True
    out = dict(d)
    out["workspace_cleaned"] = cleaned
    out.pop("run_dir", None)
    out.pop("input_dir", None)
    out.pop("output_dir", None)
    out.pop("reports_dir", None)
    return [{"json": out}]
