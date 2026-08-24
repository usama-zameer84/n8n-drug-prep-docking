#!/usr/bin/env python3
"""Extract n8n Code nodes into importable, reviewable source files."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKFLOW = ROOT / "workflow" / "drug-prep-docking.workflow.json"
MANIFEST_PATH = ROOT / "scripts" / "manifest.json"


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return slug or "code_node"


def code_nodes(workflow: dict) -> list[dict]:
    nodes = [
        node
        for node in workflow.get("nodes", [])
        if node.get("type") == "n8n-nodes-base.code"
    ]
    return sorted(nodes, key=lambda node: tuple(node.get("position", [0, 0])))


def source_for(node: dict) -> tuple[str, str, str]:
    parameters = node.get("parameters", {})
    if parameters.get("language") == "python" or "pythonCode" in parameters:
        return "python", "py", parameters.get("pythonCode", "")
    return "javascript", "js", parameters.get("jsCode", "")


def render_python(node_name: str, source: str) -> str:
    body = textwrap.indent(source.rstrip() + "\n", "    ")
    return (
        f'"""n8n Code node: {node_name}.\n\n'
        "Call run(_items) with the item structure supplied by the n8n Python runner.\n"
        '"""\n\n'
        "\n"
        "def run(_items):\n"
        f"{body}"
    )


def render_javascript(node_name: str, source: str) -> str:
    body = textwrap.indent(source.rstrip() + "\n", "  ")
    return (
        f"/** n8n Code node: {node_name}. */\n"
        "// n8n supplies globals such as $json, $items and $binary at runtime.\n\n"
        "async function run() {\n"
        f"{body}"
        "}\n\n"
        "module.exports = { run };\n"
    )


def expected_files(workflow_path: Path) -> tuple[dict[Path, str], dict]:
    workflow = json.loads(workflow_path.read_text(encoding="utf-8"))
    files: dict[Path, str] = {}
    entries = []
    for index, node in enumerate(code_nodes(workflow), 1):
        language, extension, source = source_for(node)
        relative = Path("scripts") / language / (
            f"{index:02d}_{slugify(node['name'])}.{extension}"
        )
        rendered = (
            render_python(node["name"], source)
            if language == "python"
            else render_javascript(node["name"], source)
        )
        files[ROOT / relative] = rendered
        entries.append(
            {
                "node": node["name"],
                "node_id": node.get("id"),
                "language": language,
                "file": relative.as_posix(),
                "source_sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
            }
        )
    manifest = {
        "workflow": workflow_path.relative_to(ROOT).as_posix(),
        "code_node_count": len(entries),
        "nodes": entries,
    }
    files[MANIFEST_PATH] = json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    return files, manifest


def check(files: dict[Path, str]) -> int:
    failures = []
    for path, expected in files.items():
        if not path.exists():
            failures.append(f"missing: {path.relative_to(ROOT)}")
        elif path.read_text(encoding="utf-8") != expected:
            failures.append(f"out of date: {path.relative_to(ROOT)}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"Verified {len(files) - 1} extracted Code-node files.")
    return 0


def write(files: dict[Path, str]) -> None:
    for path, content in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    print(f"Extracted {len(files) - 1} Code nodes.")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow", nargs="?", type=Path, default=DEFAULT_WORKFLOW)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    workflow_path = args.workflow.resolve()
    files, _ = expected_files(workflow_path)
    if args.check:
        return check(files)
    write(files)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
