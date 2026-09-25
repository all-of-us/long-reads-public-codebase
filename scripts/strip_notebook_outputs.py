#!/usr/bin/env python3
"""Strip Jupyter notebook output cells (and leftover execution state).

Default: rewrite notebooks in place. CI uses --check and fails if any
tracked notebook still has outputs.

Examples:
  python3 scripts/strip_notebook_outputs.py
  python3 scripts/strip_notebook_outputs.py path/to/notebook.ipynb
  python3 scripts/strip_notebook_outputs.py --check
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Optional


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def tracked_notebooks(root: Path) -> List[Path]:
    raw = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "-z", "*.ipynb"],
    )
    notebooks: List[Path] = []
    for entry in raw.split(b"\0"):
        if not entry:
            continue
        path = root / Path(entry.decode())
        if path.suffix == ".ipynb" and ".ipynb_checkpoints" not in path.parts:
            notebooks.append(path)
    return notebooks


def iter_notebooks(paths: Iterable[Path]) -> List[Path]:
    found: List[Path] = []
    for path in paths:
        if path.is_file() and path.suffix == ".ipynb":
            if ".ipynb_checkpoints" not in path.parts:
                found.append(path)
            continue
        if path.is_dir():
            for nb in sorted(path.rglob("*.ipynb")):
                if ".ipynb_checkpoints" not in nb.parts:
                    found.append(nb)
    return found


def strip_notebook(nb: dict) -> bool:
    """Clear result cells. Return True if anything changed."""
    changed = False

    metadata = nb.setdefault("metadata", {})
    if metadata.pop("widgets", None) is not None:
        changed = True

    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue

        outputs = cell.get("outputs")
        if outputs:
            cell["outputs"] = []
            changed = True
        elif "outputs" not in cell:
            cell["outputs"] = []
            changed = True

        if cell.get("execution_count") is not None:
            cell["execution_count"] = None
            changed = True

        cell_meta = cell.setdefault("metadata", {})
        for key in ("collapsed", "scrolled", "ExecuteTime", "execution", "papermill"):
            if key in cell_meta:
                del cell_meta[key]
                changed = True

    return changed


def dumps_notebook(nb: dict) -> str:
    return json.dumps(nb, indent=1, ensure_ascii=False) + "\n"


def load_notebook(text: str) -> dict:
    nb = json.loads(text)
    if not isinstance(nb, dict) or "cells" not in nb:
        raise ValueError("not a Jupyter notebook")
    return nb


def check_notebook(nb: dict, path: str) -> List[str]:
    problems: List[str] = []
    if (nb.get("metadata") or {}).get("widgets"):
        problems.append(f"{path}: notebook metadata contains widgets")
    for i, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        if cell.get("outputs"):
            problems.append(f"{path}: code cell {i} has outputs")
        if cell.get("execution_count") is not None:
            problems.append(f"{path}: code cell {i} has execution_count")
    return problems


def write_if_changed(path: Path, stripped: dict) -> bool:
    if not strip_notebook(stripped):
        return False
    path.write_text(dumps_notebook(stripped), encoding="utf-8")
    return True


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="Notebooks or directories (default: repository root)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit 1 if any notebook still has outputs; do not write",
    )
    parser.add_argument(
        "--stdin",
        action="store_true",
        help="Read one notebook from stdin and write the stripped notebook to stdout",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if args.stdin:
        original = sys.stdin.read()
        nb = load_notebook(original)
        strip_notebook(nb)
        sys.stdout.write(dumps_notebook(nb))
        return 0

    if args.paths:
        notebooks = iter_notebooks(args.paths)
    else:
        notebooks = tracked_notebooks(repo_root())
    if not notebooks:
        print("No notebooks found.", file=sys.stderr)
        return 1

    if args.check:
        problems: List[str] = []
        for path in notebooks:
            nb = load_notebook(path.read_text(encoding="utf-8"))
            problems.extend(check_notebook(nb, str(path)))
        if problems:
            print("Notebooks still contain result cells:", file=sys.stderr)
            for problem in problems:
                print(f"  {problem}", file=sys.stderr)
            print(
                "Run: python3 scripts/strip_notebook_outputs.py",
                file=sys.stderr,
            )
            return 1
        print(f"OK: {len(notebooks)} notebook(s) have no output cells.")
        return 0

    rewritten = 0
    for path in notebooks:
        nb = load_notebook(path.read_text(encoding="utf-8"))
        if write_if_changed(path, nb):
            print(f"stripped {path}")
            rewritten += 1
    if rewritten == 0:
        print(f"Already clean: {len(notebooks)} notebook(s).")
    else:
        print(f"Stripped outputs from {rewritten} notebook(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
