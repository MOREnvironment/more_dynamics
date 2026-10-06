"""Freeze every MSS / generator line and table the hydrodynamics tests cite.

Writes ``cited_lines_snapshot.json`` beside this file: the stripped text of
each ``CITED_LINES`` entry and the parsed numbers of each ``TABLES`` entry of
``tests/hydrodynamics/test_hydrodynamics_block.py``, plus the git revisions
the text came from. The tests read the snapshot, so their gates run without
an MSS or ``more_generic_models`` checkout (agents-more rule 9); with the
variables set, ``test_cited_lines_are_unchanged`` compares it with the live files.

Run from the repository root (job U3a, 2026-10-06):

    MSS_DIR=<an MSS checkout> MORE_GENERIC_MODELS_DIR=<more_generic_models repo> \\
        python tests/data/hydrodynamics/snapshot_cited_lines.py

MSS lines are MIT (T. I. Fossen, MSS); they are quoted here one line at a time
as citations.
"""

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEST_FILE = HERE.parents[1] / "hydrodynamics" / "test_hydrodynamics_block.py"


def _root(variable):
    value = os.environ.get(variable)
    if not value:
        sys.exit(f"{variable} is not set (agents-more rule 9)")
    return Path(value).expanduser()


def _revision(path):
    out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True)
    dirty = subprocess.run(["git", "-C", str(path), "status", "--porcelain"], capture_output=True, text=True)
    return {"head": out.stdout.strip(), "clean": dirty.stdout.strip() == ""}


def main():
    spec = importlib.util.spec_from_file_location("hydro_tests", TEST_FILE)
    tests = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tests)
    roots = {v: _root(v) for v in ("MSS_DIR", "MORE_GENERIC_MODELS_DIR")}

    def live(ref):
        variable, rel = ref[1:].split("/", 1)
        return (roots[variable] / rel).read_text()

    lines, tables = {}, {}
    for (ref, number), text in sorted(tests.CITED_LINES.items()):
        line = live(ref).splitlines()[number - 1].strip()
        if not line.startswith(text):
            sys.exit(f"{ref}:{number} no longer starts with {text!r}: {line!r}")
        lines.setdefault(ref, {})[str(number)] = line
    for ref, name in tests.TABLES:
        tables.setdefault(ref, {})[name] = tests._parse_table(live(ref), name).tolist()
    snapshot = {
        "made_by": "tests/data/hydrodynamics/snapshot_cited_lines.py (job U3a)",
        "revisions": {v: _revision(p) for v, p in roots.items()},
        "lines": lines,
        "tables": tables,
    }
    (HERE / "cited_lines_snapshot.json").write_text(json.dumps(snapshot, indent=1, sort_keys=True) + "\n")
    print(f"{sum(len(v) for v in lines.values())} lines, {len(tests.TABLES)} tables")


if __name__ == "__main__":
    main()
