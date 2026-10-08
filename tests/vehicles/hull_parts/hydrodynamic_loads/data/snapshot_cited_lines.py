"""Freeze every MSS / generator line and table the hydrodynamics tests cite.

Writes ``cited_lines_snapshot.json`` beside this file: the stripped text of
each ``CITED_LINES`` entry and the parsed numbers of each ``TABLES`` entry of
``tests/vehicles/hull_parts/hydrodynamic_loads/test_hydrodynamics_block.py``, plus the git revisions
the text came from. The tests read the snapshot, so their gates run without
an MSS checkout (nothing relative to a machine); with ``MSS_DIR`` set, ``test_cited_lines_are_unchanged`` compares it with the live files.

Run from the repository root (written 2026-10-06):

    MSS_DIR=<an MSS checkout> python tests/vehicles/hull_parts/hydrodynamic_loads/data/snapshot_cited_lines.py

``$TESTS/...`` pins are files of this repository (the template generators)
and need no variable.

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
TEST_FILE = HERE.parent / "test_hydrodynamics_block.py"


def _root(variable):
    value = os.environ.get(variable)
    if not value:
        sys.exit(f"{variable} is not set")
    return Path(value).expanduser()


def _revision(path):
    out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True)
    dirty = subprocess.run(["git", "-C", str(path), "status", "--porcelain"], capture_output=True, text=True)
    return {"head": out.stdout.strip(), "clean": dirty.stdout.strip() == ""}


def main():
    spec = importlib.util.spec_from_file_location("hydro_tests", TEST_FILE)
    tests = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tests)
    variables = sorted({ref[1:].split("/", 1)[0].partition("@")[0] for ref, _ in tests.CITED_LINES} - {"TESTS"})
    roots = {v: _root(v) for v in variables}
    # "$VARIABLE/path" is the file in that checkout, "$VARIABLE@<rev>/path" the
    # file at that git revision of it; the test module's own reader does both.
    live = tests._live_text

    lines, tables = {}, {}
    for (ref, number), text in sorted(tests.CITED_LINES.items()):
        line = live(ref).splitlines()[number - 1].strip()
        if not line.startswith(text):
            sys.exit(f"{ref}:{number} no longer starts with {text!r}: {line!r}")
        lines.setdefault(ref, {})[str(number)] = line
    for ref, name in tests.TABLES:
        tables.setdefault(ref, {})[name] = tests._parse_table(live(ref), name).tolist()
    snapshot = {
        "made_by": "tests/vehicles/hull_parts/hydrodynamic_loads/data/snapshot_cited_lines.py",
        "revisions": {v: _revision(p) for v, p in roots.items()},
        "lines": lines,
        "tables": tables,
    }
    (HERE / "cited_lines_snapshot.json").write_text(json.dumps(snapshot, indent=1, sort_keys=True) + "\n")
    print(f"{sum(len(v) for v in lines.values())} lines, {len(tests.TABLES)} tables")


if __name__ == "__main__":
    main()
