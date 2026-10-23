"""Record the output and the dynamics of the shared repository's jet-nozzle plugin
(``more_dynamics/plugins/force_producers/jet_nozzle.py``) on seeded states and
inputs, once, for the three parameter sets of ``jet_nozzle_test_parts.json``.

The plugin is built with rpp's own ``ComponentContext`` (as the library's other
plugin tests build theirs), its payload read through ``RppCasadiGraph``, and
every row written with 17 significant digits. Run from the library root:

    python tests/force_producers/propulsor/data/generate_jet_nozzle_reference.py

It writes ``jet_nozzle_reference.csv`` under ``OUT_DIR`` (default
``<tempdir>/more_mss_references/propulsor/``), never into this folder, so a
re-run cannot overwrite the frozen file; compare with ``cmp`` against the
frozen copy. The frozen file is the only thing the propulsor gates read: the
propulsor parts import nothing from the jet-nozzle plugin.

Columns: set, thrust_state, angle_state, desired_thrust,
desired_normalized_nozzle_angle, output_0 .. output_5 (the plugin's
``generated_force``), dot_0, dot_1 (``thrust``, ``nozzle_angle`` derivatives).

Author:    Enio Krizman
Date:      2026-10-10
"""

import csv
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIBRARY = HERE.parents[3]
SEED = 20261010
PER_SET = 200
HIS_NAMES = {"position": "location", "nozzle_velocity": "nozzle_velocity_radpersec"}


def main():
    sys.path.insert(0, str(LIBRARY / "tests" / "vehicle_models"))
    from rpp_py.data_manager import DataManager
    DataManager()  # before the generated plugin types are imported
    from more_dynamics.plugins.force_producers.jet_nozzle import JetNozzle
    from more_common.casadi_graph import RppCasadiGraph
    import vehicle_contract as vc

    sets = json.loads((HERE / "jet_nozzle_test_parts.json").read_text())["sets"]
    rng = np.random.default_rng(SEED)
    out_dir = Path(os.environ.get("OUT_DIR") or Path(tempfile.gettempdir()) / "more_mss_references") / "propulsor"
    out_dir.mkdir(parents=True, exist_ok=True)
    header = ["set", "thrust_state", "angle_state", "desired_thrust", "desired_normalized_nozzle_angle",
              *[f"output_{i}" for i in range(6)], "dot_0", "dot_1"]
    rows = []
    for name, numbers in sets.items():
        params = {HIS_NAMES.get(k, k): v for k, v in numbers.items()}
        context = vc.context_of(JetNozzle, params)
        context.initialize()
        graph = RppCasadiGraph(context.get_instance().graph())
        low, high = numbers["min_angle"], numbers["max_angle"]
        states = [[0.0, 0.0], [1.0, high], [-1.0, low]]
        inputs = [[0.0, 0.0], [1.0, 1.0], [-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0]]
        cases = [(s, i) for s in states for i in inputs]
        while len(cases) < PER_SET:  # beyond the clamps on purpose
            cases.append(([rng.uniform(-1.3, 1.3), rng.uniform(low - 0.2, high + 0.2)],
                          [rng.uniform(-1.3, 1.3), rng.uniform(-1.3, 1.3)]))
        for state, inp in cases:
            out = np.array(graph.output(state, inp), dtype=float).reshape(-1)
            dot = np.array(graph.step(state, inp), dtype=float).reshape(-1)
            rows.append([name, *state, *inp, *out, *dot])
    path = out_dir / "jet_nozzle_reference.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow([row[0], *[repr(float(v)) for v in row[1:]]])
    print(f"{len(rows)} rows -> {path}")


if __name__ == "__main__":
    main()
