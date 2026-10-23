"""The shared repository's jet-nozzle plugin run live, as an extra check beside the frozen rows.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

The frozen rows (``data/jet_nozzle_reference.csv``) are the reference the
propulsor gates read; this one module imports his plugin, only when it is
importable, and skips with a reason otherwise (his parts may be deleted later,
by agreement with their author: our parts and tests keep working without
them). It checks two things: that the frozen rows are still what his plugin
gives today (so a change in his file is seen), and that our his-forms
composition equals his live plugin on fresh seeded states and inputs.

Author:    Enio Krizman
Date:      2026-10-10
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from propulsor_contract import (EXACT_TOLERANCE, jet_forms_his, jet_parts, jet_reference, jet_values_his, max_diff,
                                propulsor_from_forms, rng, run)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vehicle_models"))
import vehicle_contract as vc  # noqa: E402

MODULE = "more_dynamics.plugins.force_producers.jet_nozzle"
HIS_NAMES = {"position": "location", "nozzle_velocity": "nozzle_velocity_radpersec"}


def _his_graph(numbers):
    vc.plugins()                                         # registers rpp's generated types (skips without rpp)
    if importlib.util.find_spec(MODULE) is None:
        pytest.skip(f"{MODULE} is not present: the shared repository's jet nozzle has been removed; the frozen rows stand")
    from more_common.casadi_graph import RppCasadiGraph
    cls = importlib.import_module(MODULE).JetNozzle
    context = vc.context_of(cls, {HIS_NAMES.get(k, k): v for k, v in numbers.items()})
    context.initialize()
    return RppCasadiGraph(context.get_instance().graph())


@pytest.mark.parametrize("name", sorted(jet_parts()))
def test_his_live_plugin_still_gives_the_frozen_rows(name):
    graph = _his_graph(jet_parts()[name])
    ref = jet_reference()[name]
    worst = 0.0
    for state, command, output, dot in list(zip(ref.state, ref.command, ref.output, ref.dot))[::5]:
        worst = max(worst, max_diff(np.array(graph.output(state, command)).reshape(-1), output),
                    max_diff(np.array(graph.step(state, command)).reshape(-1), dot))
    assert worst == 0.0, worst                           # the same function, bit for bit


@pytest.mark.parametrize("name", sorted(jet_parts()))
def test_our_his_forms_equal_his_live_plugin_on_fresh_states(name):
    numbers = jet_parts()[name]
    graph = _his_graph(numbers)
    ours = propulsor_from_forms(jet_forms_his(), jet_values_his(numbers))
    g = rng()
    worst = 0.0
    for _ in range(300):
        state = [g.uniform(-1.5, 1.5), g.uniform(numbers["min_angle"] - 0.3, numbers["max_angle"] + 0.3)]
        command = [g.uniform(-1.5, 1.5), g.uniform(-1.5, 1.5)]
        out = run(ours, command, np.zeros(6), 1025.0, state=state)
        worst = max(worst, max_diff(out["tau"], np.array(graph.output(state, command)).reshape(-1)),
                    max_diff(out["state_dot"], np.array(graph.step(state, command)).reshape(-1)))
    assert worst <= EXACT_TOLERANCE, worst
