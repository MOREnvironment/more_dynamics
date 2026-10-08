"""Vehicle gate helpers: build a ``MarineCraft6DOF`` vehicle from an ``.rppws``
tree with rpp's own ``ComponentContextBuilder`` and read its payload.

The vehicle class is data: one ``MarineCraft6DOF`` plugin receives its parts
by slot, and a vehicle is a tree of part folders with their ``parameters.py``.
REMUS 100 and the Otter differ only in their trees. The gate trees live in
``tests/data/vehicles/.rppws`` (written by ``tests/data/vehicles/make_trees.py``);
a perturbed or refused variant is written into a temporary folder from the
same generator and built in the same way. Reference numbers live under
``tests/data/vehicles``.

The trees need rpp's registry with ``more_dynamics`` registered
(``rpp library register ./more_dynamics --link``); without it the tree tests
skip with that message and the model-layer gates run as usual.

Author:    Enio Krizman
Date:      2026-10-08
"""

import sys
from pathlib import Path

import casadi as ca
import numpy as np
import pytest

DATA = Path(__file__).resolve().parents[1] / "data" / "vehicles"
GATE_TREES = DATA / ".rppws"
GATE_SCRIPT = GATE_TREES / "script_descriptions" / "vehicles.json"
SLOTS = ("site", "current", "hull_form", "rigid_body", "added_mass",
         "added_mass_coriolis", "restoring", "hydrodynamic_loads", "force_producers")
CURRENT_NAMES = ("current_speed", "current_direction", "current_vertical_speed")
CURRENT_INPUTS = tuple(f"current.{n}" for n in CURRENT_NAMES)
REGISTER_MESSAGE = ("rpp's registry has no registered more_dynamics: run "
                    "`rpp library register ./more_dynamics --link` (workspace/ as the working folder)")

sys.path.insert(0, str(DATA))
import make_trees  # noqa: E402  (the generator that wrote the gate trees)


def builder():
    """rpp's ``ComponentContextBuilder`` over the native registry, or skip."""
    try:
        from rpp_py.data_manager import DataManager
        from rpp_plugin_registrator.library_manager import LibraryManager
        manager = LibraryManager()  # puts ~/.rpp/interfaces/python on sys.path
        data_manager = DataManager(library_manager=manager)  # before the generated modules are imported
        from rpp_py.context_builder import ComponentContextBuilder
        import rpp_plugin_types.more_dynamics  # noqa: F401
        return ComponentContextBuilder(data_manager=data_manager)
    except Exception as exc:  # no rpp, no home, or more_dynamics not registered
        pytest.skip(f"{REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")


def build(configuration, script=GATE_SCRIPT):
    """The vehicle plugin of one configuration of a script description (built and initialised)."""
    context = builder().build_script_from_description_path(str(script), configuration=configuration)
    context.initialize()
    return context.get_component("vessels")[0]


def build_variant(tmp_path, tree, name="variant"):
    """Write ``tree`` (a node of ``make_trees``) as the only configuration of a workspace in ``tmp_path`` and
    build it."""
    workspace = Path(tmp_path) / ".rppws"
    make_trees.write_workspace(workspace, {name: tree}, "tests/vehicles/vehicle_contract.py", {"variant": [name]})
    return build(name, workspace / "script_descriptions" / "variant.json")


def graph_of(vessel):
    from more_common.casadi_graph import RppCasadiGraph
    return RppCasadiGraph(vessel.graph())


def input_names(vessel):
    return [d.name for d in vessel.getInputDescriptions()]


def input_size(vessel):
    """Number of stacked input values of the vehicle."""
    return sum(d.size for d in vessel.getInputDescriptions())


def entries(descriptions, vector):
    """``{name: values}`` of a stacked vector in description order."""
    out, k = {}, 0
    for d in descriptions:
        out[d.name] = np.asarray(vector[k:k + d.size]).ravel()
        k += d.size
    return out


def output_entries(vessel, x, u):
    """Every named output entry of the vehicle (default output, then the declared diagnostic outputs)."""
    graph = graph_of(vessel)
    values = np.asarray(graph.output(x, u)).ravel()
    return entries(vessel.graph().outputDescription, values)


def matrix(vector):
    """A payload matrix: its column-major ``vec`` back to 6x6."""
    return np.asarray(vector).reshape(6, 6, order="F")


def dynamics_map(vessel, n):
    """``dynamics`` mapped over ``n`` columns: ``f(x (12 x n), u (m x n)) -> xdot (12 x n)``."""
    graph = graph_of(vessel)
    return graph.step.map(n)


def mss_to_ours(x):
    """MSS state [nu; eta] to ours [eta; nu]."""
    x = np.asarray(x)
    return np.concatenate([x[..., 6:12], x[..., 0:6]], axis=-1)


ours_to_mss = mss_to_ours


def tree_node(name):
    """The generator's node of a gate tree, for a variant."""
    return make_trees.gate_trees()[name]
