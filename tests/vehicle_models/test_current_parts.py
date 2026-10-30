"""The named parts that take the ocean current in the form MSS's two vehicle models use: "REMUS 100 (MSS, ocean
current)" (``current_form = "yaw_rate_terms"``, ``remus100.m`` 118-123 @ 6e78e98) and "Otter (MSS, ocean current)"
(``current_form = "full_rotation_rate"``, ``otter_legacy.m`` 126-131 @ 6e78e98). Each is a copy of the still-water part with
that one selector changed: every parameter file equal but it, the same linked actuators; built by rpp's own
``ComponentContextBuilder`` the vehicle has the three current inputs before its actuator inputs.

Author:    Enio Krizman
Date:      2026-10-26
"""

import json
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc

LIBRARY = Path(__file__).resolve().parents[2]
PARTS = LIBRARY / ".rppws" / "parts"
PAIRS = [("SpheroidAuv", "REMUS 100 (MSS)", "REMUS 100 (MSS, ocean current)", "yaw_rate_terms",
          ["fin_angle_command"] * 4 + ["drive_command"]),
         ("Catamaran", "Otter (MSS)", "Otter (MSS, ocean current)", "full_rotation_rate",
          ["drive_command"] * 2)]
CURRENT = ["current_speed", "current_direction", "current_vertical_speed"]


def _folder(plugin, name):
    found = [d.parent for d in PARTS.glob("*/*/description.json")
             if json.loads(d.read_text()).get("PluginName") == f"more_dynamics::{plugin}"
             and json.loads(d.read_text()).get("Name") == name]
    assert len(found) == 1, (plugin, name, found)
    return found[0]


def _tree(plugin, name):
    """``({slot: parameter values}, linked part ids)`` of a named part and its own children."""
    from rpp_orchestrator.component_storage import ComponentParameterStore
    folder = _folder(plugin, name)
    record = json.loads((folder / "description.json").read_text())
    store = ComponentParameterStore()
    values, links = {"root": store.load(folder)}, []
    for slot, entry in record["Subcomponents"].items():
        for k, item in enumerate(entry if isinstance(entry, list) else [entry]):
            child = folder / "subcomponents" / item["Id"]
            if item["IsLinked"]:
                links.append(json.loads((child / "description.json").read_text())["LinkedComponentId"])
            else:
                values[slot] = store.load(child)
    return values, links


def _build(plugin, name):
    vc.plugins()
    try:
        from rpp_py.data_manager import DataManager
        from rpp_plugin_registrator.library_manager import LibraryManager
        from rpp_py.context_builder import ComponentContextBuilder
        manager = DataManager(library_manager=LibraryManager(), workspace_path=str(LIBRARY))
        context = ComponentContextBuilder(data_manager=manager).build_component_from_path(
            str(_folder(plugin, name)))
        context.initialize()
    except Exception as exc:  # more_dynamics is not registered here
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
    return context.get_instance()


@pytest.mark.parametrize("plugin, still, current, form, commands", PAIRS)
def test_the_part_differs_from_the_still_water_part_in_the_current_form_only(plugin, still, current, form, commands):
    a, a_links = _tree(plugin, still)
    b, b_links = _tree(plugin, current)
    assert a_links == b_links and a_links
    assert a.keys() == b.keys()
    different = {(slot, n) for slot in a for n in a[slot] if repr(a[slot][n]) != repr(b[slot][n])}
    assert different == {("root", "current_form")}, different
    assert a["root"]["current_form"] == "none" and b["root"]["current_form"] == form


@pytest.mark.parametrize("plugin, still, current, form, commands", PAIRS)
def test_the_part_builds_with_the_current_inputs_before_the_actuator_inputs(plugin, still, current, form, commands):
    vehicle = _build(plugin, current)
    assert [d.name for d in vehicle.graph().inputDescription] == CURRENT + commands
    assert [d.name for d in _build(plugin, still).graph().inputDescription] == commands


@pytest.mark.parametrize("plugin, still, current, form, commands", PAIRS)
def test_a_current_moves_the_part_and_zero_current_equals_the_still_water_part(plugin, still, current, form, commands):
    moving, still_part = _build(plugin, current), _build(plugin, still)
    x = np.array([1.2, 0.1, 0.05, 0.02, -0.01, 0.1, 0, 0, 0, 0.05, 0.03, 0.4])
    u = np.zeros(len(commands))
    with_zero = np.asarray(moving.step_vector(x, np.r_[0.0, 0.0, 0.0, u], 0.05))
    dry = np.asarray(still_part.step_vector(x, u, 0.05))
    assert np.max(np.abs(with_zero - dry)) <= 1e-12
    flowing = np.asarray(moving.step_vector(x, np.r_[0.5, 0.5, 0.1, u], 0.05))
    assert np.max(np.abs(flowing - dry)) > 1e-4


# ---- the unit-quaternion form of the kinematics (``kinematics_form``) -------------------------------------------------
QUATERNION_PART = "REMUS 100 (MSS, quaternion)"
BASE_PART = "REMUS 100 (MSS, ocean current)"
CURRENT_INPUTS = np.array([0.5, 0.52, 0.1])
COMMANDS = np.array([0.05, 0.05, -0.03, -0.03, 1300.0])


def _euler_state():
    return np.array([1.0, -2.0, 3.0, 0.1, 0.2, 0.5, 1.2, 0.1, 0.05, 0.02, -0.01, 0.1])


def _quaternion_state(euler_state):
    from more_transformations.more_casadi_transformations import EulerQuaternionTransforms
    q = np.asarray(EulerQuaternionTransforms.euler_to_quaternion(*euler_state[3:6])).ravel()
    return np.r_[euler_state[:3], q, euler_state[6:]]


def test_the_quaternion_part_differs_from_the_ocean_current_part_in_the_kinematics_form_only():
    a, a_links = _tree("SpheroidAuv", BASE_PART)
    b, b_links = _tree("SpheroidAuv", QUATERNION_PART)
    assert a_links == b_links and a_links and a.keys() == b.keys()
    different = {(slot, n) for slot in a for n in set(a[slot]) | set(b[slot]) if repr(a[slot].get(n)) != repr(b[slot].get(n))}
    assert different == {("root", "kinematics_form")}, different
    assert b["root"]["kinematics_form"] == "unit_quaternion" and a["root"].get("kinematics_form", "euler") == "euler"


def test_the_quaternion_vehicle_has_13_states_and_the_12_value_output_of_the_euler_one():
    from more_common.casadi_graph import RppCasadiGraph
    quaternion, euler = _build("SpheroidAuv", QUATERNION_PART), _build("SpheroidAuv", BASE_PART)
    gq, ge = RppCasadiGraph(quaternion.graph()), RppCasadiGraph(euler.graph())
    assert [d.size for d in quaternion.graph().stateDescription][-1] == 13
    assert [d.size for d in euler.graph().stateDescription][-1] == 12
    assert [d.name for d in quaternion.graph().inputDescription] == [d.name for d in euler.graph().inputDescription]
    assert [(d.name, d.size) for d in quaternion.graph().outputDescription] == \
        [(d.name, d.size) for d in euler.graph().outputDescription] == [("output", 12)]
    xe = _euler_state()
    u = np.r_[CURRENT_INPUTS, COMMANDS]
    states_q = _quaternion_state(xe)
    # fin and propulsor children carry their own states before the vessel's: take them from the vehicle's description
    pad_q = np.zeros(gq.num_states - 13) if hasattr(gq, "num_states") else np.zeros(0)
    pad_e = np.zeros(ge.num_states - 12) if hasattr(ge, "num_states") else np.zeros(0)
    out_q = np.asarray(gq.output(np.r_[pad_q, states_q], u)).ravel()
    out_e = np.asarray(ge.output(np.r_[pad_e, xe], u)).ravel()
    assert np.max(np.abs(out_q - out_e)) <= 1e-14


def test_the_quaternion_form_follows_the_euler_form_over_a_manoeuvre_and_keeps_the_quaternion_at_unit_length():
    from more_transformations.more_casadi_transformations import EulerQuaternionTransforms
    quaternion, euler = _build("SpheroidAuv", QUATERNION_PART), _build("SpheroidAuv", BASE_PART)
    xe = _euler_state()
    xq = _quaternion_state(xe)
    u = np.r_[CURRENT_INPUTS, COMMANDS]
    for _ in range(400):
        xe = np.asarray(euler.step_vector(xe, u, 0.05)).ravel()
        xq = np.asarray(quaternion.step_vector(xq, u, 0.05)).ravel()
    angles = np.array([float(a) for a in EulerQuaternionTransforms.quaternion_to_euler(xq[3:7])])
    assert abs(np.linalg.norm(xq[3:7]) - 1.0) <= 1e-15
    assert np.max(np.abs(np.r_[xq[:3], angles, xq[7:]] - xe)) <= 1e-4  # two integrations of one motion (RK4, 20 s)
    assert np.max(np.abs(np.r_[xq[:3], angles, xq[7:]] - xe)) > 0.0


def test_the_euler_form_is_the_default_of_the_vehicle_type_and_the_other_types_do_not_offer_the_option():
    """The parts that do not set ``kinematics_form`` carry the 12 Euler states; the plugin's default is ``euler``; only
    ``SpheroidAuv`` declares the option (MSS runs only the REMUS in both forms: ``otter_legacy.m`` 97 is 12-state)."""
    from more_dynamics.plugins.vehicle_models.spheroid_auv import SpheroidAuv
    from more_dynamics.plugins.vehicle_models.catamaran import Catamaran
    from more_dynamics.plugins.vehicle_models.monohull import Monohull
    from more_dynamics.plugins.vehicle_models.derivative_auv import DerivativeAuv
    assert [p.default_value for p in SpheroidAuv.PARAMETERS if p.name == "kinematics_form"] == ["euler"]
    for plugin in (Catamaran, Monohull, DerivativeAuv):
        assert "kinematics_form" not in [p.name for p in plugin.PARAMETERS], plugin.__name__
    assert [d.size for d in _build("SpheroidAuv", "REMUS 100 (MSS)").graph().stateDescription][-1] == 12


def test_an_unknown_kinematics_form_is_refused_naming_the_options():
    from more_dynamics.plugins.shared.vehicle_options import kinematics_form
    from more_dynamics.plugins.shared.vehicle_graph import CompositionError

    class Context:
        def get_parameter(self, name):
            return "rotation_vector"
    with pytest.raises(CompositionError, match="kinematics_form must be one of"):
        kinematics_form(Context(), "SpheroidAuv", CompositionError)


def test_the_quaternion_is_renormalised_in_the_vehicles_step_and_would_drift_without_it():
    """RK4 alone leaves the unit length by about 1e-7 over 1400 s of a manoeuvre (measured 8e-8 on the part under the
    commands of MSS's run); the vehicle's own step keeps it to rounding."""
    from more_common.casadi_graph import RppCasadiGraph
    vehicle = _build("SpheroidAuv", QUATERNION_PART)
    derivative = RppCasadiGraph(vehicle.graph()).step
    h = 0.05
    plain = _quaternion_state(_euler_state())
    kept = plain.copy()
    plain_drift = kept_drift = 0.0
    for k in range(4000):
        u = np.r_[CURRENT_INPUTS, 0.05 * np.sin(0.01 * k), 0.05 * np.sin(0.01 * k), -0.03, -0.03, 1300.0]
        f = lambda y: np.asarray(derivative(y, u)).ravel()
        k1 = f(plain); k2 = f(plain + h / 2 * k1); k3 = f(plain + h / 2 * k2); k4 = f(plain + h * k3)
        plain = plain + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        kept = np.asarray(vehicle.step_vector(kept, u, h)).ravel()
        plain_drift = max(plain_drift, abs(np.linalg.norm(plain[3:7]) - 1.0))
        kept_drift = max(kept_drift, abs(np.linalg.norm(kept[3:7]) - 1.0))
    assert kept_drift <= 1e-15
    assert plain_drift > 1e-10
