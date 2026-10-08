"""Assembly gates of ``MarineCraft6DOF``: couplings by name, the command list,
open parameters, diagnostic outputs, refusals, ``step`` and Luka's own parts
in the same vehicle.

The trees are built by rpp's ``ComponentContextBuilder``; the refusals that no
real part can provoke (a mis-sized input, a quantity produced twice, a state
outside a force producer) use parts written in the test with the same payload
shape and an initialised ``MarineCraft6DOF`` over a stand-in context that only
answers ``get_component`` and ``get_parameter``.

Author:    Enio Krizman
Date:      2026-10-08
"""

import copy
import csv
from pathlib import Path

import casadi as ca
import numpy as np
import pytest

import vehicle_contract as vc
from vehicle_contract import CURRENT_INPUTS, DATA, SLOTS, make_trees


def _vehicle_class():
    vc.builder()  # skips without a registered more_dynamics
    from more_dynamics.plugins.vehicles.marine_craft_6dof import CompositionError, MarineCraft6DOF
    return MarineCraft6DOF, CompositionError


# ---------------------------------------------------------------------------
# Couplings by name, the command list, diagnostic outputs
# ---------------------------------------------------------------------------

def test_slots_are_the_nine_physical_slots_in_order():
    for tree in (make_trees.torpedo(), make_trees.catamaran()):
        assert tuple(tree["children"]) == SLOTS


def test_every_input_is_connected_by_name_to_the_slot_that_produces_it():
    vessel = vc.build("remus100_hull")
    wiring = vessel.wiring
    expected = {
        "hydrodynamic_loads.0.added_mass_matrix": "slot added_mass (LambSpheroid)",
        "hydrodynamic_loads.0.rigid_body_mass_matrix": "slot rigid_body (HomogeneousSpheroid)",
        "hydrodynamic_loads.0.weight": "slot restoring (SubmergedNeutral)",
        "added_mass.water_density": "slot site (SiteAtLatitude)",
        "hydrodynamic_loads.2.section_beam": "slot hull_form (ProlateSpheroidMainDimensions)",
        "added_mass_coriolis.relative_velocity": "vehicle",
        "force_producers.0.desired_wrench": "command",
        "current.current_speed": "open input",
    }
    for key, origin in expected.items():
        assert wiring.get(key) == origin, (key, wiring.get(key))


@pytest.mark.parametrize("tree, commands", [
    ("remus100_hull", ["desired_wrench"]),
    ("remus100_full", ["fin_deflection_command", "shaft_speed_command"]),
    ("otter", []),
    ("otter_jet_nozzle", ["desired_thrust", "desired_normalized_nozzle_angle"]),
])
def test_command_list_of_each_vehicle(tree, commands):
    """A force-producer input that is neither a vehicle quantity nor an earlier part's output is a command, in
    list order, after the open inputs the tree lists."""
    vessel = vc.build(tree)
    names = vc.input_names(vessel)
    opened = list(CURRENT_INPUTS) if tree.startswith("remus") else []
    assert names == opened + commands


def test_diagnostic_outputs_follow_the_twelve_states():
    vessel = vc.build("remus100_hull")
    names = [d.name for d in vessel.graph().outputDescription]
    assert names == ["output", "mass_matrix", "relative_velocity"]
    sizes = [d.size for d in vessel.graph().outputDescription]
    assert sizes == [12, 36, 6]
    out = vc.output_entries(vessel, np.zeros(12), np.zeros(9))
    reference = DATA / "remus100" / "remus100_mass_matrix_consistent.csv"
    with open(reference, newline="") as handle:
        rows = list(csv.reader(handle))
    table = np.asarray(rows[1:], float)
    assert np.abs(vc.matrix(out["mass_matrix"]) - table[:, :6]).max() <= 1e-9
    assert np.array_equal(out["relative_velocity"], np.zeros(6))


def test_the_default_output_is_the_twelve_states_of_hull_vessel():
    vessel = vc.build("otter")
    assert [d.name for d in vessel.graph().outputDescription][0] == "output"
    assert vessel.graph().outputDescription[0].size == 12
    states = vessel.graph().stateDescription
    assert [(d.name, d.size) for d in states] == [("state", 12)]


# ---------------------------------------------------------------------------
# Open parameters
# ---------------------------------------------------------------------------

def test_gradient_through_an_open_current_parameter_equals_a_central_difference():
    vessel = vc.build("remus100_hull")
    graph = vc.graph_of(vessel)
    x = vc.mss_to_ours(np.array([1.2, 0.1, -0.03, 0.01, -0.02, 0.04, 0.0, 0.0, 2.0, 0.1, -0.2, 0.3]))
    u = np.array([0.3, 0.5, 0.0, 5.0, 0.0, -2.0, 0.0, 0.1, 0.2])
    symbol = ca.SX.sym("u", 9)
    jacobian = ca.Function("j", [symbol], [ca.jacobian(graph.step(x, symbol), symbol[0])])
    analytic = np.asarray(jacobian(u)).ravel()
    step = 1e-6
    up, down = u.copy(), u.copy()
    up[0] += step
    down[0] -= step
    numeric = (np.asarray(graph.step(x, up)).ravel() - np.asarray(graph.step(x, down)).ravel()) / (2 * step)
    assert np.abs(analytic).max() > 1e-3
    assert np.abs(analytic - numeric).max() <= 1e-6


def _remus_variant(**changes):
    tree = make_trees.torpedo("MunkCouplingsRemoved")
    for key, value in changes.items():
        tree["params"][key] = value
    return tree


@pytest.mark.parametrize("change, message", [
    ({"open_inputs": ["current.current_speed", "current.current_direction", "current.current_vertical_speed",
                      "site.nonexistent"]}, "no such open part input"),
    ({"diagnostic_outputs": ["mass_matrix", "not_a_quantity"]}, "diagnostic output 'not_a_quantity'"),
])
def test_refuses_an_open_input_or_diagnostic_output_that_names_nothing(change, message, tmp_path):
    with pytest.raises(ValueError, match=message):
        vc.build_variant(tmp_path, _remus_variant(**change))


def test_an_open_parameter_must_be_a_parameter_of_the_part(tmp_path):
    tree = _remus_variant()
    tree["children"]["current"]["params"]["open_parameters"] = ["current_speed", "not_a_parameter"]
    with pytest.raises(ValueError, match="not_a_parameter"):
        vc.build_variant(tmp_path, tree)


def test_a_part_input_left_open_without_the_vehicle_listing_it_is_refused(tmp_path):
    tree = _remus_variant(open_inputs=[])
    with pytest.raises(ValueError, match=r"slot current \(HorizontalCurrentYawRateTerms\): input 'current_speed'"):
        vc.build_variant(tmp_path, tree)


# ---------------------------------------------------------------------------
# Refusals name the slot and the part
# ---------------------------------------------------------------------------

def test_refuses_a_part_whose_input_no_slot_produces():
    with pytest.raises(ValueError) as error:
        vc.build("refused_missing_coupling")
    message = str(error.value)
    assert "slot rigid_body (HomogeneousSpheroid): input 'semi_major_axis' (1)" in message
    assert "no vehicle quantity and no earlier slot produces it" in message


def test_refuses_a_mistyped_open_input_and_quotes_it():
    with pytest.raises(ValueError) as error:
        vc.build("refused_open_input")
    message = str(error.value)
    assert "slot current (HorizontalCurrentYawRateTerms): input 'current_vertical_speed' (1)" in message
    assert "'current.current_vertical_sped'" in message


def test_refuses_a_hull_load_that_needs_a_quantity_the_restoring_part_does_not_give(tmp_path):
    tree = _remus_variant()
    tree["children"]["hydrodynamic_loads"].append(make_trees.node("SurgeResistanceIttc"))
    with pytest.raises(ValueError, match="wetted_surface"):
        vc.build_variant(tmp_path, tree)


def test_refuses_a_command_in_a_hull_load(tmp_path):
    """Only a force producer's unmatched input is a command: the same part in the hull-load list is refused."""
    tree = _remus_variant()
    tree["children"]["hydrodynamic_loads"].append(make_trees.node("PropellerLinearizedOpenWater"))
    with pytest.raises(ValueError):
        vc.build_variant(tmp_path, tree)


def test_refuses_an_empty_single_slot(tmp_path):
    tree = _remus_variant()
    del tree["children"]["restoring"]
    with pytest.raises(ValueError, match="slot 'restoring' is empty"):
        vc.build_variant(tmp_path, tree)


# --- parts written here: the payload shapes no real part gives ---------------------------------------------------

class _Context:
    """The two calls ``MarineCraft6DOF.initialize`` makes."""

    def __init__(self, children, parameters=None):
        self.children = children
        self.parameters = {"integration_max_step": 0.05, "open_inputs": [], "diagnostic_outputs": [],
                           **(parameters or {})}

    def get_component(self, slot):
        return self.children.get(slot, [])

    def get_parameter(self, name):
        return self.parameters[name]


def _real(plugin_path, class_name):
    """A real part initialised on its own defaults."""
    import importlib
    cls = getattr(importlib.import_module(plugin_path), class_name)
    defaults = {d.name: d.default_value for d in cls.PARAMETERS}
    part = cls()
    part.initialize(type("Context", (), {"get_parameter": lambda self, n: defaults[n],
                                         "get_component": lambda self, n: None})())
    return part


def _spheroid_children():
    return {
        "site": _real("more_dynamics.plugins.vehicles.hull_parts.site.site_at_latitude", "SiteAtLatitude"),
        "current": _real("more_dynamics.plugins.vehicles.hull_parts.current.no_current", "NoCurrent"),
        "hull_form": _real("more_dynamics.plugins.vehicles.hull_parts.hull_form.prolate_spheroid_main_dimensions",
                           "ProlateSpheroidMainDimensions"),
        "rigid_body": _real("more_dynamics.plugins.vehicles.hull_parts.rigid_body.homogeneous_spheroid", "HomogeneousSpheroid"),
        "added_mass": _real("more_dynamics.plugins.vehicles.hull_parts.added_mass.lamb_spheroid", "LambSpheroid"),
        "added_mass_coriolis": _real("more_dynamics.plugins.vehicles.hull_parts.added_mass_coriolis.kirchhoff_full",
                                     "KirchhoffFull"),
        "restoring": _real("more_dynamics.plugins.vehicles.hull_parts.restoring.submerged_neutral", "SubmergedNeutral"),
        "hydrodynamic_loads": [],
        "force_producers": [],
    }


class _Part:
    """A part with the payload ``builder`` makes."""

    def __init__(self, builder):
        self.builder = builder

    def graph(self):
        return self.builder()


def _hull_load(inputs, outputs, states=()):
    from rpp_plugin_types.more_dynamics import HydrodynamicsModel
    from more_dynamics.plugins.shared.payload_io import PayloadBuilder

    def build():
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        symbols = {name: io.input(name, size) for name, size in inputs}
        total = ca.vertcat(*[ca.sum1(ca.vec(s)) for s in symbols.values()]) if symbols else ca.SX(0.0)
        for name, size in outputs:
            io.output(name, ca.repmat(ca.sum1(total), size, 1) if size > 1 else ca.sum1(total))
        for name, size in states:
            state = io.state(name, size, [0.0] * size)
            io.state_dot(-state)
        return io.payload()
    return _Part(build)


def _assemble(**changes):
    cls, _ = _vehicle_class()
    children = {**_spheroid_children(), **changes}
    vehicle = cls()
    vehicle.initialize(_Context(children))
    return vehicle


def test_refuses_a_mis_sized_input():
    _, CompositionError = _vehicle_class()
    part = _hull_load([("added_mass_matrix", 9)], [("hydrodynamic_force", 6)])
    with pytest.raises(CompositionError, match=r"input 'added_mass_matrix' \(9\) but slot added_mass"):
        _assemble(hydrodynamic_loads=[part])


def test_refuses_a_quantity_two_single_slots_produce():
    _, CompositionError = _vehicle_class()
    from rpp_plugin_types.more_dynamics import CurrentModel
    from more_dynamics.plugins.shared.payload_io import PayloadBuilder

    def build():
        io = PayloadBuilder(CurrentModel.CasadyPayload())
        io.output("current_velocity", ca.SX.zeros(6))
        io.output("current_acceleration", ca.SX.zeros(6))
        io.output("water_density", ca.SX(1000.0))  # the site gives it already
        return io.payload()
    with pytest.raises(CompositionError, match="'water_density' is produced already by slot site"):
        _assemble(current=_Part(build))


def test_refuses_a_state_outside_the_force_producers():
    _, CompositionError = _vehicle_class()
    part = _hull_load([("velocity", 6)], [("hydrodynamic_force", 6)], states=[("memory", 1)])
    with pytest.raises(CompositionError, match="only force producers may declare states"):
        _assemble(hydrodynamic_loads=[part])


def test_refuses_a_first_output_that_is_not_a_six_value_force():
    _, CompositionError = _vehicle_class()
    part = _hull_load([("velocity", 6)], [("hydrodynamic_force", 3)])
    with pytest.raises(CompositionError, match="the first output must be the 6-value signed generalized force"):
        _assemble(hydrodynamic_loads=[part])


def test_a_part_that_declares_a_vehicle_quantity_as_its_own_parameter_is_only_the_site():
    """One quantity, one value: only a site part carries a water density, a gravity or a kinematic viscosity as
    its own parameter (Luka's hydrostatics and hydrodynamics parts, which carry their own, are the stated exception)."""
    vc.builder()
    import importlib
    import json
    registry = json.loads((Path(__file__).resolve().parents[2] / "plugins.json").read_text())
    for entry in registry["Plugins"]:
        module_path = entry["Path"][:-3].replace("/", ".")
        if ".hull_parts.site." in module_path or entry["Name"] in ("LinearSurfaceHydrostatics", "LinearSurfaceHydrodynamics",
                                                              "CrossflowSurfaceHydrodynamics", "JetNozzle", "HullVessel"):
            continue
        cls = getattr(importlib.import_module(module_path), entry["Name"])
        own = {d.name for d in getattr(cls, "PARAMETERS", [])}
        assert not own & {"water_density", "gravity", "kinematic_viscosity"}, (entry["Name"], own)


# ---------------------------------------------------------------------------
# step(): its own dynamics, RK4 with sub-steps, zero-order hold
# ---------------------------------------------------------------------------

def _odometry(x):
    from rpp_schema.rpp_common.Odometry3D import Odometry3D
    from more_transformations.more_casadi_transformations import EulerQuaternionTransforms
    odometry = Odometry3D()
    quaternion = np.asarray(EulerQuaternionTransforms.euler_to_quaternion(*x[3:6])).ravel()
    odometry.pose.position.x, odometry.pose.position.y, odometry.pose.position.z = x[0:3]
    (odometry.pose.orientation.w, odometry.pose.orientation.x, odometry.pose.orientation.y,
     odometry.pose.orientation.z) = quaternion
    odometry.twist.linear.x, odometry.twist.linear.y, odometry.twist.linear.z = x[6:9]
    odometry.twist.angular.x, odometry.twist.angular.y, odometry.twist.angular.z = x[9:12]
    return odometry


def test_step_equals_a_test_side_rk4_of_the_graphs_dynamics():
    from more_transformations.more_casadi_transformations import EulerQuaternionTransforms
    vessel = vc.build("otter")
    dynamics = vc.graph_of(vessel).step
    x0 = np.array([1.0, -2.0, 0.1, 0.05, -0.03, 0.7, 1.2, 0.1, 0.0, 0.02, -0.01, 0.05])
    after = vessel.step(_odometry(x0), [], 0.0, 0.5)
    x = x0.copy()
    for _ in range(10):  # 0.5 s in sub-steps of integration_max_step = 0.05 s
        k1 = np.asarray(dynamics(x, [])).ravel()
        k2 = np.asarray(dynamics(x + 0.025 * k1, [])).ravel()
        k3 = np.asarray(dynamics(x + 0.025 * k2, [])).ravel()
        k4 = np.asarray(dynamics(x + 0.05 * k3, [])).ravel()
        x = x + 0.05 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    q = after.pose.orientation
    euler = np.asarray([float(a) for a in EulerQuaternionTransforms.quaternion_to_euler([q.w, q.x, q.y, q.z])])
    got = np.r_[after.pose.position.x, after.pose.position.y, after.pose.position.z, euler,
                after.twist.linear.x, after.twist.linear.y, after.twist.linear.z,
                after.twist.angular.x, after.twist.angular.y, after.twist.angular.z]
    angle = np.angle(np.exp(1j * (got[3:6] - x[3:6])))
    others = [0, 1, 2, 6, 7, 8, 9, 10, 11]
    assert max(np.abs(got[others] - x[others]).max(), np.abs(angle).max()) <= 1e-10


def test_step_with_zero_time_returns_the_state():
    vessel = vc.build("otter")
    odometry = _odometry(np.array([1.0, -2.0, 0.1, 0.05, -0.03, 0.7, 1.2, 0.1, 0.0, 0.02, -0.01, 0.05]))
    same = vessel.step(odometry, [], 0.0, 0.0)
    assert same.twist.linear.x == pytest.approx(odometry.twist.linear.x, abs=1e-15)
    assert same.pose.position.y == pytest.approx(odometry.pose.position.y, abs=1e-15)


def test_step_refuses_a_vehicle_with_force_producer_states():
    vessel = vc.build("otter_jet_nozzle")
    odometry = _odometry(np.zeros(12))
    with pytest.raises(NotImplementedError, match="force-producer states"):
        vessel.step(odometry, [], 0.0, 0.1)


# ---------------------------------------------------------------------------
# Luka's own parts and tree in the same builder
# ---------------------------------------------------------------------------

def test_lukas_hydrostatics_and_hydrodynamics_plug_in_by_his_names():
    vc.builder()
    from more_dynamics.models.hydrostatics import (linear_surface_hydrostatics_casadi,
                                                   preprocess_linear_surface_hydrostatics)
    from more_dynamics.plugins.hydrostatics.linear_surface_hydrostatics import LinearSurfaceHydrostatics
    vessel = vc.build("otter_lukas_parts")
    assert vessel.wiring.get("restoring.pose") == "vehicle"
    assert vessel.wiring.get("hydrodynamic_loads.0.velocity") == "vehicle"
    pose = np.array([0.0, 0.0, 0.03, 0.05, -0.02, 0.4])
    out = vc.output_entries(vessel, np.r_[pose, np.zeros(6)], np.zeros(0))
    defaults = {d.name: d.default_value for d in LinearSurfaceHydrostatics.PARAMETERS}
    model = linear_surface_hydrostatics_casadi(preprocess_linear_surface_hydrostatics(**defaults))
    expected = np.asarray(model(pose)).ravel()[:6]
    assert np.abs(out["restoring_force"] - expected).max() <= 1e-12


def test_lukas_jet_nozzle_is_a_force_producer_with_two_states():
    vc.builder()
    from more_common.casadi_graph import RppCasadiGraph
    from more_dynamics.plugins.force_producers.jet_nozzle import JetNozzle
    vessel = vc.build("otter_jet_nozzle")
    assert [d.name for d in vessel.graph().stateDescription] == ["thrust", "nozzle_angle", "state"]
    nozzle = JetNozzle()
    defaults = {d.name: d.default_value for d in JetNozzle.PARAMETERS}
    nozzle.initialize(type("Context", (), {"get_parameter": lambda self, n: defaults[n]})())
    reference = RppCasadiGraph(nozzle.graph())
    state, command = np.array([0.6, 0.2]), np.array([0.8, -0.3])
    x = np.r_[0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 1.0, 0.1, 0.0, 0.0, 0.0, 0.02]
    graph = vc.graph_of(vessel)
    out = vc.output_entries(vessel, np.r_[state, x], command)
    assert np.abs(out["force_producers.0.generated_thrust"]
                  - np.asarray(reference.output(state, command)).ravel()).max() <= 1e-15
    assert np.abs(np.asarray(graph.step(np.r_[state, x], command)).ravel()[:2]
                  - np.asarray(reference.step(state, command)).ravel()).max() <= 1e-15


def test_lukas_hull_vessel_tree_still_builds_in_the_same_builder():
    script = Path(__file__).resolve().parents[2] / ".rppws" / "script_descriptions" / "hull_vessel_simulation.json"
    vessel = vc.build(None, script)
    assert type(vessel).__name__ == "HullVessel"
    assert vc.graph_of(vessel).num_states == 14
