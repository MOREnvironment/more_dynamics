"""Assembly gates of the three vehicle types (``SpheroidAuv``, ``Monohull``,
``Catamaran``): the slots of ``HullVessel``, children connected by the names
in their descriptions, the command list, refusals that name the slot and the
child, ``step``, the buoyancy methods of the submerged restoring, the
monohull, and Luka's own plugins in the same vehicles.

The refusals that no real plugin can provoke (a mis-sized input, a quantity
produced twice, a state in a hydrodynamics model) use children written in
the test with the same payload shape.

Author:    Enio Krizman
Date:      2026-10-09
"""

import casadi as ca
import numpy as np
import pytest

import vehicle_contract as vc

REMUS_STATE = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 1.2, 0.05, 0.02, 0.01, -0.02, 0.03])  # [eta; nu], a surging turn
CURRENT = list(vc.CURRENT_NAMES)


# ---------------------------------------------------------------------------
# Slots, couplings by name, the command list
# ---------------------------------------------------------------------------

def test_the_three_types_take_exactly_the_slots_of_hull_vessel_without_sensors():
    p = vc.plugins()
    from more_dynamics.plugins.vehicle_models.hull_vessel import HullVessel
    for cls in (p.SpheroidAuv, p.Monohull, p.Catamaran):
        assert cls.COMPONENTS == {key: value for key, value in HullVessel.COMPONENTS.items() if key != "sensors"}, cls


@pytest.mark.parametrize("make, names", [
    (lambda: vc.remus(), CURRENT + ["desired_wrench"]),
    (lambda: vc.remus(full=True), CURRENT + ["fin_deflection_command", "shaft_speed_command"]),
    (lambda: vc.remus(full=True, current="none"), ["fin_deflection_command", "shaft_speed_command"]),
    (lambda: vc.otter(), []),
    (lambda: vc.otter(current="full_rotation_rate"), CURRENT),
], ids=["remus_hull", "remus_full", "remus_full_no_current", "otter", "otter_current"])
def test_input_list_of_each_vehicle(make, names):
    """The current values come first when the current form asks for them, then every actuator input that is no
    vehicle quantity, in list order."""
    assert vc.input_names(make()) == names


def test_every_input_of_a_child_is_connected_by_name():
    """The hydrodynamics read the matrices of the vehicle, the weight and the centres of the hydrostatics, and the
    site's density, each by the name in their description."""
    vessel = vc.remus(coriolis="munk_couplings_removed")
    known = vessel.signals(REMUS_STATE, np.zeros(9))
    for name in ("pose", "velocity", "relative_velocity", "mass_matrix", "rigid_body_mass_matrix",
                 "added_mass_matrix", "water_density", "gravity", "weight", "center_of_buoyancy", "span",
                 "section_beam", "length", "draft", "kinematic_viscosity"):
        assert name in known, name
    assert known["water_density"][0] == pytest.approx(vc.REMUS["water_density"], rel=1e-15)



def test_the_output_is_the_twelve_states_of_hull_vessel():
    vessel = vc.otter()
    payload = vessel.graph()
    assert [(d.name, d.size) for d in payload.outputDescription] == [("output", 12)]
    assert [(d.name, d.size) for d in payload.stateDescription] == [("state", 12)]


def test_the_mass_matrix_equals_the_frozen_reference():
    import csv
    vessel = vc.remus(coriolis="munk_couplings_removed")
    mass = vc.matrix(vessel.signals(np.zeros(12), np.zeros(9))["mass_matrix"])
    with open(vc.DATA / "remus100" / "remus100_mass_matrix_consistent.csv", newline="") as handle:
        table = np.asarray(list(csv.reader(handle))[1:], float)
    assert np.abs(mass - table[:, :6]).max() <= 1e-9
    assert np.array_equal(vessel.signals(np.zeros(12), np.zeros(9))["relative_velocity"], np.zeros(6))


def test_gradient_through_a_current_input_equals_a_central_difference():
    vessel = vc.remus(coriolis="munk_couplings_removed")
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


# ---------------------------------------------------------------------------
# Options: forms, refusals
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("parameter, value", [("current_form", "full_attitude"), ("coriolis_form", "none"),
                                              ("added_mass_form", "database")])
def test_refuses_an_unknown_form_and_names_the_choices(parameter, value):
    p = vc.plugins()
    with pytest.raises(p.CompositionError, match=parameter):
        vc.remus(vehicle_params={parameter: value})


def test_the_coriolis_form_changes_only_the_added_mass_coriolis_term():
    """``munk_couplings_removed`` is the MSS shortcut, ``kirchhoff_full`` the physics form: only the acceleration
    differs, the kinematics are the same."""
    a = vc.remus(coriolis="munk_couplings_removed")
    b = vc.remus(coriolis="kirchhoff_full")
    x, u = REMUS_STATE, np.zeros(9)
    fa, fb = (np.asarray(vc.graph_of(v).step(x, u)).ravel() for v in (a, b))
    assert np.array_equal(fa[:6], fb[:6])  # eta_dot = J(eta) nu does not read the Coriolis term
    assert np.abs(fa[6:] - fb[6:]).max() > 1e-6


def test_refuses_an_empty_slot_and_names_it():
    p = vc.plugins()
    from rpp_py.context import ComponentContext
    from rpp_py.parameter_handler import ParameterHandler
    params = ParameterHandler.resolve_params(p.SpheroidAuv.PARAMETERS, {})
    context = ComponentContext(instance=p.SpheroidAuv(), params=params, spec=p.SpheroidAuv.COMPONENTS,
                               subcomponents={"hydrodynamics": [vc.child(p.AuvHullLoads)], "actuators": []})
    with pytest.raises(p.CompositionError, match="slot 'hydrostatics' is empty"):
        context.initialize()


def test_refuses_a_hydrodynamics_child_whose_input_no_quantity_gives():
    """The ITTC surge resistance needs a wetted surface; none of the three types gives one, and the message names
    the slot, the plugin and the quantity."""
    p = vc.plugins()
    with pytest.raises(p.CompositionError, match=r"slot hydrodynamics \(SurfaceHullLoads\): input 'wetted_surface'"):
        vc.otter(hydrodynamics=(p.SurfaceHullLoads, {**vc.pick(p.SurfaceHullLoads, vc.OTTER),
                                                      "surge_resistance": "ittc"}))


def test_refuses_a_command_in_the_hydrodynamics_slot():
    """Only an actuator's unmatched input is a command: the propeller in the hydrodynamics slot is refused."""
    p = vc.plugins()
    with pytest.raises(p.CompositionError, match="no vehicle quantity and no earlier slot produces it"):
        vc.build_context(p.SpheroidAuv, hydrostatics=p.SubmergedRestoring, hydrodynamics=p.Propeller)


# --- children written here: the payload shapes no real plugin gives --------------------------------------------------

class _Child:
    def __init__(self, builder):
        self.builder = builder

    def graph(self):
        return self.builder()


class _StandIn:
    """What a vehicle type asks of its context: the children by slot and the declared defaults."""

    def __init__(self, cls, children):
        self.defaults = {d.name: d.default_value for d in cls.PARAMETERS}
        self.children = children

    def get_component(self, slot):
        return self.children.get(slot, [])

    def get_parameter(self, name, default=None):
        return self.defaults.get(name, default)


def _hydrodynamics(inputs, outputs, states=()):
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
    return _Child(build)


def _assemble(hydrodynamics):
    """The Catamaran over the real surface restoring and the given hydrodynamics child."""
    p = vc.plugins()
    restoring = vc.child(p.SurfaceRestoring)
    restoring.initialize()
    vehicle = p.Catamaran()
    vehicle.initialize(_StandIn(p.Catamaran, {"hydrostatics": restoring.get_instance(),
                                              "hydrodynamics": hydrodynamics, "actuators": []}))
    return vehicle


def test_refuses_a_mis_sized_input():
    p = vc.plugins()
    child = _hydrodynamics([("rigid_body_mass_matrix", 9)], [("hydrodynamic_force", 6)])
    with pytest.raises(p.CompositionError, match=r"input 'rigid_body_mass_matrix' \(9\) but Catamaran gives 36"):
        _assemble(child)


def test_refuses_a_quantity_two_slots_produce():
    p = vc.plugins()
    child = _hydrodynamics([("velocity", 6)], [("hydrodynamic_force", 6), ("water_density", 1)])
    with pytest.raises(p.CompositionError, match="'water_density' is produced already by Catamaran"):
        _assemble(child)


def test_refuses_a_state_in_a_hydrodynamics_model():
    p = vc.plugins()
    child = _hydrodynamics([("velocity", 6)], [("hydrodynamic_force", 6)], states=[("memory", 1)])
    with pytest.raises(p.CompositionError, match="declares no states"):
        _assemble(child)


def test_refuses_a_first_output_that_is_not_a_six_value_force():
    p = vc.plugins()
    child = _hydrodynamics([("velocity", 6)], [("hydrodynamic_force", 3)])
    with pytest.raises(p.CompositionError, match="the first output must be the 6-value signed generalized force"):
        _assemble(child)


def test_only_the_site_quantities_are_the_vehicles_own_parameters():
    """One quantity, one value: no plugin but the vehicle types declares a water density, a gravity or a
    kinematic viscosity (Luka's hydrostatics and hydrodynamics plugins, which carry their own, are the stated
    exception)."""
    import importlib
    import json
    from pathlib import Path
    vc.plugins()
    registry = json.loads((Path(__file__).resolve().parents[2] / "plugins.json").read_text())
    lukas = ("LinearSurfaceHydrostatics", "LinearSurfaceHydrodynamics", "CrossflowSurfaceHydrodynamics", "JetNozzle",
             "HullVessel")
    for entry in registry["Plugins"]:
        if entry["Name"] in lukas or "/vehicle_models/" in entry["Path"]:
            continue
        cls = getattr(importlib.import_module(entry["Path"][:-3].replace("/", ".")), entry["Name"])
        own = {d.name for d in getattr(cls, "PARAMETERS", [])}
        assert not own & {"water_density", "gravity", "kinematic_viscosity", "latitude"}, (entry["Name"], own)


# ---------------------------------------------------------------------------
# Buoyancy methods of the submerged restoring
# ---------------------------------------------------------------------------

def _derivative(vessel, x):
    return np.asarray(vc.graph_of(vessel).step(x, np.zeros(vc.input_size(vessel)))).ravel()


def _restoring_vehicle(**restoring):
    p = vc.plugins()
    v = vc.remus_values()
    return vc.build_context(
        p.SpheroidAuv, {**vc.pick(p.SpheroidAuv, v), "current_form": "none"},
        hydrostatics=(p.SubmergedRestoring, restoring), hydrodynamics=(p.AuvHullLoads, vc.pick(p.AuvHullLoads, v)),
        actuators=[])


def test_neutral_is_the_default_and_keeps_a_body_at_rest_at_rest():
    vessel = _restoring_vehicle(center_of_buoyancy=vc.REMUS["center_of_buoyancy"])
    assert np.abs(_derivative(vessel, np.zeros(12))).max() <= 1e-12


def test_buoyancy_from_the_volume_of_the_mass_equals_neutral():
    """``B = rho g nabla`` with ``nabla = m / rho`` is ``B = W``: the same vehicle as ``neutral``."""
    p = vc.plugins()
    mass = 31.9
    neutral = _restoring_vehicle(center_of_buoyancy=vc.REMUS["center_of_buoyancy"])
    volume = _restoring_vehicle(buoyancy_method="from_volume", displaced_volume=mass / vc.REMUS["water_density"],
                                center_of_buoyancy=vc.REMUS["center_of_buoyancy"])
    x = REMUS_STATE
    assert np.abs(_derivative(neutral, x) - _derivative(volume, x)).max() <= 1e-12


def test_a_given_buoyancy_equal_to_the_weight_equals_neutral():
    mass = 31.9
    g = float(neutral_signals()["gravity"][0])
    neutral = _restoring_vehicle(center_of_buoyancy=vc.REMUS["center_of_buoyancy"])
    given = _restoring_vehicle(buoyancy_method="given", buoyancy=mass * g,
                               center_of_buoyancy=vc.REMUS["center_of_buoyancy"])
    x = REMUS_STATE
    assert np.abs(_derivative(neutral, x) - _derivative(given, x)).max() <= 1e-12


def neutral_signals():
    return _restoring_vehicle(center_of_buoyancy=vc.REMUS["center_of_buoyancy"]).signals(np.zeros(12), [])


def test_a_positively_buoyant_vehicle_is_pushed_up_and_trims_nose_up():
    """``B`` 5 % above ``W`` with the centre of buoyancy 0.05 m forward of the centre of gravity (which is at
    x = 0): at level attitude the net force is up (negative z in NED), ``-(B - W)``, and the moment about y is
    ``r_bb,x B - r_bg,x W`` lifting the nose (Fossen 2011, eq. 4.5, p. 60)."""
    mass = 31.9
    g = float(neutral_signals()["gravity"][0])
    weight = mass * g
    assert vc.REMUS["body_center_of_gravity"][0] == 0.0 and vc.REMUS["body_center_of_gravity"][1] == 0.0
    vessel = _restoring_vehicle(buoyancy_method="given", buoyancy=1.05 * weight, center_of_buoyancy=[0.05, 0.0, 0.0])
    force = vessel.signals(np.zeros(12), [])["hydrostatics.restoring_force"]
    assert force[2] == pytest.approx(-0.05 * weight, rel=1e-12)
    assert force[4] == pytest.approx(0.05 * 1.05 * weight, rel=1e-12)
    assert abs(force[3]) <= 1e-12 and abs(force[0]) <= 1e-12 and abs(force[1]) <= 1e-12
    assert _derivative(vessel, np.zeros(12))[6 + 2] < 0.0  # heave acceleration upwards


@pytest.mark.parametrize("method", ["from_volume", "given", "neutral"])
def test_each_method_reads_only_its_own_parameter(method):
    p = vc.plugins()
    # all three of its own numbers are declared with defaults; the vehicle builds with only the method chosen
    vessel = _restoring_vehicle(buoyancy_method=method)
    assert vessel.signals(np.zeros(12), [])["buoyancy"].shape == (1,)


def test_refuses_an_unknown_buoyancy_method():
    with pytest.raises(ValueError, match="buoyancy_method"):
        _restoring_vehicle(buoyancy_method="floating")


# ---------------------------------------------------------------------------
# The monohull
# ---------------------------------------------------------------------------

def _monohull(**vehicle_params):
    p = vc.plugins()
    return vc.build_context(p.Monohull, vehicle_params,
                            hydrostatics=(p.SurfaceRestoring, {"hull_count": 1}), hydrodynamics=p.SurfaceHullLoads)


def test_monohull_builds_on_its_defaults_and_takes_no_input():
    vessel = _monohull()
    assert vc.input_names(vessel) == []
    assert set(vc.plugins().Monohull.COMPONENTS) == {"hydrostatics", "hydrodynamics", "actuators"}


def test_monohull_at_rest_in_still_water_at_the_equilibrium_pose_stays_at_rest():
    vessel = _monohull()
    assert np.abs(_derivative(vessel, np.zeros(12))).max() <= 1e-12


def test_monohull_mass_matrix_is_symmetric_positive_definite_and_keeps_the_given_mass():
    vessel = _monohull()
    signals = vessel.signals(np.zeros(12), [])
    mass = vc.matrix(signals["mass_matrix"])
    assert np.allclose(mass, mass.T) and np.all(np.linalg.eigvalsh(mass) > 0)
    assert signals["mass"][0] == pytest.approx(806.0, rel=1e-15)


def test_monohull_restoring_equals_the_one_hull_surface_block_at_its_equilibrium_draft():
    """The draft is ``nabla / (n C_b L B)`` with ``n = 1``: the same plugin with two hulls on the same hull gives half
    the draft per hull."""
    one = _monohull().signals(np.zeros(12), [])
    p = vc.plugins()
    two = vc.build_context(p.Monohull, {}, hydrostatics=(p.SurfaceRestoring, {"hull_count": 1}),
                           hydrodynamics=p.SurfaceHullLoads).signals(np.zeros(12), [])
    assert one["draft"][0] == pytest.approx(two["draft"][0], rel=1e-15)
    volume = 806.0 / 1025.0
    assert one["draft"][0] == pytest.approx(volume / (0.233 * 2.15 * 5.2), rel=1e-12)


def test_monohull_is_a_different_vehicle_from_the_catamaran_with_the_same_hull_loads():
    mono = _monohull().signals(np.zeros(12), [])
    cat = vc.otter().signals(np.zeros(12), [])
    assert abs(mono["mass"][0] - cat["mass"][0]) > 1.0


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
    vessel = vc.otter()
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
    vessel = vc.otter()
    odometry = _odometry(np.array([1.0, -2.0, 0.1, 0.05, -0.03, 0.7, 1.2, 0.1, 0.0, 0.02, -0.01, 0.05]))
    same = vessel.step(odometry, [], 0.0, 0.0)
    assert same.twist.linear.x == pytest.approx(odometry.twist.linear.x, abs=1e-15)
    assert same.pose.position.y == pytest.approx(odometry.pose.position.y, abs=1e-15)


def test_step_refuses_a_vehicle_with_actuator_states():
    from more_dynamics.plugins.force_producers.jet_nozzle import JetNozzle
    vessel = vc.otter(actuators=[JetNozzle])
    with pytest.raises(NotImplementedError, match="actuator states"):
        vessel.step(_odometry(np.zeros(12)), [], 0.0, 0.1)


# ---------------------------------------------------------------------------
# Luka's own plugins in the same vehicles
# ---------------------------------------------------------------------------

def test_lukas_hydrostatics_and_hydrodynamics_plug_in_by_his_names():
    from more_dynamics.models.hydrostatics import (linear_surface_hydrostatics_casadi,
                                                   preprocess_linear_surface_hydrostatics)
    from more_dynamics.plugins.hydrodynamics.linear_surface_hydrodynamics import LinearSurfaceHydrodynamics
    from more_dynamics.plugins.hydrostatics.linear_surface_hydrostatics import LinearSurfaceHydrostatics
    vessel = vc.otter(hydrostatics=LinearSurfaceHydrostatics, hydrodynamics=LinearSurfaceHydrodynamics)
    pose = np.array([0.0, 0.0, 0.03, 0.05, -0.02, 0.4])
    out = vessel.signals(np.r_[pose, np.zeros(6)], [])
    defaults = {d.name: d.default_value for d in LinearSurfaceHydrostatics.PARAMETERS}
    model = linear_surface_hydrostatics_casadi(preprocess_linear_surface_hydrostatics(**defaults))
    expected = np.asarray(model(pose)).ravel()[:6]
    assert np.abs(out["hydrostatics.restoring_force"] - expected).max() <= 1e-12


def test_lukas_jet_nozzle_is_an_actuator_with_two_states():
    from more_common.casadi_graph import RppCasadiGraph
    from more_dynamics.plugins.force_producers.jet_nozzle import JetNozzle
    vessel = vc.otter(actuators=[JetNozzle])
    assert [d.name for d in vessel.graph().stateDescription] == ["thrust", "nozzle_angle", "state"]
    nozzle = vc.context_of(JetNozzle)
    nozzle.initialize()
    reference = RppCasadiGraph(nozzle.get_instance().graph())
    state, command = np.array([0.6, 0.2]), np.array([0.8, -0.3])
    x = np.r_[0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 1.0, 0.1, 0.0, 0.0, 0.0, 0.02]
    graph = vc.graph_of(vessel)
    out = vessel.signals(np.r_[state, x], command)
    assert np.abs(out["actuators.0.generated_thrust"] - np.asarray(reference.output(state, command)).ravel()).max() <= 1e-15
    assert np.abs(np.asarray(graph.step(np.r_[state, x], command)).ravel()[:2]
                  - np.asarray(reference.step(state, command)).ravel()).max() <= 1e-15
