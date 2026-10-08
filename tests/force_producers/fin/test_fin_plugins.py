"""The fin parts, the servo, the lifting fin and the force-producer set as
plugins of Luka's form, built from ``.rppws`` trees by rpp's own builder.

* Every plugin of the fin family is listed in ``plugins.json``, one class per
  file, every parameter with a default; the defaults are the printed numbers
  they cite; the two parts waiting on a source refuse to initialise saying so.
* Linked parts (``IsLinked``): a fin kept as one top-level part and linked
  into a vehicle's force-producer list gives the same vehicle as the fin
  written in place; its servo state precedes the vehicle's twelve states, its
  command is the vehicle's last input, its wrench inside the vehicle equals
  the model-layer ``lifting_fin_casadi`` assembled from the same parts and
  numbers, and its servo derivative equals the servo model's.
* The servo plugin's four settings (the lag, the rate limit and the angle
  limit on or off) equal the servo model with the same selectors; settings the
  model refuses are refused when the vehicle is built.
* The REMUS 100 rudder and stern plane as two ``LiftingFin`` compositions in a
  ``ForceProducerSet`` equal the landed ``FinPairsDeflectionOnly`` part on the
  frozen MSS rows (and the MATLAB wrench of those rows, scaled to the vehicle's
  one water density).
* A fin refuses a part that does not fit its slot, a child parameter left
  open, and a part waiting on a source; a set refuses a command map of the
  wrong shape; the plugin layer holds no trace of the retired part base.

Reference numbers: the frozen MSS rows and the printed Prestero values of
``tests/data/force_producers/fin_parts`` and the REMUS parameter file of
``tests/data/vehicles``.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc
    thesis. Eqs. 4.37, 4.41-4.43, pp. 31-33; Table A.5, p. 103.
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and
    Motion Control. Wiley. Eq. 12.226, p. 400.
[Sarhadi 2026] Sarhadi, P. (2026). Simple yet effective anti-windup techniques
    for amplitude and rate saturation: an AUV case study. arXiv:2601.01302v2,
    Fig. 4, p. 4.
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 228-254 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from fin_parts_contract import (E_X, E_Y, G1_TOLERANCE, G2_TOLERANCE, SERVO_SETTINGS, call, fin_from_forms, geometry,
                                max_diff, mss_columns, mss_constant, mss_fins, part, part_declared, prestero_value,
                                producer_set, rng, tau)

LIBRARY = Path(__file__).resolve().parents[2].parent
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vehicles"))
import vehicle_contract as vc  # noqa: E402
from vehicle_contract import make_trees  # noqa: E402

FIN_PLUGINS = {
    "Servo": "force_producers/fin/servo/servo", "FinInflowRigidPoint": "force_producers/fin/inflow/rigid_point",
    "FinInflowTranslational": "force_producers/fin/inflow/translational",
    "FinFlowAngleSmallAngle": "force_producers/fin/flow_angle/small_angle", "FinFlowAngleNone": "force_producers/fin/flow_angle/none",
    "FinInterferenceNone": "force_producers/fin/interference/none",
    "FinInterferenceSlenderBody": "force_producers/fin/interference/slender_body",
    "FinSectionQuadraticDrag": "force_producers/fin/section/quadratic_drag",
    "FinSectionLinearSection": "force_producers/fin/section/linear_section",
    "FinSectionLiftingLine": "force_producers/fin/section/lifting_line",
    "LiftingFin": "force_producers/fin/lifting_fin", "ForceProducerSet": "force_producers/force_producer_set"}
WAITING = {"FinInterferenceSlenderBody": "Pitts 1957", "FinSectionLiftingLine": "lifting-line induced drag"}
DEG = np.pi / 180.0
# the servo numbers of the plugin defaults (Sarhadi 2026, Fig. 4, p. 4)
SERVO_NUMBERS = {"time_constant": 0.1, "max_deflection": 20.0 * DEG, "max_rate": 30.0 * DEG}
RHO = make_trees.REMUS["water_density"]


def _classes():
    vc.builder()
    registry = json.loads((LIBRARY / "plugins.json").read_text())
    entries = {e["Name"]: e["Path"] for e in registry["Plugins"]}
    found = {}
    for name, path in FIN_PLUGINS.items():
        assert entries.get(name) == f"more_dynamics/plugins/{path}.py", (name, entries.get(name))
        module = importlib.import_module(entries[name][:-3].replace("/", "."))
        found[name] = getattr(module, name)
    return found


# --------------------------------------------------------------------------
# The plugin family
# --------------------------------------------------------------------------
def test_every_fin_plugin_is_listed_and_has_a_class_and_one_class_per_file():
    classes = _classes()
    assert set(classes) == set(FIN_PLUGINS)
    for name, path in FIN_PLUGINS.items():
        text = (LIBRARY / "more_dynamics" / "plugins" / f"{path}.py").read_text()
        assert text.count("\nclass ") == 1, name


def test_every_parameter_has_a_default():
    for name, cls in _classes().items():
        for description in cls.PARAMETERS:
            assert description.default_value is not None, (name, description.name)


def _default(cls, name):
    return next(d.default_value for d in cls.PARAMETERS if d.name == name)


def test_defaults_are_the_printed_numbers_they_cite():
    classes = _classes()
    fin, servo = classes["LiftingFin"], classes["Servo"]
    assert _default(fin, "fin_position") == [prestero_value("fin_position_x"), 0.0, 0.0]  # Table A.5, p. 103
    assert _default(fin, "fin_area") == prestero_value("fin_area")
    for section in ("FinSectionQuadraticDrag", "FinSectionLinearSection"):
        assert _default(classes[section], "lift_slope") == prestero_value("lift_slope")
    for name, value in SERVO_NUMBERS.items():  # Sarhadi 2026, Fig. 4, p. 4
        assert abs(_default(servo, name) - value) <= 1e-15, name
    assert [_default(servo, n) for n in ("dynamics", "rate_limit", "angle_limit")] == \
        [SERVO_SETTINGS["lag_rate_angle"][n] for n in ("dynamics", "rate_limit", "angle_limit")]


@pytest.mark.parametrize("name, source", sorted(WAITING.items()))
def test_a_part_waiting_on_a_source_refuses_to_initialise_saying_so(name, source):
    cls = _classes()[name]
    assert [d.name for d in cls.PARAMETERS] == []
    with pytest.raises(NotImplementedError, match=f"waiting on source: {source}"):
        cls().initialize(None)


def test_the_plugin_layer_holds_no_trace_of_the_retired_part_base():
    offenders = []
    for path in (LIBRARY / "more_dynamics" / "plugins").rglob("*.py"):
        text = path.read_text()
        if any(word in text for word in ("BlockPayload", ".block()", "def block(", "fin_part_plugin",
                                         "proportional_force_fin", "ProportionalForceFin")):
            offenders.append(str(path.relative_to(LIBRARY)))
    assert not offenders, offenders
    retired = [p for p in ("fin_part_plugin", "proportional_force_fin", "fin_inflow_rigid_point",
                           "fin_section_quadratic_drag") if (LIBRARY / "more_dynamics" / "plugins" /
                                                            "force_producers" / f"{p}.py").exists()]
    assert not retired, retired


# --------------------------------------------------------------------------
# Vehicle helpers
# --------------------------------------------------------------------------
def _evaluate(vessel, x, u):
    """The vehicle's state derivative and its named output entries at stacked ``x`` and ``u``."""
    graph = vc.graph_of(vessel)
    return np.asarray(graph.step(x, u)).ravel(), vc.output_entries(vessel, x, u)


def _states(vessel, vehicle_state):
    """The vehicle's stacked state: the children's states (zeros) then the twelve vehicle states."""
    n = sum(d.size for d in vessel.graph().stateDescription) - 12
    return np.concatenate([np.zeros(n), vehicle_state])


def _random_states(n, child=0):
    g = rng()
    for _ in range(n):
        yield np.concatenate([g.uniform(-0.3, 0.3, child), g.uniform(-0.2, 0.2, 6), g.uniform(-2.5, 2.5, 3),
                              g.uniform(-0.5, 0.5, 3)]), g.uniform(-0.6, 0.6)


def _fin_input(command):
    """Inputs of a hull-plus-one-fin vehicle: three current inputs (still water), the wrench, the fin command."""
    return np.concatenate([np.zeros(3), np.zeros(6), [command]])


def _reference_fin(servo="lag_rate_angle"):
    """The model-layer fin of ``make_trees.prestero_fin``: the plugin defaults' servo, rigid-point inflow, small-angle
    flow angle, no interference, quadratic-drag section."""
    declared = {"max_deflection", "max_rate", "time_constant"}
    return fin_from_forms(
        {"servo": servo, "inflow": "rigid_point", "flow_angle": "small_angle", "interference": "none",
         "section": "quadratic_drag"},
        {"servo": {k: v for k, v in SERVO_NUMBERS.items() if k in declared},
         "section": {"lift_slope": prestero_value("lift_slope")}},
        geometry([prestero_value("fin_position_x"), 0.0, 0.0], E_Y, prestero_value("fin_area"), E_X))


# --------------------------------------------------------------------------
# Linked parts (IsLinked)
# --------------------------------------------------------------------------
def test_a_linked_fin_is_a_force_producer_with_its_servo_state_before_the_vehicle_states():
    vessel = vc.build("remus100_hull_fin_linked")
    graph = vessel.graph()
    assert [(d.name, d.size) for d in graph.stateDescription] == [("servo_state", 1), ("state", 12)]
    assert vc.input_names(vessel)[-2:] == ["desired_wrench", "fin_angle_command"]
    assert vc.input_names(vessel)[-1] == "fin_angle_command"
    assert [d.name for d in graph.outputDescription][:2] == ["output", "force_producers.1.generated_force"]


def test_a_linked_fin_gives_the_same_vehicle_as_the_fin_written_in_place():
    linked, inline = vc.build("remus100_hull_fin_linked"), vc.build("remus100_hull_fin_inline")
    worst = 0.0
    for x, command in _random_states(25, child=1):
        u = _fin_input(command)
        a, b = _evaluate(linked, x, u), _evaluate(inline, x, u)
        worst = max(worst, max_diff(a[0], b[0]), max_diff(a[1]["force_producers.1.generated_force"],
                                                          b[1]["force_producers.1.generated_force"]))
    assert worst == 0.0


def test_the_fin_wrench_in_the_vehicle_is_the_model_layer_fin_on_the_same_parts():
    vessel = vc.build("remus100_hull_fin_linked")
    fin = _reference_fin()
    worst = 0.0
    for x, command in _random_states(20, child=1):
        _, out = _evaluate(vessel, x, _fin_input(command))
        reference = tau(fin, command, x[1 + 6:1 + 12], RHO, state=x[:1])  # nu_r = nu in still water
        worst = max(worst, max_diff(out["force_producers.1.generated_force"], reference))
    assert worst <= G2_TOLERANCE, worst


def test_the_servo_derivative_in_the_vehicle_is_the_servo_models():
    vessel = vc.build("remus100_hull_fin_linked")
    servo = part("servo", "lag_rate_angle", SERVO_NUMBERS)
    worst = 0.0
    for x, command in _random_states(20, child=1):
        xdot, out = _evaluate(vessel, x, _fin_input(command))
        reference = call(servo, command=command, servo_state=x[:1])
        worst = max(worst, abs(xdot[0] - reference["servo_state_dot"][0]),
                    abs(out["force_producers.1.deflection"][0] - reference["deflection"][0]))
    assert worst <= 1e-14, worst


# --------------------------------------------------------------------------
# The servo's switches
# --------------------------------------------------------------------------
def _vehicle_with_servo(tmp_path, **settings):
    fin = make_trees.prestero_fin()
    fin["children"]["servo"] = make_trees.node("Servo", settings)
    tree = make_trees.torpedo("MunkCouplingsRemoved", diagnostics=make_trees.FIN_DIAGNOSTICS,
                              producers=[make_trees.node("PrescribedWrench"), fin])
    return vc.build_variant(tmp_path, tree)


@pytest.mark.parametrize("name", sorted(SERVO_SETTINGS))
def test_each_servo_setting_equals_the_servo_model_with_the_same_selectors(name, tmp_path):
    settings = SERVO_SETTINGS[name]
    vessel = _vehicle_with_servo(tmp_path, **settings)
    lag = settings["dynamics"] == "first_order_lag"
    declared = {d.name for d in part_declared("servo", name)}
    servo = part("servo", name, {k: v for k, v in SERVO_NUMBERS.items() if k in declared})
    assert [d.name for d in vessel.graph().stateDescription] == (["servo_state", "state"] if lag else ["state"])
    worst = 0.0
    for x, command in _random_states(20, child=1 if lag else 0):
        command *= 2.0  # beyond the 20 degree limit now and then
        xdot, out = _evaluate(vessel, x, _fin_input(command))
        reference = call(servo, command=command, servo_state=x[:1] if lag else [])
        worst = max(worst, abs(out["force_producers.1.deflection"][0] - reference["deflection"][0]))
        if lag:
            worst = max(worst, abs(xdot[0] - reference["servo_state_dot"][0]))
    assert worst <= 1e-14, worst


@pytest.mark.parametrize("settings, word", [
    ({"dynamics": "none", "rate_limit": True}, "rate_limit"),
    ({"dynamics": "none", "rate_limit": False, "angle_limit": "on_output"}, "angle_limit"),
    ({"angle_limit": "on_both"}, "angle_limit"),
])
def test_a_servo_setting_the_model_refuses_is_refused_when_the_vehicle_is_built(settings, word, tmp_path):
    with pytest.raises(Exception, match=word):
        _vehicle_with_servo(tmp_path, **settings)


# --------------------------------------------------------------------------
# A fin set against the landed fin pairs, on the frozen MSS rows
# --------------------------------------------------------------------------
def _rows():
    ref = mss_fins()
    return ref, np.column_stack([ref["ui1"], ref["ui2"]]), mss_columns(ref, "nu_r", 6), mss_columns(ref, "tau", 6)


def _forces(vessel, commands, nu_r):
    """The first force producer's wrench of a hull-plus-fins vehicle in still water, one row per command."""
    key = "force_producers.0.generated_force"
    return np.array([_evaluate(vessel, np.concatenate([np.zeros(6), v]), np.concatenate([np.zeros(3), c]))[1][key]
                     for c, v in zip(commands, nu_r)])


def test_the_rudder_and_stern_plane_as_lifting_fins_in_a_set_equal_the_fin_pairs_part_on_the_mss_rows():
    """The gate trees' REMUS 100 (one geometry for the whole vehicle): a set of two ``LiftingFin`` equals
    ``FinPairsDeflectionOnly`` on the commands and velocities of the frozen MSS rows."""
    pairs, fin_set = vc.build("remus100_fin_pairs"), vc.build("remus100_fin_set")
    assert vc.input_names(pairs)[-1] == "fin_deflection_command"
    assert vc.input_names(fin_set)[-1] == "force_producer_command"
    _, commands, nu_r, _ = _rows()
    assert max_diff(_forces(pairs, commands, nu_r), _forces(fin_set, commands, nu_r)) <= G2_TOLERANCE


def test_the_fin_set_equals_the_matlab_wrench_of_the_mss_rows_with_the_mss_fin_positions(tmp_path):
    """The same two trees with the fin positions of the MATLAB rows (``remus100.m`` computes them from its own
    hull length; the vehicle's consistent set uses ``-L/2``): the wrench equals ``tau`` of every row."""
    ref, commands, nu_r, expected = _rows()
    values = {**make_trees.REMUS, "rudder_position": mss_constant(ref, "x_r"),
              "stern_plane_position": mss_constant(ref, "x_s")}
    assert values["water_density"] == mss_constant(ref, "rho")
    pairs = make_trees.torpedo("MunkCouplingsRemoved", values=values, diagnostics=("force_producers.0.generated_force",),
                               producers=[make_trees.remus_part("FinPairsDeflectionOnly", {
                                   k: values[k] for k in ("rudder_area", "stern_plane_area", "rudder_lift_coefficient",
                                                          "stern_plane_lift_coefficient", "rudder_position",
                                                          "stern_plane_position", "max_deflection")})])
    fin_set = make_trees.torpedo("MunkCouplingsRemoved", values=values,
                                 diagnostics=("force_producers.0.generated_force",),
                                 producers=[make_trees.fin_pair_set(values)])
    for name, tree in (("pairs", pairs), ("set", fin_set)):
        vessel = vc.build_variant(tmp_path / name, tree)
        assert max_diff(_forces(vessel, commands, nu_r), expected) <= G1_TOLERANCE, name


def test_a_set_of_lag_fins_stacks_their_states_and_equals_the_model_layer_set(tmp_path):
    lag = make_trees.prestero_fin()
    second = make_trees.prestero_fin()
    second["params"].update({"fin_position": [-0.7, 0.0, 0.0], "lift_axis": [0.0, 0.0, -1.0]})
    tree = make_trees.torpedo("MunkCouplingsRemoved", diagnostics=("force_producers.0.generated_force",), producers=[
        make_trees.node("ForceProducerSet", {"command_count": 2, "command_map": [[1.0, 0.0], [0.0, 1.0]]}, [],
                        {"producers": [lag, second]})])
    vessel = vc.build_variant(tmp_path, tree)
    assert [(d.name, d.size) for d in vessel.graph().stateDescription] == [("producer_state", 2), ("state", 12)]
    forms = {"servo": "lag_rate_angle", "inflow": "rigid_point", "flow_angle": "small_angle",
             "interference": "none", "section": "quadratic_drag"}
    values = {"servo": SERVO_NUMBERS, "section": {"lift_slope": prestero_value("lift_slope")}}
    a = fin_from_forms(forms, values, geometry([prestero_value("fin_position_x"), 0, 0], E_Y,
                                               prestero_value("fin_area"), E_X))
    b = fin_from_forms(forms, values, geometry([-0.7, 0, 0], [0.0, 0.0, -1.0], prestero_value("fin_area"), E_X))
    reference = producer_set([a, b], np.eye(2))
    g = rng()
    worst = 0.0
    for _ in range(20):
        x = np.concatenate([g.uniform(-0.3, 0.3, 2), g.uniform(-0.2, 0.2, 6), g.uniform(-2.5, 2.5, 3),
                            g.uniform(-0.5, 0.5, 3)])
        command = g.uniform(-0.6, 0.6, 2)
        u = np.concatenate([np.zeros(3), command])
        xdot, out = _evaluate(vessel, x, u)
        ref = call(reference, command=command, state=x[:2], nu_r=x[2 + 6:2 + 12], water_density=RHO)
        worst = max(worst, max_diff(out["force_producers.0.generated_force"], ref["tau"]),
                    max_diff(xdot[:2], ref["state_dot"]))
    assert worst <= G2_TOLERANCE, worst


# --------------------------------------------------------------------------
# Refusals
# --------------------------------------------------------------------------
def _vehicle_with_fin(tmp_path, change):
    fin = make_trees.prestero_fin()
    change(fin)
    tree = make_trees.torpedo("MunkCouplingsRemoved", producers=[make_trees.node("PrescribedWrench"), fin])
    return vc.build_variant(tmp_path, tree)


def test_a_fin_refuses_a_part_that_does_not_fit_its_slot(tmp_path):
    with pytest.raises(Exception, match="slot servo"):
        _vehicle_with_fin(tmp_path, lambda fin: fin["children"].update(servo=make_trees.node("FinInflowTranslational")))


def test_a_fin_refuses_a_child_parameter_left_open(tmp_path):
    def open_slope(fin):
        fin["children"]["section"]["params"]["open_parameters"] = ["lift_slope"]
    with pytest.raises(Exception, match="open"):
        _vehicle_with_fin(tmp_path, open_slope)


@pytest.mark.parametrize("slot, plugin, source", [
    ("section", "FinSectionLiftingLine", "lifting-line induced drag"),
    ("interference", "FinInterferenceSlenderBody", "Pitts 1957"),
])
def test_a_fin_with_a_part_waiting_on_a_source_is_refused_saying_so(slot, plugin, source, tmp_path):
    with pytest.raises(Exception, match=f"waiting on source: {source}"):
        _vehicle_with_fin(tmp_path, lambda fin: fin["children"].update({slot: make_trees.node(plugin)}))


def test_a_set_refuses_a_command_map_of_the_wrong_shape(tmp_path):
    tree = make_trees.torpedo("MunkCouplingsRemoved", producers=[make_trees.node(
        "ForceProducerSet", {"command_count": 2, "command_map": [[1.0, 0.0, 0.0]]}, [],
        {"producers": [make_trees.prestero_fin(), make_trees.prestero_fin()]})])
    with pytest.raises(Exception, match="command_map"):
        vc.build_variant(tmp_path, tree)
