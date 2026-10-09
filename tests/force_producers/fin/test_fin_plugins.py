"""The fin, the fin pairs and the propeller as plugins of Luka's form
(``ForceProducer``), built with rpp's own ``ComponentContext`` inside a REMUS
100 vehicle.

* Every force-producer plugin of ours is listed in ``plugins.json``, one class
  per file, every parameter with a default; the defaults are the printed
  numbers they cite; the two forms waiting on a source refuse to initialise
  saying so.
* The ``Fin`` plugin inside a vehicle: its servo state precedes the vehicle's
  twelve states, its command is the vehicle's last input, its wrench equals the
  model-layer ``lifting_fin_casadi`` assembled from the same forms and
  numbers, and its servo derivative equals the servo model's.
* The servo's four settings (the lag, the rate limit and the angle limit on or
  off) equal the servo model with the same selectors; settings the model
  refuses are refused when the vehicle is built.
* The REMUS 100 rudder and stern plane as two ``Fin`` plugins equal the landed
  ``FinPairsDeflectionOnly`` plugin on the frozen MSS rows (and the MATLAB
  wrench of those rows).
* A fin refuses an unknown form and a form waiting on a source; the plugin
  layer holds no trace of the retired part base.

Reference numbers: the frozen MSS rows and the printed Prestero values of
``tests/force_producers/fin/data`` and the REMUS parameter file of
``tests/vehicle_models/data``.

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

LIBRARY = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vehicle_models"))
import vehicle_contract as vc  # noqa: E402

FORCE_PRODUCER_PLUGINS = {"Fin": "force_producers/fin", "FinPairsDeflectionOnly": "force_producers/fin_pairs_deflection_only",
                          "Propeller": "force_producers/propeller", "PrescribedWrench": "force_producers/prescribed_wrench"}
DEG = np.pi / 180.0
# the servo numbers of the plugin defaults (Sarhadi 2026, Fig. 4, p. 4)
SERVO_NUMBERS = {"time_constant": 0.1, "max_deflection": 20.0 * DEG, "max_rate": 30.0 * DEG}
RHO = vc.REMUS["water_density"]


def _classes():
    vc.plugins()
    registry = json.loads((LIBRARY / "plugins.json").read_text())
    entries = {e["Name"]: e["Path"] for e in registry["Plugins"]}
    found = {}
    for name, path in FORCE_PRODUCER_PLUGINS.items():
        assert entries.get(name) == f"more_dynamics/plugins/{path}.py", (name, entries.get(name))
        module = importlib.import_module(entries[name][:-3].replace("/", "."))
        found[name] = getattr(module, name)
    return found


# --------------------------------------------------------------------------
# The plugins
# --------------------------------------------------------------------------
def test_every_force_producer_plugin_is_listed_and_has_a_class_and_one_class_per_file():
    classes = _classes()
    assert set(classes) == set(FORCE_PRODUCER_PLUGINS)
    for name, path in FORCE_PRODUCER_PLUGINS.items():
        text = (LIBRARY / "more_dynamics" / "plugins" / f"{path}.py").read_text()
        assert text.count("\nclass ") == 1, name


def test_every_parameter_has_a_default():
    for name, cls in _classes().items():
        for description in cls.PARAMETERS:
            assert description.default_value is not None, (name, description.name)


def _default(cls, name):
    return next(d.default_value for d in cls.PARAMETERS if d.name == name)


def test_defaults_are_the_printed_numbers_they_cite():
    fin = _classes()["Fin"]
    assert _default(fin, "fin_position") == [prestero_value("fin_position_x"), 0.0, 0.0]  # Table A.5, p. 103
    assert _default(fin, "fin_area") == prestero_value("fin_area")
    assert _default(fin, "lift_slope") == prestero_value("lift_slope")
    for name, value in SERVO_NUMBERS.items():  # Sarhadi 2026, Fig. 4, p. 4
        assert abs(_default(fin, name) - value) <= 1e-15, name
    assert [_default(fin, f"servo_{n}") for n in ("dynamics", "rate_limit", "angle_limit")] == \
        [SERVO_SETTINGS["lag_rate_angle"][n] for n in ("dynamics", "rate_limit", "angle_limit")]


@pytest.mark.parametrize("slot, form, source", [
    ("section", "lifting_line", "lifting-line induced drag"),
    ("interference", "slender_body", "Pitts 1957"),
])
def test_a_form_waiting_on_a_source_refuses_to_initialise_saying_so(slot, form, source):
    fin = _classes()["Fin"]
    with pytest.raises(NotImplementedError, match=f"waits on source: {source}"):
        vc.context_of(fin, {slot: form}).initialize()


@pytest.mark.parametrize("slot, form", [("inflow", "wake"), ("flow_angle", "large_angle"), ("section", "stall"),
                                        ("interference", "none_such")])
def test_a_fin_refuses_an_unknown_form_and_names_the_choices(slot, form):
    fin = _classes()["Fin"]
    with pytest.raises(ValueError, match=f"Fin: {slot} must be one of"):
        vc.context_of(fin, {slot: form}).initialize()


def test_the_plugin_layer_holds_no_trace_of_the_retired_part_base():
    offenders = []
    for path in (LIBRARY / "more_dynamics" / "plugins").rglob("*.py"):
        text = path.read_text()
        if any(word in text for word in ("BlockPayload", ".block()", "def block(", "fin_part_plugin",
                                         "proportional_force_fin", "ProportionalForceFin")):
            offenders.append(str(path.relative_to(LIBRARY)))
    assert not offenders, offenders
    plugins = LIBRARY / "more_dynamics" / "plugins"
    assert not (plugins / "vehicles").exists() and not (plugins / "force_producers" / "fin").exists()


# --------------------------------------------------------------------------
# Vehicle helpers: the REMUS hull with a wrench and one Fin
# --------------------------------------------------------------------------
def _hull_with_fin(fin_params=None):
    p = vc.plugins()
    return vc.remus(coriolis="munk_couplings_removed", current="none",
                    actuators=[p.PrescribedWrench, (p.Fin, fin_params or {})])


def _evaluate(vessel, x, u):
    """The vehicle's state derivative and its named signals at stacked ``x`` and ``u``."""
    graph = vc.graph_of(vessel)
    return np.asarray(graph.step(x, u)).ravel(), vessel.signals(x, u)


def _random_states(n, child=0):
    g = rng()
    for _ in range(n):
        yield np.concatenate([g.uniform(-0.3, 0.3, child), g.uniform(-0.2, 0.2, 6), g.uniform(-2.5, 2.5, 3),
                              g.uniform(-0.5, 0.5, 3)]), g.uniform(-0.6, 0.6)


def _fin_input(command):
    """Inputs of a hull-plus-one-fin vehicle: the wrench, the fin command."""
    return np.concatenate([np.zeros(6), [command]])


def _reference_fin(servo="lag_rate_angle"):
    """The model-layer fin of the plugin defaults: the lag servo, rigid-point inflow, small-angle flow angle, no
    interference, quadratic-drag section."""
    declared = {"max_deflection", "max_rate", "time_constant"}
    return fin_from_forms(
        {"servo": servo, "inflow": "rigid_point", "flow_angle": "small_angle", "interference": "none",
         "section": "quadratic_drag"},
        {"servo": {k: v for k, v in SERVO_NUMBERS.items() if k in declared},
         "section": {"lift_slope": prestero_value("lift_slope")}},
        geometry([prestero_value("fin_position_x"), 0.0, 0.0], E_Y, prestero_value("fin_area"), E_X))


# --------------------------------------------------------------------------
# The fin in a vehicle
# --------------------------------------------------------------------------
def test_a_fin_is_an_actuator_with_its_servo_state_before_the_vehicle_states():
    vessel = _hull_with_fin()
    graph = vessel.graph()
    assert [(d.name, d.size) for d in graph.stateDescription] == [("servo_state", 1), ("state", 12)]
    assert vc.input_names(vessel) == ["desired_wrench", "fin_angle_command"]


def test_the_fin_wrench_in_the_vehicle_is_the_model_layer_fin_on_the_same_forms():
    vessel = _hull_with_fin()
    fin = _reference_fin()
    worst = 0.0
    for x, command in _random_states(20, child=1):
        _, out = _evaluate(vessel, x, _fin_input(command))
        reference = tau(fin, command, x[1 + 6:1 + 12], RHO, state=x[:1])  # nu_r = nu in still water
        worst = max(worst, max_diff(out["actuators.1.generated_force"], reference))
    assert worst <= G2_TOLERANCE, worst


def test_the_servo_derivative_in_the_vehicle_is_the_servo_models():
    vessel = _hull_with_fin()
    servo = part("servo", "lag_rate_angle", SERVO_NUMBERS)
    worst = 0.0
    for x, command in _random_states(20, child=1):
        xdot, out = _evaluate(vessel, x, _fin_input(command))
        reference = call(servo, command=command, servo_state=x[:1])
        worst = max(worst, abs(xdot[0] - reference["servo_state_dot"][0]),
                    abs(out["actuators.1.deflection"][0] - reference["deflection"][0]))
    assert worst <= 1e-14, worst


# --------------------------------------------------------------------------
# The servo's switches
# --------------------------------------------------------------------------
def _servo_settings(settings):
    return {f"servo_{key}": value for key, value in settings.items()}


@pytest.mark.parametrize("name", sorted(SERVO_SETTINGS))
def test_each_servo_setting_equals_the_servo_model_with_the_same_selectors(name):
    settings = SERVO_SETTINGS[name]
    vessel = _hull_with_fin(_servo_settings(settings))
    lag = settings["dynamics"] == "first_order_lag"
    declared = {d.name for d in part_declared("servo", name)}
    servo = part("servo", name, {k: v for k, v in SERVO_NUMBERS.items() if k in declared})
    assert [d.name for d in vessel.graph().stateDescription] == (["servo_state", "state"] if lag else ["state"])
    worst = 0.0
    for x, command in _random_states(20, child=1 if lag else 0):
        command *= 2.0  # beyond the 20 degree limit now and then
        xdot, out = _evaluate(vessel, x, _fin_input(command))
        reference = call(servo, command=command, servo_state=x[:1] if lag else [])
        worst = max(worst, abs(out["actuators.1.deflection"][0] - reference["deflection"][0]))
        if lag:
            worst = max(worst, abs(xdot[0] - reference["servo_state_dot"][0]))
    assert worst <= 1e-14, worst


@pytest.mark.parametrize("settings, word", [
    ({"dynamics": "none", "rate_limit": True}, "rate_limit"),
    ({"dynamics": "none", "rate_limit": False, "angle_limit": "on_output"}, "angle_limit"),
    ({"angle_limit": "on_both"}, "angle_limit"),
])
def test_a_servo_setting_the_model_refuses_is_refused_when_the_vehicle_is_built(settings, word):
    with pytest.raises(Exception, match=word):
        _hull_with_fin(_servo_settings(settings))


# --------------------------------------------------------------------------
# Two fins against the landed fin pairs, on the frozen MSS rows
# --------------------------------------------------------------------------
def _rows():
    ref = mss_fins()
    return ref, np.column_stack([ref["ui1"], ref["ui2"]]), mss_columns(ref, "nu_r", 6), mss_columns(ref, "tau", 6)


def _two_ideal_fins(values):
    """The rudder and the stern plane of REMUS 100 as two ``Fin`` plugins (ideal servo, translational inflow, no
    flow angle, no interference, quadratic-drag section): the deflection-only fin pairs of ``remus100.m`` 228-245."""
    def fin(prefix, lift_axis):
        return {"fin_position": [values[f"{prefix}_position"], 0.0, 0.0], "chord_axis": [1.0, 0.0, 0.0],
                "lift_axis": lift_axis, "fin_area": values[f"{prefix}_area"], "servo_dynamics": "none",
                "servo_rate_limit": False, "servo_angle_limit": "on_command",
                "max_deflection": values["max_deflection"], "inflow": "translational", "flow_angle": "none",
                "interference": "none", "section": "quadratic_drag",
                "lift_slope": values[f"{prefix}_lift_coefficient"]}
    p = vc.plugins()
    return [(p.Fin, fin("rudder", [0.0, -1.0, 0.0])), (p.Fin, fin("stern_plane", [0.0, 0.0, -1.0]))]


def _pairs(values, **choices):
    p = vc.plugins()
    return [(p.FinPairsDeflectionOnly, {**vc.pick(p.FinPairsDeflectionOnly, values), **choices})]


def _wrench(vessel, commands, nu_r, keys):
    """The summed wrench of the named actuators of a hull-plus-fins vehicle in still water, one row per command."""
    rows = []
    for c, v in zip(commands, nu_r):
        signals = vessel.signals(np.concatenate([np.zeros(6), v]), c)
        rows.append(sum(signals[key] for key in keys))
    return np.array(rows)


def test_the_rudder_and_stern_plane_as_two_fins_equal_the_fin_pairs_plugin_on_the_mss_rows():
    values = vc.remus_values()
    pairs = vc.remus(coriolis="munk_couplings_removed", current="none", actuators=_pairs(values))
    fins = vc.remus(coriolis="munk_couplings_removed", current="none", actuators=_two_ideal_fins(values))
    assert vc.input_names(pairs) == ["fin_deflection_command"] and vc.input_names(fins) == [
        "fin_angle_command", "fin_angle_command"]
    _, commands, nu_r, _ = _rows()
    a = _wrench(pairs, commands, nu_r, ["actuators.0.generated_force"])
    b = _wrench(fins, commands, nu_r, ["actuators.0.generated_force", "actuators.1.generated_force"])
    assert max_diff(a, b) <= G2_TOLERANCE


def test_the_two_fins_equal_the_matlab_wrench_of_the_mss_rows_with_the_mss_fin_positions():
    """The same compositions with the fin positions of the MATLAB rows (``remus100.m`` computes them from its own
    hull length; the vehicle's consistent set uses ``-L/2``): the wrench equals ``tau`` of every row."""
    ref, commands, nu_r, expected = _rows()
    values = {**vc.remus_values(), "rudder_position": mss_constant(ref, "x_r"),
              "stern_plane_position": mss_constant(ref, "x_s")}
    assert values["water_density"] == mss_constant(ref, "rho")
    pairs = vc.remus(coriolis="munk_couplings_removed", current="none", values=values,
                     actuators=_pairs(values, fin_positions_method="given"))  # the rows' own positions, not -a
    fins = vc.remus(coriolis="munk_couplings_removed", current="none", values=values,
                    actuators=_two_ideal_fins(values))
    assert max_diff(_wrench(pairs, commands, nu_r, ["actuators.0.generated_force"]), expected) <= G1_TOLERANCE
    assert max_diff(_wrench(fins, commands, nu_r, ["actuators.0.generated_force", "actuators.1.generated_force"]),
                    expected) <= G1_TOLERANCE


def test_two_lag_fins_stack_their_states_and_equal_the_model_layer_set():
    p = vc.plugins()
    first = {}
    second = {"fin_position": [-0.7, 0.0, 0.0], "lift_axis": [0.0, 0.0, -1.0]}
    vessel = vc.remus(coriolis="munk_couplings_removed", current="none", actuators=[(p.Fin, first), (p.Fin, second)])
    assert [(d.name, d.size) for d in vessel.graph().stateDescription] == [
        ("servo_state", 1), ("servo_state", 1), ("state", 12)]
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
        xdot, out = _evaluate(vessel, x, command)
        ref = call(reference, command=command, state=x[:2], nu_r=x[2 + 6:2 + 12], water_density=RHO)
        worst = max(worst, max_diff(out["actuators.0.generated_force"] + out["actuators.1.generated_force"], ref["tau"]),
                    max_diff(xdot[:2], ref["state_dot"]))
    assert worst <= G2_TOLERANCE, worst
