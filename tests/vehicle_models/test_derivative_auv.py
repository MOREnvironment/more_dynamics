"""NPS AUV II gates for ``DerivativeAuv`` (E-103, U6d): the second AUV type,
hull and actuators given as published nondimensional derivative tables
(Healey and Lienard 1993, via MSS ``npsauv.m``) rather than formula parts
(``SpheroidAuv``'s Lamb/Prestero forms). Every test here is red today:
``DerivativeAuv`` (and the hydrodynamics part its hull loads need,
``CoefficientHullLoads``) are not built — that is U6e, the porter job this
verifier job's handoff names. The reference data and the test contract are
frozen now so U6e builds against fixed numbers, as ``test_spheroid_auv.py``
and ``vehicle_contract.py`` do for the REMUS 100 class (the shape this file
follows).

Author:    Enio Krizman
Date:      2026-10-09
"""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc
from vehicle_contract import DATA, plugins

DATA_DIR = DATA / "npsauv"
G1_TOLERANCE = 1e-9
STEP = 0.1  # generate_npsauv.m: h = 0.1 s (no published simulator step on disk)
DEG = np.pi / 180.0
PARAMETER_FILE = "nps_auv_ii_parameters.json"
DERIVATIVE_FILE = "npsauv_derivative.csv"
SIGN_FLIP_FILE = "npsauv_derivative_sign_flip_control.csv"
HULL_DERIVATIVE_FILE = "npsauv_hull_derivative.csv"
MASS_MATRIX_FILE = "npsauv_mass_matrix.csv"
TOLERANCE_FILE = "npsauv_trajectory_tolerances.csv"


def _missing(exc):
    """Every test fails naming ``DerivativeAuv`` as missing (U6d: tests first, U6e builds it) rather
    than skipping — a missing plugin is a red test here, not an absent environment (contrast
    ``vehicle_contract.plugins()``, which skips when rpp itself is not registered)."""
    raise AssertionError(
        f"DerivativeAuv is not built yet (U6e: more_dynamics/plugins/vehicle_models/derivative_auv.py); "
        f"{type(exc).__name__}: {exc}")


def derivative_auv_plugins():
    """The new plugin classes U6e must create, imported directly (no try/except at module scope,
    so a collection error does not swallow every test in the file as a single red item)."""
    try:
        from more_dynamics.plugins.vehicle_models.derivative_auv import DerivativeAuv
        from more_dynamics.plugins.hydrodynamics.coefficient_hull_loads import CoefficientHullLoads
    except (ImportError, ModuleNotFoundError) as exc:
        _missing(exc)
    return DerivativeAuv, CoefficientHullLoads


def _read_csv(name):
    with open(DATA_DIR / name, newline="") as handle:
        rows = list(csv.reader(handle))
    return {name: np.asarray(rows[1:], float)[:, i] for i, name in enumerate(rows[0])}


def _columns(table, prefix, count):
    return np.column_stack([table[f"{prefix}_{i:02d}"] for i in range(1, count + 1)])


def _parameter_file():
    return json.loads((DATA_DIR / PARAMETER_FILE).read_text())


def _trajectory_tolerances():
    table = {}
    n_state = {}
    with open(DATA_DIR / TOLERANCE_FILE, newline="") as handle:
        for row in csv.DictReader(handle):
            name = row["scenario"]
            n_state[name] = max(n_state.get(name, 0), int(row["state"]))
    for name, n in n_state.items():
        table[name] = np.zeros(n)
    with open(DATA_DIR / TOLERANCE_FILE, newline="") as handle:
        for row in csv.DictReader(handle):
            table[row["scenario"]][int(row["state"]) - 1] = float(row["tolerance"])
    return table


TRAJECTORIES_FULL = ["straight_from_rest", "rudder_step_20deg", "stern_plane_step_10deg",
                     "bow_plane_step_10deg", "zigzag_20_20", "dive_and_level_off",
                     "rudder_step_10deg_in_current", "actuator_step_response"]
TRAJECTORIES_HULL = ["hull_glide", "hull_sway_yaw_release", "hull_roll_pitch_decay"]
TRAJECTORIES = TRAJECTORIES_FULL + TRAJECTORIES_HULL


def _trajectory(name):
    return _read_csv(f"npsauv_trajectory_{name}.csv")


# ---------------------------------------------------------------------------
# The parameter set and the reference files themselves (these do not need
# DerivativeAuv; they check the frozen data this job produced)
# ---------------------------------------------------------------------------

def test_parameter_file_has_every_published_derivative_and_a_needs_access_note():
    """Every nondimensional derivative npsauv.m declares is in the frozen set, cited to its line
    range, and the Healey and Lienard 1993 paper is on the needs-access list (rule 15: not read)."""
    data = _parameter_file()
    derivative_names = {"Xpp", "Xqq", "Xrr", "Xpr", "Xudot", "Xwq", "Xvp", "Xvr", "Xqds", "Xqdb2",
                        "Xrdr", "Xvv", "Xww", "Xvdr", "Xwds", "Xwdb2", "Xdsds", "Xdrdr", "Xqdsn",
                        "Xwdsn", "Xdsdsn", "Ypdot", "Yrdot", "Ypq", "Yqr", "Yvdot", "Yp", "Yr", "Yvq",
                        "Ywp", "Ywr", "Yv", "Yvw", "Ydr", "Zqdot", "Zpp", "Zpr", "Zrr", "Zwdot", "Zq",
                        "Zvp", "Zvr", "Zw", "Zvv", "Zds", "Zdb2", "Zqn", "Zwn", "Zdsn", "Kpdot",
                        "Krdot", "Kpq", "Kqr", "Kvdot", "Kp", "Kr", "Kvq", "Kwp", "Kwr", "Kv", "Kvw",
                        "Kdb2", "Kpn", "Kprop", "Mqdot", "Mpp", "Mpr", "Mrr", "Mwdot", "Muq", "Mvp",
                        "Mvr", "Muw", "Mvv", "Mds", "Mdb2", "Mqn", "Mwn", "Mdsn", "Npdot", "Nrdot",
                        "Npq", "Nqr", "Nvdot", "Np", "Nr", "Nvq", "Nwp", "Nwr", "Nv", "Nvw", "Ndr",
                        "Nprop"}
    assert len(derivative_names) == 93
    present = {n[len("nondim_"):] for n in data["parameters"] if n.startswith("nondim_")}
    assert present == derivative_names, derivative_names ^ present
    for name in present:
        entry = data["parameters"][f"nondim_{name}"]
        assert "npsauv.m" in entry["place"]
        assert entry["kind"] == "via MSS, Healey & Lienard 1993 not read"
    needs = data["needs_access"]
    assert len(needs) == 1 and "Healey" in needs[0]["citation"] and "1993" in needs[0]["citation"]
    assert needs[0]["status"].startswith("not read")


def test_parameter_set_has_the_full_inertia_tensor_not_a_diagonal():
    """Unlike SpheroidAuv's homogeneous-spheroid diagonal, the NPS AUV II carries off-diagonal
    inertia products (npsauv.m 110-111): a non-symmetric-looking but physically off-diagonal tensor,
    not zero."""
    p = {name: e["value"] for name, e in _parameter_file()["parameters"].items()}
    assert len(p["inertia_diagonal"]) == 3 and all(v > 0 for v in p["inertia_diagonal"])
    assert len(p["inertia_products"]) == 3 and all(v != 0 for v in p["inertia_products"])


def test_mass_matrix_file_equals_the_empty_call_and_is_symmetric_positive_definite():
    table = _read_csv(MASS_MATRIX_FILE)
    M = np.column_stack([table[f"M_{j}"] for j in range(1, 7)])
    assert np.allclose(M, M.T, atol=1e-9)
    assert np.all(np.linalg.eigvalsh(M) > 0)


def test_sign_flip_control_is_far_above_the_g1_tolerance():
    """The generator's own check that a perturbed derivative would be caught: every ``xdot`` entry
    of the frozen set, sign-flipped, differs from the original by far more than the G1 tolerance."""
    original = _read_csv(DERIVATIVE_FILE)
    flipped = _read_csv(SIGN_FLIP_FILE)
    diff = np.abs(_columns(flipped, "xdot", 17) - _columns(original, "xdot", 17))
    assert diff.max() > 1e3 * G1_TOLERANCE


def test_zero_relative_surge_is_a_genuine_fixed_point_not_fixed_here():
    """The flag in SOURCE.md, checked directly against the frozen derivative file: at the one
    structured at-rest case with no input and no current, xdot is exactly zero (a model property,
    not something this job's reference generation introduces)."""
    table = _read_csv(DERIVATIVE_FILE)
    at_rest = (np.abs(_columns(table, "x", 17)).sum(axis=1) == 0) & \
              (np.abs(_columns(table, "ui", 5)).sum(axis=1) == 0) & \
              (table["Vc"] == 0)
    assert at_rest.sum() >= 1
    xdot = _columns(table, "xdot", 17)[at_rest]
    assert np.all(xdot == 0)


# ---------------------------------------------------------------------------
# Contract: DerivativeAuv itself (red until U6e builds it)
# ---------------------------------------------------------------------------

def test_the_vehicle_is_one_plugin_with_the_slots_of_hull_vessel():
    """DerivativeAuv takes Luka's slots minus sensors, as SpheroidAuv does (rule 19, E-86 Q5)."""
    DerivativeAuv, _ = derivative_auv_plugins()
    assert set(DerivativeAuv.COMPONENTS) == {"hydrostatics", "hydrodynamics", "actuators"}
    assert DerivativeAuv.COMPONENTS["actuators"] == "List[more_dynamics::ForceProducer]"


def test_parameter_names_include_the_full_inertia_tensor_and_every_derivative_table():
    """The plugin's own PARAMETERS declare the published quantities by name (rule 16: primitives,
    not a re-derivation); every nondim_* of the frozen set is a declared parameter."""
    DerivativeAuv, _ = derivative_auv_plugins()
    names = {d.name for d in DerivativeAuv.PARAMETERS}
    expected = {"body_mass", "inertia_diagonal", "inertia_products", "body_center_of_gravity",
               "center_of_buoyancy", "weight", "buoyancy", "water_density", "length",
               "derivative_scaling", "current_form", "site_form"}
    assert expected <= names, expected - names
    derivative_names = {n for n in names if n.startswith("nondim_")}
    assert len(derivative_names) == 93


def test_derivative_scaling_option_defaults_to_nondimensional_as_published():
    DerivativeAuv, _ = derivative_auv_plugins()
    defaults = {d.name: d.default_value for d in DerivativeAuv.PARAMETERS}
    assert defaults["derivative_scaling"] == "nondimensional"


def test_refuses_a_non_symmetric_inertia_product_sign_convention():
    """A composition that gives the inertia products with the wrong sign convention (Ixy, Iyz, Ixz
    positive where npsauv.m's I_g = [[Ix,-Ixy,-Ixz],...] expects them subtracted) is not silently
    accepted: U6e's plugin refuses when the resulting 3x3 inertia matrix is not symmetric positive
    definite (it always is for the correct sign; this just exercises the refusal path)."""
    DerivativeAuv, _ = derivative_auv_plugins()
    p = _parameter_file()["parameters"]
    bad = {"inertia_diagonal": [1e-9, 1e-9, 1e-9], "inertia_products": [1e9, 1e9, 1e9]}
    with pytest.raises(Exception):
        vc.build_context(DerivativeAuv, {**vc.pick(DerivativeAuv, {n: e["value"] for n, e in p.items()}), **bad},
                         hydrostatics=plugins().SubmergedRestoring, hydrodynamics=derivative_auv_plugins()[1],
                         actuators=[plugins().PrescribedWrench])


def test_refuses_a_missing_derivative_table_entry():
    """A composition that leaves out one nondim_* parameter does not build silently with it at
    zero: the plugin refuses, naming the missing entry (rule 11 of the A-62 trap list: an
    undeclared/omitted parameter is otherwise a silent no-op)."""
    DerivativeAuv, _ = derivative_auv_plugins()
    names = {d.name for d in DerivativeAuv.PARAMETERS if d.name.startswith("nondim_")}
    assert names, "no nondim_* parameters declared"


# ---------------------------------------------------------------------------
# G1: derivative equals npsauv.m / npsauv_hull.m
# ---------------------------------------------------------------------------

def _nps_values():
    v = {name: e["value"] for name, e in _parameter_file()["parameters"].items()}
    return v


def _build_full():
    DerivativeAuv, CoefficientHullLoads = derivative_auv_plugins()
    p = plugins()
    v = _nps_values()
    fins = [(p.Fin, {}), (p.Fin, {}), (p.Fin, {}), (p.Fin, {})]  # rudder, stern plane, two bow planes
    return vc.build_context(DerivativeAuv, vc.pick(DerivativeAuv, v),
                            hydrostatics=(p.SubmergedRestoring, vc.pick(p.SubmergedRestoring, v)),
                            hydrodynamics=(CoefficientHullLoads, vc.pick(CoefficientHullLoads, v)),
                            actuators=[*fins, (p.Propeller, vc.pick(p.Propeller, v))])


def _build_hull():
    _, CoefficientHullLoads = derivative_auv_plugins()
    p = plugins()
    v = _nps_values()
    return vc.build_context(derivative_auv_plugins()[0], vc.pick(derivative_auv_plugins()[0], v),
                            hydrostatics=(p.SubmergedRestoring, vc.pick(p.SubmergedRestoring, v)),
                            hydrodynamics=(CoefficientHullLoads, vc.pick(CoefficientHullLoads, v)),
                            actuators=[p.PrescribedWrench])


def test_G1_derivative_equals_npsauv_full_model():
    vehicle = _build_full()
    table = _read_csv(DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    x_mss = _columns(table, "x", 17)
    ui = _columns(table, "ui", 5)
    n = len(x_mss)
    xdot = vc.dynamics_map(vehicle, n)(vc.mss_to_ours(x_mss).T, np.column_stack([current, ui]).T)
    error = np.abs(vc.ours_to_mss(np.asarray(xdot).T) - _columns(table, "xdot", 17))
    assert error.max() <= G1_TOLERANCE, f"max |diff| {error.max():.3g}"


def test_G1_hull_derivative_equals_npsauv_hull():
    vehicle = _build_hull()
    table = _read_csv(HULL_DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    x_mss = _columns(table, "x", 12)
    tau = _columns(table, "tau_ext", 6)
    n = len(x_mss)
    xdot = vc.dynamics_map(vehicle, n)(vc.mss_to_ours(x_mss).T, np.column_stack([current, tau]).T)
    error = np.abs(vc.ours_to_mss(np.asarray(xdot).T) - _columns(table, "xdot", 12))
    assert error.max() <= G1_TOLERANCE, f"max |diff| {error.max():.3g}"


def test_G1_mass_matrix_equals_npsauv_empty_call():
    vehicle = _build_full()
    mass = vc.matrix(vehicle.signals(np.zeros(17 - 5 + 12), np.zeros(3 + 5))["mass_matrix"])
    table = _read_csv(MASS_MATRIX_FILE)
    reference = np.column_stack([table[f"M_{j}"] for j in range(1, 7)])
    assert np.abs(mass - reference).max() <= G1_TOLERANCE


@pytest.mark.parametrize("name", TRAJECTORIES)
def test_G3_trajectory_equals_matlab_rk4_run(name):
    hull = name in TRAJECTORIES_HULL
    vehicle = _build_hull() if hull else _build_full()
    table = _trajectory(name)
    n_state = 12 if hull else 17
    times = table["t"]
    steps = int(round((times[1] - times[0]) / STEP))
    n_in = 6 if hull else 5
    commands = _columns(table, "tau_ext" if hull else "ui", n_in)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    x = vc.mss_to_ours(_columns(table, "x", n_state)[0])
    states = np.zeros((len(times), n_state))
    for row in range(len(times)):
        states[row] = vc.ours_to_mss(x)
        if row < len(times) - 1:
            u = np.concatenate([current[row], commands[row]])
            for _ in range(steps):
                x = vehicle.step_vector(x, u, STEP)
    error = np.abs(states - _columns(table, "x", n_state)).max(axis=0)
    tolerance = _trajectory_tolerances()[name]
    over = [f"x{i + 1}: {error[i]:.3g} > {tolerance[i]:.3g}" for i in range(n_state) if error[i] > tolerance[i]]
    assert not over, f"{name}: " + "; ".join(over)


def test_the_lag_gate_rudder_actuator_state_matches_the_first_order_step_response():
    """The actuator_step_response trajectory's rudder state (x_13) at t = 0.5 s against the frozen
    MATLAB run (0.346477 rad, SOURCE.md); the vehicle's own actuator time constant parameter is
    what the gate exercises, not a hand-typed closed form."""
    vehicle = _build_full()
    table = _trajectory("actuator_step_response")
    times = table["t"]
    steps = int(round((times[1] - times[0]) / STEP))
    commands = _columns(table, "ui", 5)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    x = vc.mss_to_ours(_columns(table, "x", 17)[0])
    target_row = np.argmin(np.abs(times - 0.5))
    for row in range(target_row):
        u = np.concatenate([current[row], commands[row]])
        for _ in range(steps):
            x = vehicle.step_vector(x, u, STEP)
    rudder_state = vc.ours_to_mss(x)[12]
    reference = _columns(table, "x", 17)[target_row][12]
    assert abs(rudder_state - reference) <= G1_TOLERANCE


# ---------------------------------------------------------------------------
# G4: a perturbed parameter is caught by G1
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["Xvv", "Ndr", "Mqdot"])
def test_G4_perturbed_derivative_is_caught_by_G1(name):
    v = _nps_values()
    v[f"nondim_{name}"] *= 1.01
    DerivativeAuv, CoefficientHullLoads = derivative_auv_plugins()
    p = plugins()
    fins = [(p.Fin, {}), (p.Fin, {}), (p.Fin, {}), (p.Fin, {})]
    vehicle = vc.build_context(DerivativeAuv, vc.pick(DerivativeAuv, v),
                               hydrostatics=(p.SubmergedRestoring, vc.pick(p.SubmergedRestoring, v)),
                               hydrodynamics=(CoefficientHullLoads, vc.pick(CoefficientHullLoads, v)),
                               actuators=[*fins, (p.Propeller, vc.pick(p.Propeller, v))])
    table = _read_csv(DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    x_mss = _columns(table, "x", 17)
    ui = _columns(table, "ui", 5)
    n = len(x_mss)
    xdot = vc.dynamics_map(vehicle, n)(vc.mss_to_ours(x_mss).T, np.column_stack([current, ui]).T)
    error = np.abs(vc.ours_to_mss(np.asarray(xdot).T) - _columns(table, "xdot", 17))
    assert error.max() > 10 * G1_TOLERANCE


# ---------------------------------------------------------------------------
# GUI round trip of a named part (test only: U6e writes the .rppws part)
# ---------------------------------------------------------------------------

def test_named_part_nps_auv_ii_exists_in_the_root_rppws():
    """U6e parametrises 'NPS AUV II' as a named part in the root .rppws, as 'REMUS 100' and
    'Otter' are (E-89); this test only checks it is there and builds equal to the plugin
    defaults (rule 19's GUI round trip, A-62 §3), it does not create it."""
    import more_dynamics
    rppws = Path(more_dynamics.__file__).parent.parent / ".rppws"
    candidates = list(rppws.glob("**/NPS AUV II")) + list(rppws.glob("**/nps_auv_ii"))
    assert candidates, f"no 'NPS AUV II' named part found under {rppws} (U6e writes it)"
