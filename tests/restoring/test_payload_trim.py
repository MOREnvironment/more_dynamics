"""Gate tests for the payload trim of the surface restoring (``payload_trim_form``).

A point payload off the centre of flotation shifts the equilibrium of the linear restoring, ``g = G (eta - eta_0)``
(``models/restoring/payload_trim.py``). ``z = 0`` is the unloaded waterline, so the payload weight is a wrench on it
(heave and moment about the origin). Forms: ``weight_wrench`` (the default; that wrench, kinematics at the actual
attitude), ``lever_arm_and_heave_shift`` (MSS ``otter_legacy.m`` 264-275 as it is: the same wrench, the kinematics at
the shifted attitude) and ``lever_arm_moment`` (legacy, origin-dependent: no heave term; PT1, E-211 withdrew it as the
default, ``test_payload_trim_origin.py`` shows why).

Gates (the test names carry them): G1 the vehicle against MATLAB running MSS (frozen CSV
``data/payload_trim_mss_aff76cf.csv`` for the physics form, ``data/payload_trim_mss_6e78e98_beam_strips.csv`` the earlier record for the
MSS form, ``data/generate_payload_trim_mss.m``): the shift, the restoring vector and the
whole derivative of the Otter with two thrusters, for the MSS form against the unmodified ``otter_legacy.m`` and for the
physics form against the same call with the kinematics at the actual attitude (written out in the generator; equal to
``otter_legacy.m`` at d78d97e, which repaired that line, to 3e-15), 1e-9; G2 properties the equations imply (the heave
balance, the shift shared by both forms, the loaded displacement); G4 a perturbed payload or the other form is detected.

Hidden assumptions: the physics form is a second implementation written from ``otter_legacy.m``'s own workspace and
checked against its later revision, not a published result (its physics is the argument of
``models/restoring/payload_trim.py``); the reference is MSS, not the boat; the shifted-attitude kinematics of the MSS
form need the vehicle's ``kinematic_attitude_offset``.

Author:    Enio Krizman
Date:      2026-10-26
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "vehicle_models"))
import vehicle_contract as vc  # noqa: E402

DATA = Path(__file__).resolve().parent / "data"
G1_TOLERANCE = 1e-9
FORMS = ("lever_arm_and_heave_shift", "weight_wrench")
SUFFIX = {"lever_arm_and_heave_shift": "mss", "weight_wrench": "corr"}


def _read(name):
    data = np.genfromtxt(DATA / name, delimiter=",", names=True)
    return {name: np.atleast_1d(data[name]) for name in data.dtype.names}


def _table():
    """The earlier record of the unmodified function (kinematics at the shifted attitude): the history the
    ``lever_arm_and_heave_shift`` form follows. It is not re-made: that form reproduces a state of the source model
    that its author repaired."""
    return _read("payload_trim_mss_6e78e98_beam_strips.csv")


def _table_current():
    """The unmodified function at the pinned revision (kinematics at the actual attitude)."""
    return _read("payload_trim_mss_aff76cf.csv")


def _cols(table, prefix, count):
    return np.column_stack([table[f"{prefix}_{i:02d}"] for i in range(1, count + 1)])


def _otter(form, payload_mass, payload_position, **options):
    """The Otter on the frozen MSS set with two thrusters (the named parts' values), the current of ``otter.m``
    113-118, and the payload trim ``form``."""
    p = vc.plugins()
    thrusters = [(p.Propulsor, vc.pick(p.Propulsor, vc.reference_part(name)))
                 for name in ("Otter port thruster", "Otter starboard thruster")]
    hydrostatics = (p.SurfaceRestoring, {**vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2,
                                         "payload_trim_form": form, **options})
    return vc.otter(current="full_rotation_rate", hydrostatics=hydrostatics, actuators=thrusters,
                    vehicle_params={"payload_mass": float(payload_mass), "payload_position": list(payload_position)})


def _cases(table):
    for k in range(len(table["case_id"])):
        yield k, float(table["mp"][k]), _cols(table, "rp", 3)[k]


def _evaluate(table, form):
    """Per case: the vehicle's shift, restoring vector and derivative at the case's state, commands and current."""
    worst = {"eta_0": 0.0, "restoring": 0.0, "xdot": 0.0}
    suffix = SUFFIX[form]
    reference = {"lever_arm_and_heave_shift": ("eta0", "g_mss", "xdot_mss"),
                 "weight_wrench": ("eta0c", "g_corr", "xdot_corr")}[form]
    vehicles = {}
    for k, mp, rp in _cases(table):
        key = (mp, tuple(rp))
        if key not in vehicles:
            vehicles[key] = _otter(form, mp, rp)
        vehicle = vehicles[key]
        x = _cols(table, "x", 12)[k]
        u = np.r_[table["Vc"][k], table["beta_c"][k], 0.0, _cols(table, "n", 2)[k]]
        ours = vc.mss_to_ours(x)
        signals = vehicle.signals(ours, u)
        eta_0 = np.asarray(signals["payload_trim_offset"]).ravel()
        restoring = -np.asarray(signals["hydrostatics.restoring_force"]).ravel()
        xdot = vc.ours_to_mss(np.asarray(vc.graph_of(vehicle).step(ours, u)).ravel())
        worst["eta_0"] = max(worst["eta_0"], np.abs(eta_0 - _cols(table, reference[0], 6)[k]).max())
        worst["restoring"] = max(worst["restoring"], np.abs(restoring - _cols(table, reference[1], 6)[k]).max())
        worst["xdot"] = max(worst["xdot"], np.abs(xdot - _cols(table, reference[2], 12)[k]).max())
    return worst


@pytest.mark.parametrize("form", FORMS)
def test_G1_the_vehicle_equals_matlab_in_shift_restoring_and_derivative(form):
    """The shift ``eta_0``, the restoring vector and the whole derivative of the Otter at 48 cases (payload masses 0,
    10, 25, 50 kg at two positions, seeded states and commands and a current), against the unmodified ``otter.m`` for
    the MSS form and against the corrected form of the generator."""
    worst = _evaluate(_table() if form == "lever_arm_and_heave_shift" else _table_current(), form)
    assert max(worst.values()) <= G1_TOLERANCE, worst


def test_G1_the_physics_form_equals_the_unmodified_function_at_the_pinned_revision():
    """The source model repaired the kinematics line: its unmodified ``otter_legacy.m`` now equals the physics form
    (``weight_wrench``), in the frozen file (the generator's corrected column against its own unmodified column) and in
    the vehicle against that column."""
    table = _table_current()
    for i in range(1, 13):
        assert np.abs(table[f"xdot_mss_{i:02d}"] - table[f"xdot_corr_{i:02d}"]).max() <= 1e-13
    worst = 0.0
    vehicles = {}
    for k, mp, rp in _cases(table):
        key = (mp, tuple(rp))
        if key not in vehicles:
            vehicles[key] = _otter("weight_wrench", mp, rp)
        x = vc.mss_to_ours(_cols(table, "x", 12)[k])
        u = np.r_[table["Vc"][k], table["beta_c"][k], 0.0, _cols(table, "n", 2)[k]]
        xdot = vc.ours_to_mss(np.asarray(vc.graph_of(vehicles[key]).step(x, u)).ravel())
        worst = max(worst, np.abs(xdot - _cols(table, "xdot_mss", 12)[k]).max())
    assert worst <= G1_TOLERANCE, worst


def test_G4_each_form_is_told_from_the_other_by_the_reference():
    """A perturbed model is detected: the physics form misses the MSS derivative (its kinematics are at the actual
    attitude) and the legacy form's shift misses the reference ``eta_0`` (no heave term)."""
    table = _table()
    mss = _evaluate(table, "lever_arm_and_heave_shift")
    physics_against_mss = {"eta_0": 0.0, "xdot": 0.0}
    legacy_against_mss = 0.0
    vehicles = {}
    for k, mp, rp in _cases(table):
        if mp == 0.0:
            continue
        key = (mp, tuple(rp))
        vehicles.setdefault(key, (_otter("weight_wrench", mp, rp), _otter("lever_arm_moment", mp, rp)))
        physics, legacy = vehicles[key]
        x = vc.mss_to_ours(_cols(table, "x", 12)[k])
        u = np.r_[table["Vc"][k], table["beta_c"][k], 0.0, _cols(table, "n", 2)[k]]
        reference = _cols(table, "eta0", 6)[k]
        physics_against_mss["eta_0"] = max(physics_against_mss["eta_0"], np.abs(
            np.asarray(physics.signals(x, u)["payload_trim_offset"]).ravel() - reference).max())
        legacy_against_mss = max(legacy_against_mss, np.abs(
            np.asarray(legacy.signals(x, u)["payload_trim_offset"]).ravel() - reference).max())
        xdot = vc.ours_to_mss(np.asarray(vc.graph_of(physics).step(x, u)).ravel())
        physics_against_mss["xdot"] = max(physics_against_mss["xdot"], np.abs(
            xdot - _cols(table, "xdot_mss", 12)[k]).max())
    assert max(mss.values()) <= G1_TOLERANCE
    assert physics_against_mss["eta_0"] <= G1_TOLERANCE, physics_against_mss  # the same shift in both forms
    assert physics_against_mss["xdot"] > 1e-3, physics_against_mss  # the kinematics are the difference
    assert legacy_against_mss > 1e-3, legacy_against_mss  # the heave shift is missing


def test_G2_without_a_payload_both_forms_equal_the_restoring_without_a_trim():
    table = _table()
    x = _cols(table, "x", 12)[0]
    u = np.r_[0.0, 0.0, 0.0, 0.0, 0.0]
    outs = [vc.otter(current="full_rotation_rate", hydrostatics=(
        vc.plugins().SurfaceRestoring, {**vc.pick(vc.plugins().SurfaceRestoring, vc.OTTER), "hull_count": 2,
                                        "payload_trim_form": form}),
        vehicle_params={"payload_mass": 0.0, "payload_position": [0.05, 0.0, -0.35]}).signals(vc.mss_to_ours(x), u[:3])
        for form in FORMS]
    for out in outs:
        assert np.array_equal(np.asarray(out["payload_trim_offset"]).ravel(), np.zeros(6))
    assert np.array_equal(np.asarray(outs[0]["hydrostatics.restoring_force"]),
                          np.asarray(outs[1]["hydrostatics.restoring_force"]))


def _weight_balance(form, longitudinal_center_of_flotation):
    """Rows of ``G eta_0`` at the case of the J8 note (25 kg at [0.05, 0, -0.35]), still water, and ``g_0``'s moment."""
    p = vc.plugins()
    vehicle = vc.otter(hydrostatics=(p.SurfaceRestoring, {
        **vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2, "payload_trim_form": form,
        "longitudinal_center_of_flotation": longitudinal_center_of_flotation}))
    signals = vehicle.signals(np.zeros(12), [])
    eta_0 = np.asarray(signals["payload_trim_offset"]).ravel()
    return vc.matrix(signals["restoring_matrix"]) @ eta_0, eta_0, 25.0 * 9.81


def test_G2_the_payload_weight_is_a_heave_wrench_in_the_physics_and_the_mss_form_and_absent_in_the_legacy_form():
    """``G eta_0`` in heave is ``g_0(3)``: the payload weight ``m_P g`` for ``weight_wrench`` and the MSS form (``z = 0``
    is the unloaded waterline, the weight is balanced by the waterplane stiffness), zero for the legacy form."""
    for form in ("weight_wrench", "lever_arm_and_heave_shift"):
        balance, _, weight = _weight_balance(form, -0.2)
        assert balance[2] == pytest.approx(weight, rel=1e-12)
    balance, _, weight = _weight_balance("lever_arm_moment", -0.2)
    assert abs(balance[2]) <= 1e-9 * weight


def test_G2_the_physics_form_and_the_mss_form_have_the_same_shift():
    """Both forms solve ``G eta_0 = g_0`` with the same ``g_0``: the heave shift is the payload weight over the waterplane
    stiffness ``m_P g / (rho g A_wp)`` (``otter_legacy.m`` 187, 190) when ``G`` is uncoupled, and the shifts are equal
    with the Otter's coupled ``G`` too; only the kinematics differ."""
    for lcf in (0.0, -0.2):
        _, physics, weight = _weight_balance("weight_wrench", lcf)
        _, mss, _ = _weight_balance("lever_arm_and_heave_shift", lcf)
        assert np.array_equal(physics, mss)
    _, physics, weight = _weight_balance("weight_wrench", 0.0)
    assert physics[2] == pytest.approx(weight / (1025.0 * 9.81 * 2 * 0.75 * 2.0 * 0.25), rel=1e-12)


def test_G2_the_legacy_form_differs_from_the_physics_form_by_the_heave_term():
    """With the centre of flotation aft of the origin heave and pitch are coupled in ``G``, so dropping the heave term
    also moves the pitch shift (by the weight's moment about the origin instead of the centre of flotation)."""
    _, physics, _ = _weight_balance("weight_wrench", -0.2)
    _, legacy, _ = _weight_balance("lever_arm_moment", -0.2)
    assert abs(physics[4] - legacy[4]) > 1e-3
    assert abs(legacy[2]) < abs(physics[2])


def test_the_kinematics_are_taken_at_the_shifted_attitude_in_the_mss_form_only():
    """``kinematic_attitude_offset`` is ``eta_0`` for the MSS form (``otter_legacy.m`` 266, 269) and zero for the other
    forms."""
    p = vc.plugins()
    for form, expect_shifted in (("lever_arm_and_heave_shift", True), ("weight_wrench", False),
                                 ("lever_arm_moment", False)):
        vehicle = vc.otter(hydrostatics=(p.SurfaceRestoring, {
            **vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2, "payload_trim_form": form}))
        signals = vehicle.signals(np.zeros(12), [])
        offset = np.asarray(signals["kinematic_attitude_offset"]).ravel()
        eta_0 = np.asarray(signals["payload_trim_offset"]).ravel()
        assert (np.array_equal(offset, eta_0) and offset.any()) if expect_shifted else not offset.any()


def test_G4_a_moved_payload_moves_the_trim():
    table = _table()
    x, u = vc.mss_to_ours(_cols(table, "x", 12)[0]), np.zeros(3)
    base = _otter("weight_wrench", 25.0, [0.05, 0.0, -0.35]).signals(x, np.r_[u, 0.0, 0.0])["payload_trim_offset"]
    moved = _otter("weight_wrench", 25.0, [0.35, 0.0, -0.35]).signals(x, np.r_[u, 0.0, 0.0])["payload_trim_offset"]
    assert np.abs(np.asarray(base) - np.asarray(moved)).max() > 1e-3


def test_an_unknown_payload_trim_form_is_refused_naming_the_choices():
    p = vc.plugins()
    with vc.refusing(ValueError, match="payload_trim_form must be one of"):
        vc.otter(hydrostatics=(p.SurfaceRestoring, {**vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2,
                                                    "payload_trim_form": "heave_only"}))


def test_a_monohull_with_a_payload_carries_the_physics_trim_and_reports_it():
    """Single hull, the default form: the shift is in the vehicle's signals (the 2026 trim check of the monohull waits
    on a measured payload and stays skipped in ``test_monohull_grethe.py``)."""
    import test_monohull_grethe as grethe
    vessel = grethe._grethe(vehicle_params={"payload_mass": 25.0, "payload_position": [0.5, 0.2, -0.2]})
    out = vessel.signals(np.zeros(12), [])
    offset = np.asarray(out["payload_trim_offset"]).ravel()
    assert offset[3:5].any() and not np.asarray(out["kinematic_attitude_offset"]).any()
    assert (vc.matrix(out["restoring_matrix"]) @ offset)[2] == pytest.approx(25.0 * 9.81, rel=1e-12)  # the weight, in heave


# ---------------------------------------------------------------------------------------------------------------
# The commanded Otter runs: the physics form against the pinned revision
# ---------------------------------------------------------------------------------------------------------------
OTTER_DATA = Path(__file__).resolve().parents[1] / "vehicle_models" / "data" / "otter"
OTTER_RUNS = ("astern_from_speed", "commands_beyond_the_limits", "differential_turn", "spin_in_place",
              "straight_from_rest_thrust_step", "zigzag_differential_steps")
STEP = 0.05


def _run_table(name, tag):
    data = np.genfromtxt(OTTER_DATA / f"otter_trajectory_{name}_mss_{tag}.csv", delimiter=",", names=True)
    return {key: np.atleast_1d(data[key]) for key in data.dtype.names}


def _tolerances(tag):
    import csv
    table = {}
    with open(OTTER_DATA / f"otter_trajectory_tolerances_mss_{tag}.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            table.setdefault(row["scenario"], np.zeros(12))[int(row["state"]) - 1] = float(row["tolerance"])
    return table


def _otter_with_thrusters(form):
    p = vc.plugins()
    thrusters = [(p.Propulsor, vc.pick(p.Propulsor, vc.reference_part(name)))
                 for name in ("Otter port thruster", "Otter starboard thruster")]
    hydrostatics = (p.SurfaceRestoring, {**vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2,
                                         "payload_trim_form": form})
    return vc.otter(hydrostatics=hydrostatics, actuators=thrusters,
                    vehicle_params={"payload_mass": 25.0, "payload_position": [0.05, 0.0, -0.35]})


def _rk4_run(vehicle, table):
    """The frozen run stepped by the vehicle's own RK4 at the frozen step (MSS state order), commands held per row."""
    inputs = np.column_stack([table["ui_01"], table["ui_02"]])
    steps = int(round((table["t"][1] - table["t"][0]) / STEP))
    x = vc.mss_to_ours(np.array([table[f"x_{i:02d}"][0] for i in range(1, 13)]))
    states = np.zeros((len(table["t"]), 12))
    for row in range(len(table["t"])):
        states[row] = vc.ours_to_mss(x)
        if row < len(table["t"]) - 1:
            for _ in range(steps):
                x = vehicle.step_vector(x, inputs[row], STEP)
    return states


def _over(form, tag, name):
    table = _run_table(name, tag)
    reference = np.column_stack([table[f"x_{i:02d}"] for i in range(1, 13)])
    return np.abs(_rk4_run(_otter_with_thrusters(form), table) - reference).max(axis=0) / _tolerances(tag)[name]


@pytest.mark.parametrize("name", OTTER_RUNS)
def test_G3_the_physics_form_reproduces_the_commanded_otter_run_of_the_pinned_revision(name):
    """Within the frozen tolerance of every state (the run of the unmodified function at the pinned revision, kinematics
    at the actual attitude)."""
    assert (_over("weight_wrench", "aff76cf", name) <= 1.0).all()


@pytest.mark.parametrize("name", ["differential_turn", "zigzag_differential_steps"])
def test_G4_the_shifted_attitude_form_misses_the_pinned_revision(name):
    """The gate can fail: the form with the kinematics at the shifted attitude leaves the run of the pinned revision by
    far more than the tolerance."""
    assert (_over("lever_arm_and_heave_shift", "aff76cf", name) > 10.0).any()
