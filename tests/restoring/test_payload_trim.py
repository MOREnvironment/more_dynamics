"""Gate tests for the payload trim of the surface restoring (``payload_trim_form``).

A point payload off the centre of flotation shifts the equilibrium of the linear restoring, ``g = G (eta - eta_0)``
(``models/restoring/payload_trim.py``). Two forms: ``lever_arm_moment`` (the default: the payload weight is in the
loaded displacement, only the moment of its lever arm shifts the equilibrium) and ``lever_arm_and_heave_shift`` (MSS
``otter.m`` 259-269, a comparison form that counts the weight a second time in heave and takes the kinematics at the
shifted attitude).

Gates (the test names carry them): G1 the vehicle against MATLAB running MSS (frozen CSV
``data/payload_trim_mss_ed2ce7b.csv``, ``data/generate_payload_trim_mss.m``): the shift, the restoring vector and the
whole derivative of the Otter with two thrusters, for the MSS form against the unmodified ``otter.m`` and for the
corrected form against the same call with its two terms changed (written out in the generator), 1e-9; G2 properties the
equations imply (the heave balance, the lever-arm moment shared by both forms, the loaded displacement); G4 a
perturbed payload or the other form is detected.

Hidden assumptions: the corrected form is a second implementation written from ``otter.m``'s own workspace, not a
published result (its physics is the argument of ``models/restoring/payload_trim.py``); the reference is MSS, not the
boat; the shifted-attitude kinematics of the MSS form need the vehicle's ``kinematic_attitude_offset``.

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
FORMS = ("lever_arm_and_heave_shift", "lever_arm_moment")
SUFFIX = {"lever_arm_and_heave_shift": "mss", "lever_arm_moment": "corr"}


def _table():
    data = np.genfromtxt(DATA / "payload_trim_mss_ed2ce7b.csv", delimiter=",", names=True)
    return {name: np.atleast_1d(data[name]) for name in data.dtype.names}


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
                 "lever_arm_moment": ("eta0c", "g_corr", "xdot_corr")}[form]
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
    worst = _evaluate(_table(), form)
    assert max(worst.values()) <= G1_TOLERANCE, worst


def test_G4_each_form_is_told_from_the_other_by_the_reference():
    """A perturbed model is detected: the shift of the other form misses the reference by more than the gate."""
    table = _table()
    mss = _evaluate(table, "lever_arm_and_heave_shift")
    corrected_against_mss = {"eta_0": 0.0}
    vehicles = {}
    for k, mp, rp in _cases(table):
        if mp == 0.0:
            continue
        key = (mp, tuple(rp))
        vehicles.setdefault(key, _otter("lever_arm_moment", mp, rp))
        x = vc.mss_to_ours(_cols(table, "x", 12)[k])
        u = np.r_[table["Vc"][k], table["beta_c"][k], 0.0, _cols(table, "n", 2)[k]]
        eta_0 = np.asarray(vehicles[key].signals(x, u)["payload_trim_offset"]).ravel()
        corrected_against_mss["eta_0"] = max(corrected_against_mss["eta_0"],
                                             np.abs(eta_0 - _cols(table, "eta0", 6)[k]).max())
    assert max(mss.values()) <= G1_TOLERANCE
    assert corrected_against_mss["eta_0"] > 1e-3, corrected_against_mss  # the heave shift is the difference


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


def test_G2_the_corrected_form_leaves_no_unbalanced_heave_force_and_the_mss_form_leaves_the_payload_weight():
    """``G eta_0`` in heave is ``g_0(3)``: zero for ``lever_arm_moment`` (buoyancy already carries the weight, the loaded
    displacement), the payload weight ``m_P g`` for the MSS form (the double count of F-MSS-31)."""
    balance, _, weight = _weight_balance("lever_arm_moment", -0.2)
    assert abs(balance[2]) <= 1e-9 * weight
    balance, _, weight = _weight_balance("lever_arm_and_heave_shift", -0.2)
    assert balance[2] == pytest.approx(weight, rel=1e-12)


def test_G2_the_two_forms_shift_roll_and_pitch_equally_when_heave_and_pitch_are_uncoupled():
    """With the centre of flotation at the origin ``G`` has no heave-pitch coupling: the moment of the lever arm gives
    the same roll and pitch shift in both forms, only the heave shift differs (it is the payload weight over the
    waterplane stiffness in the MSS form and zero in the other)."""
    _, corrected, weight = _weight_balance("lever_arm_moment", 0.0)
    _, mss, _ = _weight_balance("lever_arm_and_heave_shift", 0.0)
    assert np.allclose(corrected[3:5], mss[3:5], rtol=1e-12, atol=0.0)
    assert corrected[2] == pytest.approx(0.0, abs=1e-15) and mss[2] > 0.0
    assert mss[2] == pytest.approx(weight / (1025.0 * 9.81 * 2 * 0.75 * 2.0 * 0.25), rel=1e-12)  # m_P g / (rho g A_wp), otter.m 187, 190


def test_G2_a_trim_with_a_coupled_restoring_moves_pitch_with_heave():
    """With the centre of flotation aft of the origin (the Otter's -0.2 m) heave and pitch are coupled in ``G``, so the
    MSS form's heave shift also moves its pitch: the two forms' pitch shifts differ, by the heave weight's moment."""
    _, corrected, _ = _weight_balance("lever_arm_moment", -0.2)
    _, mss, _ = _weight_balance("lever_arm_and_heave_shift", -0.2)
    assert abs(mss[4] - corrected[4]) > 1e-3
    assert abs(corrected[2]) < abs(mss[2])


def test_the_kinematics_are_taken_at_the_shifted_attitude_in_the_mss_form_only():
    """``kinematic_attitude_offset`` is ``eta_0`` for the MSS form (``otter.m`` 266, 269) and zero for the other."""
    p = vc.plugins()
    for form, expect_shifted in (("lever_arm_and_heave_shift", True), ("lever_arm_moment", False)):
        vehicle = vc.otter(hydrostatics=(p.SurfaceRestoring, {
            **vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2, "payload_trim_form": form}))
        signals = vehicle.signals(np.zeros(12), [])
        offset = np.asarray(signals["kinematic_attitude_offset"]).ravel()
        eta_0 = np.asarray(signals["payload_trim_offset"]).ravel()
        assert (np.array_equal(offset, eta_0) and offset.any()) if expect_shifted else not offset.any()


def test_G4_a_moved_payload_moves_the_trim():
    table = _table()
    x, u = vc.mss_to_ours(_cols(table, "x", 12)[0]), np.zeros(3)
    base = _otter("lever_arm_moment", 25.0, [0.05, 0.0, -0.35]).signals(x, np.r_[u, 0.0, 0.0])["payload_trim_offset"]
    moved = _otter("lever_arm_moment", 25.0, [0.35, 0.0, -0.35]).signals(x, np.r_[u, 0.0, 0.0])["payload_trim_offset"]
    assert np.abs(np.asarray(base) - np.asarray(moved)).max() > 1e-3


def test_an_unknown_payload_trim_form_is_refused_naming_the_choices():
    p = vc.plugins()
    with pytest.raises(ValueError, match="payload_trim_form must be one of"):
        vc.otter(hydrostatics=(p.SurfaceRestoring, {**vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2,
                                                    "payload_trim_form": "heave_only"}))


def test_a_monohull_with_a_payload_carries_the_corrected_trim_and_reports_it():
    """Single hull, the default form: the shift is in the vehicle's signals (the 2026 trim check of the monohull waits
    on a measured payload and stays skipped in ``test_monohull_grethe.py``)."""
    import test_monohull_grethe as grethe
    vessel = grethe._grethe(vehicle_params={"payload_mass": 25.0, "payload_position": [0.5, 0.2, -0.2]})
    out = vessel.signals(np.zeros(12), [])
    offset = np.asarray(out["payload_trim_offset"]).ravel()
    assert offset[3:5].any() and not np.asarray(out["kinematic_attitude_offset"]).any()
    assert abs((vc.matrix(out["restoring_matrix"]) @ offset)[2]) <= 1e-9  # no payload weight in heave
