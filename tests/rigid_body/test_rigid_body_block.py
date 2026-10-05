"""Gate tests for the rigid-body + added-mass block (Otter-based catamaran).

Written before the block existed (job A-4a, 2026-10-05); the block was ported
by job A-4b. Repaired by job A-4c (2026-10-05) for owner decision E-11 c: the
block has a ``legacy_otter_inertia`` flag, default False = the corrected MSS
inertia. The legacy path (flag True) is tested against the frozen reference and
the numpy source (G1, G2); the default path against MSS ``otter.m`` (G5); G4 runs
on both paths.

Contract the porter must provide
--------------------------------
``more_dynamics.models.rigid_body`` exports:

* ``preprocess_rigid_body(**PARAMETERS, legacy_otter_inertia=False)
  -> RigidBodyConstants`` (numpy). The keyword names are the keys of
  ``PARAMETERS`` below plus the flag (E-11 c).
* ``RigidBodyConstants``: frozen dataclass with at least
  ``rigid_body_mass_matrix`` (M_RB), ``added_mass_matrix`` (M_A) and
  ``total_mass_matrix`` (M = M_RB + M_A), each a 6x6 ``np.ndarray``, and
  ``inertia`` (3x3) and ``center_of_gravity`` (3,), checked on the default path.
* ``rigid_body_casadi(constants) -> ca.Function`` with inputs named
  ``["nu", "nu_r"]`` (6x1 each) and outputs named ``["M", "C_RB", "C_A"]``
  (6x6 each): ``C_RB(nu)`` and ``C_A(nu_r)``.

Map rows covered (``agents-more/MIGRATION_MAP.md``), numpy source read in full
(paths relative to ``more_generic_models/more_generic_models/``):

* ``dynamics/plant/asv_catamaran/asv_catamaran.py``: ``compute_constant_values``
  121-153 (M_RB, M_A, M), ``get_M_RB`` 203, ``get_M_A`` 206,
  ``get_C_RB`` 221-222, ``get_C_A`` 224-225.
* ``dynamics/plant/matrices/rigid_body_kinetics_surface_vessel.py``:
  ``get_params_catamaran`` 17-51, ``M_RB`` 96-121, ``C_RB_runtime`` 147-170.
* ``dynamics/plant/matrices/added_mass.py``: ``M_A_6dof`` 24-41,
  ``C_A_6dof_lagrangian`` 116-122, ``get_added_mass_derivates_catamaran``
  269-289, ``added_mass_surge`` 313-335, ``get_C_RB_lagrangian`` 362-388.
* ``config/dataclass/plant/asv_catamaran_params.py`` 8-31 (parameter values).
* ``rigid_body_kinetics.py::RigidBody6DOF`` is not used by the catamaran
  (``added_mass.py`` imports it, line 3, but never calls it).

Conventions found in the source (hidden assumptions, see the ledger)
---------------------------------------------------------------------
* Body frame, z down; CO is the body origin; ``r_g`` is CO -> CG.
* ``H(r) = [[I, S(r)^T], [0, I]]``; ``M_RB = H^T diag(m I, I_o) H`` with
  ``m = m_hull + m_payload`` and ``r_g = (m_hull r_hull + m_p r_p) / m``.
* Legacy path: ``I_o = I_CG - m_hull S(r_g)^2 - m_p S(r_p)^2`` is already
  about the CO and is shifted a second time by ``H`` (as in the reference
  generator). Default path: the inertia about the combined CG, as MSS
  ``otter.m`` lines 123-129 since its revision of 2026-04-20 (line 78).
* ``C_RB = H^T diag(m S(w), -S(I_o w)) H`` depends only on ``w = nu[3:]``.
* ``M_A = -diag(c * [A11, m_hull, m_hull, I_o[0,0], I_o[1,1], I_o[2,2]])``,
  ``A11 = 2.7 rho (m_hull/rho)^(5/3) / L^2``; added mass uses the hull mass
  only, and the rotational terms use the CO inertia.
* ``C_A = m2c(M_A, nu_r)`` (Fossen 2021 Theorem 3.2 form), no Munk-moment
  cancellation.

Frozen reference: ``tests/data/rigid_body/`` (``SOURCE.md`` names origin,
revisions and layout).
"""

import importlib
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "rigid_body"
CONTRACT_MODULE = "more_dynamics.models.rigid_body"

G1_TOLERANCE = 1e-9   # 30_checks/README.md, gate G1
G2_TOLERANCE = 1e-10  # 30_checks/README.md, gate G2
G4_FACTOR = 10.0      # 30_checks/README.md, gate G4
SEED = 20261005
N_RANDOM_STATES = 1000

# Values from config/dataclass/plant/asv_catamaran_params.py (line per key);
# they equal test_dynamics_consistency.m lines 12-13, 138-140 and 165-167.
PARAMETERS = {
    "length": 2.0,                                   # l, line 8
    "beam": 1.08,                                    # b, line 9
    "water_density": 1025.0,                         # rho, line 10
    "hull_mass": 55.0,                               # m_hull, line 12
    "payload_mass": 25.0,                            # m_payload, line 13
    "hull_center_of_gravity": [0.2, 0.0, -0.2],      # r_hull, line 14
    "payload_position": [0.05, 0.0, -0.35],          # r_payload, line 15
    "added_mass_coefficients": [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7],  # 21-23
    "radii_of_gyration": [0.4, 0.25, 0.25],          # R_456_scale, line 28
}

# Current of the reference generator, test_dynamics_consistency.m lines 14-15.
CURRENT_SPEED = 0.3
CURRENT_DIRECTION = np.deg2rad(30.0)

# Owner decision E-11 c (agents-more/80_owner/inbox/E-11_otter_inertia.md): the
# frozen reference and the numpy source are the legacy inertia path.
LEGACY = {"legacy_otter_inertia": True}
INERTIA_PATHS = {"legacy": True, "default": False}


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _contract():
    try:
        return importlib.import_module(CONTRACT_MODULE)
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed (owner decision E-7): {exc}")
        pytest.fail(f"block not ported yet: {exc}")


def _build(**overrides):
    block = _contract()
    params = {**PARAMETERS, **overrides}
    constants = block.preprocess_rigid_body(**params)
    function = block.rigid_body_casadi(constants)
    return constants, function


def _evaluate(function, nu, nu_r):
    out = function(nu=np.asarray(nu, float), nu_r=np.asarray(nu_r, float))
    return {name: np.array(out[name], dtype=float) for name in ("M", "C_RB", "C_A")}


def _load_csv(name):
    path = DATA_DIR / name
    header = path.read_text().splitlines()[0].split(",")
    values = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return header, values


def _reference_cases():
    """The 50 frozen cases: nu, nu_r and the reference matrices."""
    in_header, inputs = _load_csv("inputs.csv")
    ref_header, ref = _load_csv("matlab_reference.csv")
    assert in_header[:12] == [f"x{i}" for i in range(1, 13)]

    def matrix(prefix):
        cols = [ref_header.index(f"{prefix}_{i:02d}") for i in range(1, 37)]
        return ref[:, cols].reshape(-1, 6, 6)  # row-major, SOURCE.md

    nu = inputs[:, 0:6]
    psi = inputs[:, 11]
    nu_c = np.zeros_like(nu)
    nu_c[:, 0] = CURRENT_SPEED * np.cos(CURRENT_DIRECTION - psi)
    nu_c[:, 1] = CURRENT_SPEED * np.sin(CURRENT_DIRECTION - psi)
    return {
        "nu": nu,
        "nu_r": nu - nu_c,
        "M_RB": matrix("M_RB"),
        "M_A": matrix("M_A"),
        "C_RB": matrix("C_RB"),
        "C_A": matrix("C_A"),
        "C_total": matrix("C_total"),
    }


def _random_states():
    rng = np.random.default_rng(SEED)
    nu = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    nu_r = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    return nu, nu_r


def _numpy_source():
    try:
        from more_generic_models.dynamics.plant.asv_catamaran.asv_catamaran import (
            ASVCatamaran,
        )
    except ImportError as exc:
        pytest.skip(
            "numpy source not importable; install it into the venv with "
            "pip install -e <more>/more_generic_models (job A-4c)"
            f" — {exc}"
        )
    return ASVCatamaran()


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


# --------------------------------------------------------------------------
# Reference sanity (no port needed)
# --------------------------------------------------------------------------
def test_reference_shapes_and_current_model():
    ref = _reference_cases()
    assert ref["nu"].shape == (50, 6)
    for key in ("M_RB", "M_A", "C_RB", "C_A", "C_total"):
        assert ref[key].shape == (50, 6, 6), key
    np.testing.assert_allclose(
        ref["C_total"], ref["C_RB"] + ref["C_A"], atol=G1_TOLERANCE, rtol=0.0
    )


# --------------------------------------------------------------------------
# G1 — numpy source vs MATLAB reference
# --------------------------------------------------------------------------
def test_G1_numpy_source_matches_matlab_reference():
    vessel = _numpy_source()
    ref = _reference_cases()
    for k in range(len(ref["nu"])):
        nu, nu_r = ref["nu"][k], ref["nu_r"][k]
        checks = {
            "M_RB": (vessel.get_M_RB(), ref["M_RB"][k]),
            "M_A": (vessel.get_M_A(), ref["M_A"][k]),
            "C_RB": (vessel.get_C_RB(nu), ref["C_RB"][k]),
            "C_A": (vessel.get_C_A(nu_r), ref["C_A"][k]),
            "C_total": (vessel.get_C_total(nu, nu_r), ref["C_total"][k]),
        }
        for name, (value, expected) in checks.items():
            assert _max_diff(value, expected) <= G1_TOLERANCE, (name, k)


def test_G1_block_matches_matlab_reference():
    constants, function = _build(**LEGACY)
    ref = _reference_cases()
    np.testing.assert_allclose(
        constants.rigid_body_mass_matrix, ref["M_RB"][0], atol=G1_TOLERANCE, rtol=0.0
    )
    np.testing.assert_allclose(
        constants.added_mass_matrix, ref["M_A"][0], atol=G1_TOLERANCE, rtol=0.0
    )
    for k in range(len(ref["nu"])):
        out = _evaluate(function, ref["nu"][k], ref["nu_r"][k])
        expected = {
            "M": ref["M_RB"][k] + ref["M_A"][k],
            "C_RB": ref["C_RB"][k],
            "C_A": ref["C_A"][k],
        }
        for name, value in expected.items():
            assert _max_diff(out[name], value) <= G1_TOLERANCE, (name, k)
        assert _max_diff(out["C_RB"] + out["C_A"], ref["C_total"][k]) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G2 — CasADi block vs numpy source
# --------------------------------------------------------------------------
def test_G2_function_signature():
    constants, function = _build()
    assert function.name_in() == ["nu", "nu_r"]
    assert function.name_out() == ["M", "C_RB", "C_A"]
    for i in range(2):
        assert function.size_in(i) == (6, 1)
    for i in range(3):
        assert function.size_out(i) == (6, 6)
    np.testing.assert_allclose(
        constants.total_mass_matrix,
        constants.rigid_body_mass_matrix + constants.added_mass_matrix,
        atol=G2_TOLERANCE,
        rtol=0.0,
    )


def _assert_block_equals_source(function, constants, vessel, nu_cases, nu_r_cases):
    np.testing.assert_allclose(
        constants.rigid_body_mass_matrix, vessel.get_M_RB(), atol=G2_TOLERANCE, rtol=0.0
    )
    np.testing.assert_allclose(
        constants.added_mass_matrix, vessel.get_M_A(), atol=G2_TOLERANCE, rtol=0.0
    )
    for k, (nu, nu_r) in enumerate(zip(nu_cases, nu_r_cases)):
        out = _evaluate(function, nu, nu_r)
        assert _max_diff(out["M"], vessel.get_M_total()) <= G2_TOLERANCE, ("M", k)
        assert _max_diff(out["C_RB"], vessel.get_C_RB(nu)) <= G2_TOLERANCE, ("C_RB", k)
        assert _max_diff(out["C_A"], vessel.get_C_A(nu_r)) <= G2_TOLERANCE, ("C_A", k)


def test_G2_block_matches_numpy_source_on_reference_cases():
    constants, function = _build(**LEGACY)
    vessel = _numpy_source()
    ref = _reference_cases()
    _assert_block_equals_source(function, constants, vessel, ref["nu"], ref["nu_r"])


def test_G2_block_matches_numpy_source_on_random_states():
    constants, function = _build(**LEGACY)
    vessel = _numpy_source()
    nu, nu_r = _random_states()
    _assert_block_equals_source(function, constants, vessel, nu, nu_r)


# --------------------------------------------------------------------------
# G4 — a broken model must fail
# --------------------------------------------------------------------------
def _perturbations():
    hull_cg = np.array(PARAMETERS["hull_center_of_gravity"], dtype=float)
    added = np.array(PARAMETERS["added_mass_coefficients"], dtype=float)
    added_sway_zero = added.copy()
    added_sway_zero[1] = 0.0
    return {
        "hull_mass_plus_1_percent": {"hull_mass": PARAMETERS["hull_mass"] * 1.01},
        "hull_cg_x_plus_1_cm": {"hull_center_of_gravity": hull_cg + [0.01, 0.0, 0.0]},
        "added_mass_sway_times_zero": {"added_mass_coefficients": added_sway_zero},
    }


def _path_reference(legacy):
    """Per reference case: M, C_RB and C_A the unperturbed block must reproduce.

    Legacy path: the frozen MATLAB reference. Default path: M and C_RB re-derived
    from MSS otter.m (``_otter_m_current_matrices``); otter.m has no frozen C_A,
    so C_A is not compared on that path (None).
    """
    ref = _reference_cases()
    if legacy:
        return [
            {"M": ref["M_RB"][k] + ref["M_A"][k], "C_RB": ref["C_RB"][k], "C_A": ref["C_A"][k]}
            for k in range(len(ref["nu"]))
        ]
    mrb, ma, ig, rg_total = _otter_m_current_matrices()
    return [
        {"M": mrb + ma, "C_RB": _otter_m_crb(ig, rg_total, ref["nu"][k][3:]), "C_A": None}
        for k in range(len(ref["nu"]))
    ]


def _worst_diff(function, legacy):
    ref = _reference_cases()
    worst = 0.0
    for k, expected in enumerate(_path_reference(legacy)):
        out = _evaluate(function, ref["nu"][k], ref["nu_r"][k])
        for name, value in expected.items():
            if value is not None:
                worst = max(worst, _max_diff(out[name], value))
    return worst


@pytest.mark.parametrize("path", sorted(INERTIA_PATHS))
@pytest.mark.parametrize("name", sorted(_perturbations()))
def test_G4_perturbed_block_is_detected(name, path):
    legacy = INERTIA_PATHS[path]
    flag = {"legacy_otter_inertia": legacy}
    # Control: the unperturbed block matches this path's reference, so a
    # detection below is caused by the perturbation, not by the path.
    _, control = _build(**flag)
    assert _worst_diff(control, legacy) <= G1_TOLERANCE, (path, "control")
    _, function = _build(**flag, **_perturbations()[name])
    worst = _worst_diff(function, legacy)
    assert worst > G4_FACTOR * G1_TOLERANCE, (name, path, worst)


# --------------------------------------------------------------------------
# G5 — Otter constants, MSS otter.m
# (source-sim/MSS/CRAFT/USV/models/otter.m, MSS 99bf0b3, file e1dff2a)
# --------------------------------------------------------------------------
OTTER_M = Path(__file__).resolve().parents[4] / "source-sim/MSS/CRAFT/USV/models/otter.m"
OTTER_LINES = {  # line number -> text the constants below were taken from
    91: "rho = 1025;",
    92: "L = 2.0;",
    93: "B = 1.08;",
    94: "m = 55.0;",
    95: "rg = [0.2 0 -0.2]';",
    96: "R44 = 0.4 * B;",
    97: "R55 = 0.25 * L;",
    98: "R66 = 0.25 * L;",
    123: "Ig_CG = m * diag([R44^2, R55^2, R66^2]);",
    124: "rg_hull = rg;",
    125: "rg_total = (m*rg_hull + mp*rp)/(m+mp);",
    126: "r_hull = rg_hull - rg_total;",
    127: "r_payload = rp - rg_total;",
    128: "Ig = Ig_CG - m * Smtrx(r_hull)^2 - mp * Smtrx(r_payload)^2;",
    129: "rg = rg_total;",
    143: "MRB_CG = [ (m+mp) * I3  O3",
    145: "CRB_CG = [ (m+mp) * Smtrx(nu2)         O3",
    146: "O3               -Smtrx(Ig*nu2)  ];",
    148: "H = Hmtrx(rg);",
    153: "Xudot = -addedMassSurge(m,L,rho);",
    154: "Yvdot = -1.5 * m;",
    155: "Zwdot = -1.0 * m;",
    156: "Kpdot = -0.2 * Ig(1,1);",
    157: "Mqdot = -0.8 * Ig(2,2);",
    158: "Nrdot = -1.7 * Ig(3,3);",
}
OTTER = {"rho": 1025.0, "L": 2.0, "B": 1.08, "m": 55.0, "rg": [0.2, 0.0, -0.2],
         "R": [0.4, 0.25, 0.25], "added": [1.5, 1.0, 0.2, 0.8, 1.7]}
# Payload of the reference generator (test_dynamics_consistency.m lines 12-13);
# otter.m takes mp and rp as arguments.
OTTER_MP, OTTER_RP = 25.0, [0.05, 0.0, -0.35]


def test_G5_cited_otter_lines_are_unchanged():
    if not OTTER_M.exists():
        pytest.skip(f"MSS checkout not found at {OTTER_M}")
    lines = OTTER_M.read_text().splitlines()
    for number, text in OTTER_LINES.items():
        assert lines[number - 1].strip().startswith(text), (number, lines[number - 1])


def _otter_inputs_as_parameters():
    return {
        **PARAMETERS,
        "length": OTTER["L"],
        "beam": OTTER["B"],
        "water_density": OTTER["rho"],
        "hull_mass": OTTER["m"],
        "payload_mass": OTTER_MP,
        "hull_center_of_gravity": OTTER["rg"],
        "payload_position": OTTER_RP,
        "radii_of_gyration": OTTER["R"],
        "added_mass_coefficients": [-1.0] + [-c for c in OTTER["added"]],
    }


def test_G5_otter_translational_constants():
    constants, _ = _build(**_otter_inputs_as_parameters())
    m, mp, L, rho = OTTER["m"], OTTER_MP, OTTER["L"], OTTER["rho"]
    # otter.m 143: MRB_CG(1:3,1:3) = (m+mp)*I3; H keeps this block (Hmtrx.m 17-18)
    np.testing.assert_allclose(
        constants.rigid_body_mass_matrix[:3, :3], (m + mp) * np.eye(3),
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m line 143",
    )
    # otter.m 153 + addedMassSurge.m 33-34: A11 = 2.7*rho*nabla^(5/3)/L^2, nabla = m/rho
    a11 = 2.7 * rho * (m / rho) ** (5.0 / 3.0) / L**2
    # otter.m 154-155, 160: MA = -diag(Xudot, Yvdot, Zwdot, ...)
    expected = [a11, OTTER["added"][0] * m, OTTER["added"][1] * m]
    np.testing.assert_allclose(
        np.diag(constants.added_mass_matrix)[:3], expected,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 153-155, 160",
    )


def _skew(v):
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def _otter_m_current_inertia():
    """otter.m lines 123-129 (revision 2026-04-20): inertia about the combined CG."""
    m, mp = OTTER["m"], OTTER_MP
    rg, rp = np.array(OTTER["rg"]), np.array(OTTER_RP)
    radii = np.array(OTTER["R"]) * np.array([OTTER["B"], OTTER["L"], OTTER["L"]])
    ig_cg = m * np.diag(radii**2)
    rg_total = (m * rg + mp * rp) / (m + mp)
    r_hull, r_payload = rg - rg_total, rp - rg_total
    ig = ig_cg - m * _skew(r_hull) @ _skew(r_hull) - mp * _skew(r_payload) @ _skew(r_payload)
    return ig, rg_total


def test_G5_otter_rotational_inertia_current_mss():
    # Default path (E-11 c). Until A-4c this was a strict expected failure.
    constants, _ = _build(**_otter_inputs_as_parameters())
    assert constants.legacy_otter_inertia is False
    ig, rg_total = _otter_m_current_inertia()
    m_total = OTTER["m"] + OTTER_MP
    h = np.block([[np.eye(3), _skew(rg_total).T], [np.zeros((3, 3)), np.eye(3)]])
    mrb_cg = np.block([[m_total * np.eye(3), np.zeros((3, 3))], [np.zeros((3, 3)), ig]])
    # otter.m lines 143-144, 148-149
    np.testing.assert_allclose(
        constants.rigid_body_mass_matrix, h.T @ mrb_cg @ h,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 123-129, 143-149",
    )
    # otter.m lines 156-158, 160
    expected = [c * ig[i, i] for i, c in enumerate(OTTER["added"][2:])]
    np.testing.assert_allclose(
        np.diag(constants.added_mass_matrix)[3:], expected,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 156-158, 160",
    )


def _otter_m_current_matrices():
    """MRB and MA of otter.m (current), on the frozen parameters.

    MRB: lines 143-144 (MRB_CG), 148-149 (H = Hmtrx(rg), MRB = H' MRB_CG H) with
    rg = rg_total (line 129). MA: lines 153-158, 160; Xudot from line 153 and
    addedMassSurge.m lines 33-34 (as in ``test_G5_otter_translational_constants``).
    """
    m, mp, L, rho = OTTER["m"], OTTER_MP, OTTER["L"], OTTER["rho"]
    ig, rg_total = _otter_m_current_inertia()
    h = np.block([[np.eye(3), _skew(rg_total).T], [np.zeros((3, 3)), np.eye(3)]])
    mrb_cg = np.block([[(m + mp) * np.eye(3), np.zeros((3, 3))], [np.zeros((3, 3)), ig]])
    a11 = 2.7 * rho * (m / rho) ** (5.0 / 3.0) / L**2
    y_v, z_w, k_p, m_q, n_r = OTTER["added"]
    ma = np.diag([a11, y_v * m, z_w * m, k_p * ig[0, 0], m_q * ig[1, 1], n_r * ig[2, 2]])
    return h.T @ mrb_cg @ h, ma, ig, rg_total


def _otter_m_crb(ig, rg_total, nu2):
    """otter.m lines 145-146 (CRB_CG), 148, 150 (CRB = H' CRB_CG H)."""
    m_total = OTTER["m"] + OTTER_MP
    h = np.block([[np.eye(3), _skew(rg_total).T], [np.zeros((3, 3)), np.eye(3)]])
    crb_cg = np.block([
        [m_total * _skew(nu2), np.zeros((3, 3))],
        [np.zeros((3, 3)), -_skew(ig @ nu2)],
    ])
    return h.T @ crb_cg @ h


def test_G5_default_path_inertia_cg_and_rotational_added_mass():
    """Default path (E-11 c) against MSS otter.m, computed here from the lines.

    Every constant comes from a line pinned by ``test_G5_cited_otter_lines_are_unchanged``
    (91-98 data, 123-129 inertia, 156-158 rotational added mass); the payload
    from the reference generator (``OTTER_MP``, ``OTTER_RP``).
    """
    constants, _ = _build(**_otter_inputs_as_parameters())
    assert constants.legacy_otter_inertia is False
    m, mp = OTTER["m"], OTTER_MP
    rg_hull, rp = np.array(OTTER["rg"]), np.array(OTTER_RP)  # lines 95, 124; argument rp
    # line 96-98: R44 = 0.4*B, R55 = 0.25*L, R66 = 0.25*L
    r44, r55, r66 = OTTER["R"][0] * OTTER["B"], OTTER["R"][1] * OTTER["L"], OTTER["R"][2] * OTTER["L"]
    ig_cg = m * np.diag([r44**2, r55**2, r66**2])                 # line 123
    rg_total = (m * rg_hull + mp * rp) / (m + mp)                  # line 125
    r_hull = rg_hull - rg_total                                    # line 126
    r_payload = rp - rg_total                                      # line 127
    s_hull, s_payload = _skew(r_hull), _skew(r_payload)
    ig = ig_cg - m * s_hull @ s_hull - mp * s_payload @ s_payload  # line 128 (Smtrx(.)^2)
    np.testing.assert_allclose(
        constants.inertia, ig, atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 123-128"
    )
    np.testing.assert_allclose(  # line 129: rg = rg_total, the CG used by H (line 148)
        constants.center_of_gravity, rg_total,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 125, 129",
    )
    # lines 156-158: Kpdot = -0.2 Ig(1,1), Mqdot = -0.8 Ig(2,2), Nrdot = -1.7 Ig(3,3);
    # line 160: MA = -diag([..., Kpdot, Mqdot, Nrdot])
    k_pdot = -OTTER["added"][2] * ig[0, 0]
    m_qdot = -OTTER["added"][3] * ig[1, 1]
    n_rdot = -OTTER["added"][4] * ig[2, 2]
    np.testing.assert_allclose(
        np.diag(constants.added_mass_matrix)[3:], [-k_pdot, -m_qdot, -n_rdot],
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 156-158, 160",
    )
    # the default path must differ from the legacy one (otherwise the flag is dead)
    legacy, _ = _build(**_otter_inputs_as_parameters(), **LEGACY)
    assert _max_diff(legacy.inertia, constants.inertia) > G4_FACTOR * G1_TOLERANCE
