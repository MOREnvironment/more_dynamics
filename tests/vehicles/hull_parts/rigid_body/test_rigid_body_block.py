"""Gate tests for the rigid-body + added-mass block (Otter-based catamaran).

The block has no ``legacy_otter_inertia`` flag: it has one path, the inertia
about the combined CG of current ``otter.m``. G1 and G4 read the MATLAB
reference recomputed by current MSS on 2026-10-05
(``matlab_reference_mss_current.csv``); the legacy ``matlab_reference.csv`` stays
on disk and no test reads it. The ``otter.m`` transcriptions (G5) are the
second check of the whole matrices.

Contract of the block
---------------------
``more_dynamics.models.vehicles.hull_parts.rigid_body`` exports (all CasADi, no numpy in the
block):

* ``rigid_body_parameters(mass_properties="hull_with_payload")``: the
  declared parameters; their names are the keys of ``PARAMETERS`` below.
* ``rigid_body_casadi(*, mass_properties, coriolis,
  stabilize_added_mass_coriolis) -> ca.Function`` with inputs ``nu``,
  ``nu_r`` (6x1 each) and one input per declared parameter, by name, and
  outputs by name: ``M``, ``C_RB``, ``C_A``, ``M_RB``, ``M_A`` (6x6 each;
  ``C_RB(nu)``, ``C_A(nu_r)``, ``M = M_RB + M_A``), ``mass``,
  ``center_of_gravity`` (3x1, CO -> CG), ``inertia`` (3x3, about the CG). No
  ``legacy_otter_inertia`` keyword, input or output.
* ``check_rigid_body_values(values, *, mass_properties)``: the numbers
  checked against the declaration (and ``M`` symmetric positive definite)
  where they enter, returned as ``ca.DM``.

Every test calls the block with numbers, ``rigid_body_casadi(...)(nu=...,
nu_r=..., **values)``, and reads the named outputs; ``_build`` returns them
(the state-independent ones as a namespace with the attribute names the
assertions use) together with the block and its checked values.

Conventions (hidden assumptions)
--------------------------------
* Body frame, z down; CO is the body origin; ``r_g`` is CO -> CG.
* ``H(r) = [[I, S(r)^T], [0, I]]``; ``M_RB = H^T diag(m I, I_o) H`` with
  ``m = m_hull + m_payload`` and ``r_g = (m_hull r_hull + m_p r_p) / m``.
* The template generator of ``matlab_reference.csv`` uses ``I_o = I_CG -
  m_hull S(r_g)^2 - m_p S(r_p)^2``, already about the CO, and shifts it a second
  time by ``H``. The block uses the inertia about the combined CG, as MSS ``otter.m``
  lines 122-128 since its revision of 2026-04-20 (line 77).
* ``C_RB = H^T diag(m S(w), -S(I_o w)) H`` depends only on ``w = nu[3:]``.
* ``M_A = -diag(c * [A11, m_hull, m_hull, I_o[0,0], I_o[1,1], I_o[2,2]])``,
  ``A11 = 2.7 rho (m_hull/rho)^(5/3) / L^2``; added mass uses the hull mass
  only, and the rotational terms use the CO inertia.
* ``C_A = m2c(M_A, nu_r)`` (Fossen 2021 Theorem 3.2 form), no Munk-moment
  cancellation.

Frozen reference: ``tests/data/rigid_body/matlab_reference_mss_current.csv``
(``SOURCE.md`` names origin, revisions and layout).
"""

import importlib
import inspect
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[3] / "data" / "rigid_body"
CONTRACT_MODULE = "more_dynamics.models.vehicles.hull_parts.rigid_body"

G1_TOLERANCE = 1e-9   # G1: block vs MATLAB running MSS, absolute
G2_TOLERANCE = 1e-10  # G2: block vs a transcription or its own frozen form, absolute
G4_FACTOR = 10.0      # G4: a perturbed model must differ by more than 10 x G1
SEED = 20261005
N_RANDOM_STATES = 1000

# Values of MSS otter.m (lines at cc07579, pinned in OTTER_LINES) and the payload
# of the catamaran reference (tests/data/hydrodynamics/generate_catamaran_mss.m 61-62).
PARAMETERS = {
    "length": 2.0,                                   # L, otter.m 92
    "beam": 1.08,                                    # B, otter.m 93
    "water_density": 1025.0,                         # rho, otter.m 91
    "hull_mass": 55.0,                               # m, otter.m 94
    "payload_mass": 25.0,                            # mp, generate_catamaran_mss.m 61
    "hull_center_of_gravity": [0.2, 0.0, -0.2],      # rg, otter.m 95
    "payload_position": [0.05, 0.0, -0.35],          # rp, generate_catamaran_mss.m 62
    "added_mass_coefficients": [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7],  # otter.m 153-158 (A11 scale 1)
    "radii_of_gyration": [0.4, 0.25, 0.25],          # R44/B, R55/L, R66/L, otter.m 96-98
}

# MATLAB reference recomputed by current MSS on 2026-10-05 (SOURCE.md). The
# legacy matlab_reference.csv is not read.
REFERENCE_CSV = "matlab_reference_mss_current.csv"


# Current of the reference generator, test_dynamics_consistency.m lines 14-15.
CURRENT_SPEED = 0.3
CURRENT_DIRECTION = np.deg2rad(30.0)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _contract():
    try:
        return importlib.import_module(CONTRACT_MODULE)
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed: {exc}")
        pytest.fail(f"block not ported yet: {exc}")


def _build(**overrides):
    """The block of the default form with the checked numbers, and its
    state-independent named outputs under the attribute names the assertions
    below use (``rigid_body_mass_matrix`` = output ``M_RB``, ...)."""
    block = _contract()
    values = block.check_rigid_body_values({**PARAMETERS, **overrides})
    function = block.rigid_body_casadi()
    out = function(nu=np.zeros(6), nu_r=np.zeros(6), **values)
    constants = SimpleNamespace(
        rigid_body_mass_matrix=np.array(out["M_RB"], dtype=float),
        added_mass_matrix=np.array(out["M_A"], dtype=float),
        total_mass_matrix=np.array(out["M"], dtype=float),
        inertia=np.array(out["inertia"], dtype=float),
        center_of_gravity=np.array(out["center_of_gravity"], dtype=float).ravel(),
        mass=float(out["mass"]),
    )
    return constants, (function, values)


def _evaluate(function, nu, nu_r):
    block, values = function
    out = block(nu=np.asarray(nu, float), nu_r=np.asarray(nu_r, float), **values)
    return {name: np.array(out[name], dtype=float) for name in ("M", "C_RB", "C_A")}


def _load_csv(name):
    path = DATA_DIR / name
    header = path.read_text().splitlines()[0].split(",")
    values = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return header, values


def _reference_cases():
    """The 50 frozen cases: nu, nu_r and the reference matrices."""
    in_header, inputs = _load_csv("inputs.csv")
    ref_header, ref = _load_csv(REFERENCE_CSV)
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
# G1 — block vs MATLAB reference
# --------------------------------------------------------------------------


def test_G1_block_matches_matlab_reference():
    """The block (one path) vs MATLAB running current MSS otter.m (2026-10-05)."""
    constants, function = _build()
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
# G2 — CasADi block: signature and frozen form
# --------------------------------------------------------------------------
def test_G2_function_signature():
    """The block takes the state and every declared parameter by name; with
    the numbers frozen in, it is a function of the state only."""
    from more_transformations.more_casadi_transformations import freeze

    constants, (block, values) = _build()
    declared = _contract().rigid_body_parameters()
    assert [d.name for d in declared] == list(PARAMETERS)
    assert block.name_in() == ["nu", "nu_r", *PARAMETERS]
    assert block.name_out() == ["M", "C_RB", "C_A", "M_RB", "M_A", "mass", "center_of_gravity", "inertia"]
    for d in declared:
        assert block.size_in(d.name) == d.shape, d.name
    function = freeze(block, declared, PARAMETERS)
    assert function.name_in() == ["nu", "nu_r"]
    assert function.name_out()[:3] == ["M", "C_RB", "C_A"]
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


def test_G2_block_matches_otter_m_transcription_on_random_states():
    """Second check of the whole matrices, the rotational block included:
    ``otter.m`` lines 122-159 transcribed below (``_otter_m_current_matrices``,
    ``_otter_m_crb``) on the 1000 seeded states."""
    _, function = _build()
    mrb, ma, ig, rg_total = _otter_m_current_matrices()
    nu, nu_r = _random_states()
    for k in range(len(nu)):
        out = _evaluate(function, nu[k], nu_r[k])
        assert _max_diff(out["M"], mrb + ma) <= G2_TOLERANCE, ("M", k)
        assert _max_diff(out["C_RB"], _otter_m_crb(ig, rg_total, nu[k][3:])) <= G2_TOLERANCE, ("C_RB", k)


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


def _worst_diff(function):
    """Largest deviation of M, C_RB, C_A from the current-MSS MATLAB reference."""
    ref = _reference_cases()
    worst = 0.0
    for k in range(len(ref["nu"])):
        out = _evaluate(function, ref["nu"][k], ref["nu_r"][k])
        expected = {"M": ref["M_RB"][k] + ref["M_A"][k], "C_RB": ref["C_RB"][k], "C_A": ref["C_A"][k]}
        for name, value in expected.items():
            worst = max(worst, _max_diff(out[name], value))
    return worst


@pytest.mark.parametrize("name", sorted(_perturbations()))
def test_G4_perturbed_block_is_detected(name):
    # Control: the unperturbed block matches the reference, so a detection
    # below is caused by the perturbation.
    _, control = _build()
    assert _worst_diff(control) <= G1_TOLERANCE, "control"
    _, function = _build(**_perturbations()[name])
    worst = _worst_diff(function)
    assert worst > G4_FACTOR * G1_TOLERANCE, (name, worst)


# --------------------------------------------------------------------------
# The temporary inertia flag is gone from the contract
# --------------------------------------------------------------------------
def test_contract_has_no_legacy_otter_inertia_flag():
    """``legacy_otter_inertia`` is no keyword, input or output of the block."""
    block = _contract()
    signature = inspect.signature(block.rigid_body_casadi)
    assert "legacy_otter_inertia" not in signature.parameters, "keyword still in the contract"
    _, (function, _) = _build()
    assert "legacy_otter_inertia" not in function.name_in() + function.name_out(), "still on the block"
    assert "legacy_otter_inertia" not in [d.name for d in block.rigid_body_parameters()]


# --------------------------------------------------------------------------
# G5 — Otter constants, MSS otter.m
# (<MSS_DIR>/CRAFT/USV/models/otter.m, MSS cc07579, 2026-10-07; every line is
# one higher than at 72656d1, release 2.0.2)
# --------------------------------------------------------------------------
OTTER_M_RELATIVE = "CRAFT/USV/models/otter.m"  # under MSS_DIR (nothing relative to one machine)
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
    mss = os.environ.get("MSS_DIR")
    if not mss:
        pytest.skip("MSS_DIR is not set (path to the MSS checkout); otter.m line check skipped")
    otter_m = Path(mss) / OTTER_M_RELATIVE
    if not otter_m.exists():
        pytest.fail(f"MSS_DIR={mss}: {OTTER_M_RELATIVE} not found")
    lines = otter_m.read_text().splitlines()
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
    # otter.m 142: MRB_CG(1:3,1:3) = (m+mp)*I3; H keeps this block (Hmtrx.m 17-18)
    np.testing.assert_allclose(
        constants.rigid_body_mass_matrix[:3, :3], (m + mp) * np.eye(3),
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m line 142",
    )
    # otter.m 152 + addedMassSurge.m 32-33: A11 = 2.7*rho*nabla^(5/3)/L^2, nabla = m/rho
    a11 = 2.7 * rho * (m / rho) ** (5.0 / 3.0) / L**2
    # otter.m 153-154, 159: MA = -diag(Xudot, Yvdot, Zwdot, ...)
    expected = [a11, OTTER["added"][0] * m, OTTER["added"][1] * m]
    np.testing.assert_allclose(
        np.diag(constants.added_mass_matrix)[:3], expected,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 152-154, 159",
    )


def _skew(v):
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def _otter_m_current_inertia():
    """otter.m lines 122-128 (revision 2026-04-20): inertia about the combined CG."""
    m, mp = OTTER["m"], OTTER_MP
    rg, rp = np.array(OTTER["rg"]), np.array(OTTER_RP)
    radii = np.array(OTTER["R"]) * np.array([OTTER["B"], OTTER["L"], OTTER["L"]])
    ig_cg = m * np.diag(radii**2)
    rg_total = (m * rg + mp * rp) / (m + mp)
    r_hull, r_payload = rg - rg_total, rp - rg_total
    ig = ig_cg - m * _skew(r_hull) @ _skew(r_hull) - mp * _skew(r_payload) @ _skew(r_payload)
    return ig, rg_total


def test_G5_otter_rotational_inertia_current_mss():
    # The block's one path (the inertia about the combined CG).
    constants, _ = _build(**_otter_inputs_as_parameters())
    ig, rg_total = _otter_m_current_inertia()
    m_total = OTTER["m"] + OTTER_MP
    h = np.block([[np.eye(3), _skew(rg_total).T], [np.zeros((3, 3)), np.eye(3)]])
    mrb_cg = np.block([[m_total * np.eye(3), np.zeros((3, 3))], [np.zeros((3, 3)), ig]])
    # otter.m lines 142-143, 147-148
    np.testing.assert_allclose(
        constants.rigid_body_mass_matrix, h.T @ mrb_cg @ h,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 122-128, 142-148",
    )
    # otter.m lines 155-157, 159
    expected = [c * ig[i, i] for i, c in enumerate(OTTER["added"][2:])]
    np.testing.assert_allclose(
        np.diag(constants.added_mass_matrix)[3:], expected,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 155-157, 159",
    )


def _otter_m_current_matrices():
    """MRB and MA of otter.m (current), on the frozen parameters.

    MRB: lines 142-143 (MRB_CG), 147-148 (H = Hmtrx(rg), MRB = H' MRB_CG H) with
    rg = rg_total (line 128). MA: lines 152-157, 159; Xudot from line 152 and
    addedMassSurge.m lines 32-33 (as in ``test_G5_otter_translational_constants``).
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
    """otter.m lines 144-145 (CRB_CG), 147, 149 (CRB = H' CRB_CG H)."""
    m_total = OTTER["m"] + OTTER_MP
    h = np.block([[np.eye(3), _skew(rg_total).T], [np.zeros((3, 3)), np.eye(3)]])
    crb_cg = np.block([
        [m_total * _skew(nu2), np.zeros((3, 3))],
        [np.zeros((3, 3)), -_skew(ig @ nu2)],
    ])
    return h.T @ crb_cg @ h


def test_G5_default_path_inertia_cg_and_rotational_added_mass():
    """The block (the corrected inertia) against MSS otter.m, computed here from the lines.

    Every constant comes from a line pinned by ``test_G5_cited_otter_lines_are_unchanged``
    (91-98 data, 123-129 inertia, 156-158 rotational added mass); the payload
    from the reference generator (``OTTER_MP``, ``OTTER_RP``).
    """
    constants, _ = _build(**_otter_inputs_as_parameters())
    m, mp = OTTER["m"], OTTER_MP
    rg_hull, rp = np.array(OTTER["rg"]), np.array(OTTER_RP)  # lines 94, 123; argument rp
    # lines 95-97: R44 = 0.4*B, R55 = 0.25*L, R66 = 0.25*L
    r44, r55, r66 = OTTER["R"][0] * OTTER["B"], OTTER["R"][1] * OTTER["L"], OTTER["R"][2] * OTTER["L"]
    ig_cg = m * np.diag([r44**2, r55**2, r66**2])                 # line 122
    rg_total = (m * rg_hull + mp * rp) / (m + mp)                  # line 124
    r_hull = rg_hull - rg_total                                    # line 125
    r_payload = rp - rg_total                                      # line 126
    s_hull, s_payload = _skew(r_hull), _skew(r_payload)
    ig = ig_cg - m * s_hull @ s_hull - mp * s_payload @ s_payload  # line 127 (Smtrx(.)^2)
    np.testing.assert_allclose(
        constants.inertia, ig, atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 122-127"
    )
    np.testing.assert_allclose(  # line 128: rg = rg_total, the CG used by H (line 147)
        constants.center_of_gravity, rg_total,
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 124, 128",
    )
    # lines 155-157: Kpdot = -0.2 Ig(1,1), Mqdot = -0.8 Ig(2,2), Nrdot = -1.7 Ig(3,3);
    # line 159: MA = -diag([..., Kpdot, Mqdot, Nrdot])
    k_pdot = -OTTER["added"][2] * ig[0, 0]
    m_qdot = -OTTER["added"][3] * ig[1, 1]
    n_rdot = -OTTER["added"][4] * ig[2, 2]
    np.testing.assert_allclose(
        np.diag(constants.added_mass_matrix)[3:], [-k_pdot, -m_qdot, -n_rdot],
        atol=G1_TOLERANCE, rtol=0.0, err_msg="otter.m lines 155-157, 159",
    )
