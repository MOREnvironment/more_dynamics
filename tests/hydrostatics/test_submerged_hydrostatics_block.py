"""Gate tests for the submerged hydrostatics block (restoring forces g(eta)).

Written before the block exists (job A-22, 2026-10-05); the block is ported by
job A-23. Every block test fails today with "block not ported yet".

Rewired by job A-30 (2026-10-05): the block's G1 and G4 read the MATLAB
reference regenerated with current MSS (job A-26,
``spheroid_matlab_reference_mss_current.csv``; byte-identical to the legacy
file, which the numpy source's G1 keeps reading). Added for A-27 pass A
finding 3: ``g`` is the left-hand-side restoring vector; the force on the
vehicle is ``-g`` (``test_restoring_sign_positive_buoyancy_lifts_the_vehicle``).

Completed by job U2a (2026-10-06), additions only: G5 against the printed
``gvect`` example of Fossen (2011) p. 61; the rotation comes from L0
``more_casadi_transformations`` (owner E-25, ADR 0003 §4.4), no local trig;
the block's docstring states the ``-g`` step (ADR 0003 U2 done-criterion,
A-27 pass A finding 3). The last two fail until the porter (U2b) changes
``submerged.py``.

Contract the porter must provide
--------------------------------
``more_dynamics.models.hydrostatics.submerged`` exports:

* ``preprocess_submerged_hydrostatics(weight, buoyancy, center_of_gravity,
  center_of_buoyancy) -> SubmergedHydrostaticsConstants`` (numpy). Weight and
  buoyancy in N; both centres in m, body frame, measured from the CO. No
  keyword has a default (generic by construction, owner E-16).
* ``SubmergedHydrostaticsConstants``: frozen dataclass with at least the four
  inputs as fields (``center_of_gravity`` and ``center_of_buoyancy`` as (3,)
  arrays).
* ``submerged_hydrostatics_casadi(constants) -> ca.Function`` with input
  ``["eta"]`` (6x1, ``[x y z phi theta psi]``, NED, rad) and output ``["g"]``
  (6x1, N and N m), the restoring vector of MSS ``gRvect.m`` with ``R`` the
  zyx Euler rotation of ``Rzyx.m`` (body to NED).

No departure from MSS exists in this block, so it has no flag.

Map rows covered (``agents-more/MIGRATION_MAP.md``), numpy source read in full
(paths relative to ``more_generic_models/more_generic_models/``):

* ``dynamics/plant/matrices/restoring_forces.py``: ``g_restoring_submerged``
  19-66 (``gvect.m`` form, Euler angles), ``g_restoring_neutral`` 71-117
  (W = B), ``g_restoring_R`` 121-167 (``gRvect.m`` form, row 3 of R).
* ``dynamics/plant/auv_spheroid/auv_spheroid.py``: ``get_restoring`` 241-261
  (passes ``r_cg``, not ``r_bg``, and ``R = euler_tf.R_bn``),
  ``compute_constant_values`` 193-194 (``W = mass * gravity``, ``B = W``).
* ``config/dataclass/plant/auv_spheroid_params.py`` 27-28 (``r_cg``, ``r_cb``).

Conventions found (hidden assumptions, see the ledger)
------------------------------------------------------
* NED, z down; ``eta[3:6]`` are zyx Euler angles; ``R(3,:) = [-s(th),
  c(th)s(phi), c(th)c(phi)]``; yaw does not enter g.
* ``r_bG``, ``r_bB`` are measured from the CO in BODY; the source vehicle
  passes ``r_cg`` and ``r_cb`` directly (it equals ``r_bg`` only because
  ``r_cb = 0``).
* ``W = m g(mu)`` with the WGS-84 latitude gravity of ``gravity.m`` at
  63.446827 deg; the stored reference is neutrally buoyant.

* ``g`` enters the equation of motion as ``M nu_dot + ... + g(eta) = tau``
  (``remus100.m`` line 256: ``... - D * nu_r - g``), so the
  applied force is ``-g``: for B > W at level attitude ``g_z = B - W > 0`` and
  the force ``-g_z`` points up (negative z in NED).

Frozen reference: ``tests/data/hydrostatics/`` (``SOURCE.md``).
"""

import ast
import importlib
import inspect
import os
import re
import sys
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "hydrostatics"
# Outside this repo only through environment variables (agents-more rule 9):
# MSS_DIR = an MSS checkout, MORE_GENERIC_MODELS_DIR = the more_generic_models
# repo root. Unset -> the tests that need them skip; the gates on the frozen
# CSVs in tests/data run everywhere.
CONTRACT_MODULE = "more_dynamics.models.hydrostatics.submerged"
# Current-MSS MATLAB reference (job A-26) for the block; legacy for the source.
REFERENCE_CSV = "spheroid_matlab_reference_mss_current.csv"
LEGACY_CSV = "spheroid_matlab_reference.csv"

G1_TOLERANCE = 1e-9   # 30_checks/README.md, gate G1
G2_TOLERANCE = 1e-10  # 30_checks/README.md, gate G2
G4_FACTOR = 10.0      # 30_checks/README.md, gate G4
SEED = 20261005
N_RANDOM_STATES = 1000

VEHICLE_NAMES = ("remus", "otter", "grethe", "marie", "hugin", "lauv",
                 "mariner", "torqeedo", "cybership", "prestero")

# MSS lines the constants and equations below are read from: (file, line) ->
# how the stripped line starts. ``test_cited_mss_lines_are_unchanged`` pins them.
REMUS = "CRAFT/AUV/models/remus100.m"
CITED_LINES = {
    (REMUS, 97): "mu = deg2rad(63.446827);",
    (REMUS, 132): "L_auv = 1.6;",
    (REMUS, 133): "D_auv = 0.19;",
    (REMUS, 135): "a = 1.0096 * L_auv/2;",
    (REMUS, 136): "b = 1.0096 * D_auv/2;",
    (REMUS, 138): "r_bG = [ 0 0 0.02 ]';",
    (REMUS, 139): "r_bB = [ 0 0 0 ]';",
    (REMUS, 212): "m = MRB(1,1); W = m * g_mu; B = W;",
    (REMUS, 229): "g = gRvect(W,B,R,r_bG,r_bB);",
    (REMUS, 256): "(tau + tau_liftdrag + tau_crossflow - C * nu_r - D * nu_r  - g)",
    ("LIBRARY/modeling/spheroid.m", 35): "rho = 1025;",
    ("LIBRARY/modeling/spheroid.m", 36): "m = 4/3 * pi * rho * a * b^2;",
    ("LIBRARY/modeling/gRvect.m", 27): "-(W-B) * R(3,1)",
    ("LIBRARY/modeling/gRvect.m", 28): "-(W-B) * R(3,2)",
    ("LIBRARY/modeling/gRvect.m", 29): "-(W-B) * R(3,3)",
    ("LIBRARY/modeling/gRvect.m", 30): "-(r_bG(2)*W - r_bB(2)*B) * R(3,3) + (r_bG(3)*W - r_bB(3)*B) * R(3,2)",
    ("LIBRARY/modeling/gRvect.m", 31): "-(r_bG(3)*W - r_bB(3)*B) * R(3,1) + (r_bG(1)*W - r_bB(1)*B) * R(3,3)",
    ("LIBRARY/modeling/gRvect.m", 32): "-(r_bG(1)*W - r_bB(1)*B) * R(3,2) + (r_bG(2)*W - r_bB(2)*B) * R(3,1) ];",
    ("LIBRARY/modeling/gvect.m", 23): "(W-B) * sth",
    ("LIBRARY/modeling/gvect.m", 24): "-(W-B) * cth * sphi",
    ("LIBRARY/modeling/gvect.m", 25): "-(W-B) * cth * cphi",
    ("LIBRARY/modeling/gvect.m", 26): "-(r_bG(2)*W-r_bB(2)*B) * cth * cphi + (r_bG(3)*W-r_bB(3)*B) * cth * sphi",
    ("LIBRARY/modeling/gvect.m", 27): "(r_bG(3)*W-r_bB(3)*B) * sth        + (r_bG(1)*W-r_bB(1)*B) * cth * cphi",
    ("LIBRARY/modeling/gvect.m", 28): "-(r_bG(1)*W-r_bB(1)*B) * cth * sphi - (r_bG(2)*W-r_bB(2)*B) * sth       ];",
    ("LIBRARY/kinematics/Rzyx.m", 19): "-sth      cth*sphi                  cth*cphi ];",
    ("INS/functions/gravity.m", 11): "g = 9.7803253359 * ( 1 + 0.001931850400 * sin(mu)^2 ) /...",
    ("INS/functions/gravity.m", 12): "sqrt( 1 - 0.006694384442 * sin(mu)^2 );",
}


# --------------------------------------------------------------------------
# Helpers: MSS text, contract, data
# --------------------------------------------------------------------------
def _env_dir(variable):
    value = os.environ.get(variable)
    if not value:
        pytest.skip(f"{variable} is not set (agents-more rule 9); set it to run this check")
    path = Path(value).expanduser()
    if not path.is_dir():
        pytest.skip(f"{variable}={value} is not a directory")
    return path


def _mss_line(rel, number):
    path = _env_dir("MSS_DIR") / rel
    if not path.exists():
        pytest.skip(f"{rel} not found under MSS_DIR")
    return path.read_text().splitlines()[number - 1].strip()


def _mss_value(rel, number, env=None):
    """Right-hand side of ``name = expr;`` on a pinned MSS line, evaluated."""
    line = _mss_line(rel, number)
    assert line.startswith(CITED_LINES[(rel, number)]), (rel, number, line)
    expr = line.split("=", 1)[1].split(";")[0].strip()
    if expr.startswith("["):
        return np.array([float(t) for t in expr.strip("[]' ").split()])
    expr = expr.replace("^", "**")
    names = {"__builtins__": {}, "deg2rad": np.deg2rad, "pi": np.pi, "sqrt": np.sqrt}
    return float(eval(expr, names, dict(env or {})))  # pinned MSS arithmetic only


def _gravity(mu):
    """INS/functions/gravity.m lines 11-12 (WGS-84 latitude gravity)."""
    for number in (11, 12):
        line = _mss_line("INS/functions/gravity.m", number)
        assert line.startswith(CITED_LINES[("INS/functions/gravity.m", number)]), line
    line = _mss_line("INS/functions/gravity.m", 11)
    nums = [float(t) for t in line.replace("(", " ").replace(")", " ").replace("*", " ").split()
            if t[0].isdigit()]
    g0, k1 = nums[0], nums[2]
    e2 = float(_mss_line("INS/functions/gravity.m", 12).split("-")[1].split("*")[0])
    return g0 * (1 + k1 * np.sin(mu) ** 2) / np.sqrt(1 - e2 * np.sin(mu) ** 2)


def _contract():
    try:
        return importlib.import_module(CONTRACT_MODULE)
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed (owner decision E-7): {exc}")
        pytest.fail(f"block not ported yet: {exc}")


def _build(params):
    block = _contract()
    constants = block.preprocess_submerged_hydrostatics(**params)
    return constants, block.submerged_hydrostatics_casadi(constants)


def _evaluate(function, eta):
    return np.array(function(eta=np.asarray(eta, float))["g"], dtype=float).reshape(-1)


def _load_csv(name):
    path = DATA_DIR / name
    header = path.read_text().splitlines()[0].split(",")
    return header, np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def _reference_cases(csv=REFERENCE_CSV):
    in_header, inputs = _load_csv("inputs.csv")
    ref_header, ref = _load_csv(csv)
    assert in_header[6:12] == [f"x{i}" for i in range(7, 13)]
    return {
        "eta": inputs[:, 6:12],
        "mass": ref[:, ref_header.index("M_RB_01")],
        "g": ref[:, [ref_header.index(f"g_{i:02d}") for i in range(1, 7)]],
    }


def _random_eta():
    rng = np.random.default_rng(SEED)
    eta = rng.uniform(-np.pi, np.pi, size=(N_RANDOM_STATES, 6))
    eta[:, :3] = rng.uniform(-50.0, 50.0, size=(N_RANDOM_STATES, 3))
    return eta


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _gravity_l0(mu_deg):
    """L0 gravity (owner E-22, E-25): the WGS-84 formula of gravity.m."""
    from more_transformations.ecef_ned_transforms import ECEFNEDtransform
    return ECEFNEDtransform.gravity(mu_deg)


# --------------------------------------------------------------------------
# Parameter sets
# --------------------------------------------------------------------------
# MSS remus100.m / spheroid.m constants, typed with their lines so the gates
# run without an MSS checkout; test_typed_constants_equal_the_cited_mss_lines
# reads them back when MSS_DIR is set.
REMUS_CONSTANTS = {
    "mu_deg": 63.446827,          # remus100.m line 97
    "L_auv": 1.6,                 # remus100.m line 132
    "D_auv": 0.19,                # remus100.m line 133
    "r_bG": [0.0, 0.0, 0.02],     # remus100.m line 138
    "r_bB": [0.0, 0.0, 0.0],      # remus100.m line 139
    "rho": 1025.0,                # spheroid.m line 35
}
def _remus_like_parameters():
    """MSS remus100.m: m from spheroid.m 35-36 with a, b from lines 132-136;
    W = B = m g(mu) (line 212, mu line 97); centres lines 138-139."""
    c = REMUS_CONSTANTS
    a = 1.0096 * c["L_auv"] / 2   # remus100.m line 135
    b = 1.0096 * c["D_auv"] / 2   # remus100.m line 136
    m = 4 / 3 * np.pi * c["rho"] * a * b ** 2   # spheroid.m line 36
    weight = m * _gravity_l0(c["mu_deg"])
    return {
        "weight": weight,
        "buoyancy": weight,
        "center_of_gravity": np.array(c["r_bG"]),
        "center_of_buoyancy": np.array(c["r_bB"]),
    }


def _scaled_parameters():
    """Second vehicle of the same class (test construction, not reference data):
    lengths x2, weight x8, 2 % positively buoyant, CG and CB moved off the axis."""
    base = _remus_like_parameters()
    return {
        "weight": 8.0 * base["weight"],
        "buoyancy": 8.0 * base["weight"] * 1.02,
        "center_of_gravity": 2.0 * base["center_of_gravity"] + np.array([0.02, -0.01, 0.0]),
        "center_of_buoyancy": 2.0 * base["center_of_buoyancy"] + np.array([0.01, 0.005, -0.01]),
    }


PARAMETER_SETS = {"remus_like": _remus_like_parameters, "scaled": _scaled_parameters}


def _matlab_parameters():
    """Weight exactly as the reference generator built it (SOURCE.md)."""
    ref = _reference_cases()
    assert np.ptp(ref["mass"]) == 0.0
    params = _remus_like_parameters()
    weight = ref["mass"][0] * _gravity_l0(REMUS_CONSTANTS["mu_deg"])
    return {**params, "weight": weight, "buoyancy": weight}


# --------------------------------------------------------------------------
# MSS transcriptions (each line cited and pinned)
# --------------------------------------------------------------------------
def _rzyx_row3(phi, theta):
    """Rzyx.m line 19: third row of R."""
    return np.array([-np.sin(theta), np.cos(theta) * np.sin(phi), np.cos(theta) * np.cos(phi)])


def _g_rvect(p, eta):
    """gRvect.m lines 26-32 with R from Rzyx.m."""
    w, b = p["weight"], p["buoyancy"]
    rg, rb = np.asarray(p["center_of_gravity"]), np.asarray(p["center_of_buoyancy"])
    r3 = _rzyx_row3(eta[3], eta[4])
    return np.array([
        -(w - b) * r3[0],
        -(w - b) * r3[1],
        -(w - b) * r3[2],
        -(rg[1] * w - rb[1] * b) * r3[2] + (rg[2] * w - rb[2] * b) * r3[1],
        -(rg[2] * w - rb[2] * b) * r3[0] + (rg[0] * w - rb[0] * b) * r3[2],
        -(rg[0] * w - rb[0] * b) * r3[1] + (rg[1] * w - rb[1] * b) * r3[0],
    ])


def _g_vect(p, eta):
    """gvect.m lines 19-28 (Euler-angle form)."""
    w, b = p["weight"], p["buoyancy"]
    rg, rb = np.asarray(p["center_of_gravity"]), np.asarray(p["center_of_buoyancy"])
    sth, cth, sphi, cphi = np.sin(eta[4]), np.cos(eta[4]), np.sin(eta[3]), np.cos(eta[3])
    return np.array([
        (w - b) * sth,
        -(w - b) * cth * sphi,
        -(w - b) * cth * cphi,
        -(rg[1] * w - rb[1] * b) * cth * cphi + (rg[2] * w - rb[2] * b) * cth * sphi,
        (rg[2] * w - rb[2] * b) * sth + (rg[0] * w - rb[0] * b) * cth * cphi,
        -(rg[0] * w - rb[0] * b) * cth * sphi - (rg[1] * w - rb[1] * b) * sth,
    ])


# --------------------------------------------------------------------------
# Numpy source
# --------------------------------------------------------------------------
def _numpy_source():
    root = str(_env_dir("MORE_GENERIC_MODELS_DIR"))
    if root not in sys.path:
        sys.path.insert(0, root)
    try:
        from more_generic_models.dynamics.plant.auv_spheroid.auv_spheroid import AUVSpheroid
        from more_generic_models.dynamics.plant.matrices.restoring_forces import RestoringForces
    except ImportError as exc:
        pytest.skip(f"numpy source not importable from MORE_GENERIC_MODELS_DIR={root}: {exc}")
    return AUVSpheroid(), RestoringForces


# --------------------------------------------------------------------------
# Reference sanity and pinned lines (no port needed)
# --------------------------------------------------------------------------
def test_reference_shapes_and_neutral_buoyancy():
    ref = _reference_cases()
    assert ref["eta"].shape == (50, 6) and ref["g"].shape == (50, 6)
    # B = W in the generator (line 261), so the three force rows are zero
    np.testing.assert_allclose(ref["g"][:, :3], 0.0, atol=G1_TOLERANCE, rtol=0.0)


def test_current_mss_reference_is_byte_identical_to_legacy():
    """SOURCE.md (A-26): current MSS changed nothing on this path."""
    assert (DATA_DIR / REFERENCE_CSV).read_bytes() == (DATA_DIR / LEGACY_CSV).read_bytes()


def test_typed_constants_equal_the_cited_mss_lines():
    """REMUS_CONSTANTS are typed with their lines (agents-more rules 4, 9);
    with MSS_DIR set they are read back from MSS, and the L0 gravity equals the
    transcription of gravity.m lines 11-12."""
    c = REMUS_CONSTANTS
    assert _mss_value(REMUS, 132) == c["L_auv"] and _mss_value(REMUS, 133) == c["D_auv"]
    assert _mss_value("LIBRARY/modeling/spheroid.m", 35) == c["rho"]
    assert abs(_mss_value(REMUS, 97) - np.deg2rad(c["mu_deg"])) == 0.0
    assert list(_mss_value(REMUS, 138)) == c["r_bG"] and list(_mss_value(REMUS, 139)) == c["r_bB"]
    for prefix, key in (("a = 1.0096", 135), ("b = 1.0096", 136)):
        assert _mss_line(REMUS, key).startswith(prefix)
    mu = np.deg2rad(c["mu_deg"])
    assert abs(_gravity(mu) - _gravity_l0(c["mu_deg"])) <= 1e-12


def test_cited_mss_lines_are_unchanged():
    for (rel, number), text in CITED_LINES.items():
        line = _mss_line(rel, number)
        assert line.startswith(text), (rel, number, line)


def test_MSS_transcriptions_agree_with_each_other():
    """gRvect.m (R row 3) and gvect.m (Euler) are the same vector; checks the
    transcriptions before they are used as references."""
    for name, make in PARAMETER_SETS.items():
        p = make()
        for eta in _random_eta()[:200]:
            assert _max_diff(_g_rvect(p, eta), _g_vect(p, eta)) <= G2_TOLERANCE, name


# --------------------------------------------------------------------------
# G1 — numpy source and block vs the stored MATLAB reference
# --------------------------------------------------------------------------
def test_G1_numpy_source_matches_matlab_reference():
    vehicle, _ = _numpy_source()
    ref = _reference_cases(LEGACY_CSV)   # the generator the source was checked on
    for k, eta in enumerate(ref["eta"]):
        assert _max_diff(vehicle.get_restoring(eta), ref["g"][k]) <= G1_TOLERANCE, k


def test_G1_block_matches_matlab_reference():
    """The block vs MATLAB running current MSS ``remus100.m`` / ``gRvect.m`` (A-26)."""
    _, function = _build(_matlab_parameters())
    ref = _reference_cases()
    for k, eta in enumerate(ref["eta"]):
        assert _max_diff(_evaluate(function, eta), ref["g"][k]) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G2 — CasADi block vs numpy source (both parameter sets)
# --------------------------------------------------------------------------
def test_G2_function_signature():
    constants, function = _build(_remus_like_parameters())
    assert function.name_in() == ["eta"]
    assert function.name_out() == ["g"]
    assert function.size_in(0) == (6, 1)
    assert function.size_out(0) == (6, 1)
    for field in ("weight", "buoyancy", "center_of_gravity", "center_of_buoyancy"):
        assert hasattr(constants, field), field


@pytest.mark.parametrize("parameter_set", sorted(PARAMETER_SETS))
def test_G2_block_matches_numpy_source(parameter_set):
    vehicle, restoring = _numpy_source()
    p = PARAMETER_SETS[parameter_set]()
    _, function = _build(p)
    w, b = p["weight"], p["buoyancy"]
    rg, rb = p["center_of_gravity"], p["center_of_buoyancy"]
    states = np.vstack([_reference_cases()["eta"], _random_eta()])
    for k, eta in enumerate(states):
        out = _evaluate(function, eta)
        r_bn = vehicle.euler_tf.R_bn(*eta[3:6])
        assert _max_diff(out, restoring.g_restoring_R(w, b, eta, rg, rb, r_bn)) <= G2_TOLERANCE, k
        assert _max_diff(out, restoring.g_restoring_submerged(w, b, eta, rg, rb)) <= G2_TOLERANCE, k
        if w == b:
            assert _max_diff(out, restoring.g_restoring_neutral(w, eta, rg, rb)) <= G2_TOLERANCE, k


# --------------------------------------------------------------------------
# G4 — a broken model must fail
# --------------------------------------------------------------------------
def _perturbations():
    return {
        "buoyancy_plus_1_percent": lambda p: {**p, "buoyancy": 1.01 * p["buoyancy"]},
        "cg_z_plus_1_cm": lambda p: {
            **p, "center_of_gravity": np.asarray(p["center_of_gravity"]) + [0.0, 0.0, 0.01]},
    }


def _worst_vs_reference(function):
    ref = _reference_cases()
    return max(_max_diff(_evaluate(function, eta), ref["g"][k]) for k, eta in enumerate(ref["eta"]))


@pytest.mark.parametrize("name", sorted(_perturbations()))
def test_G4_perturbed_block_is_detected(name):
    p = _matlab_parameters()
    _, control = _build(p)
    assert _worst_vs_reference(control) <= G1_TOLERANCE, "control"
    _, function = _build(_perturbations()[name](p))
    worst = _worst_vs_reference(function)
    assert worst > G4_FACTOR * G1_TOLERANCE, (name, worst)


# --------------------------------------------------------------------------
# Headline (owner E-18, E-20 Q2 a): the default block equals MSS
# --------------------------------------------------------------------------
@pytest.mark.parametrize("parameter_set", sorted(PARAMETER_SETS))
def test_MSS_default_block_equals_gRvect(parameter_set):
    p = PARAMETER_SETS[parameter_set]()
    _, function = _build(p)
    for k, eta in enumerate(_random_eta()):
        expected = _g_rvect(p, eta)  # gRvect.m 26-32, Rzyx.m 19
        assert _max_diff(_evaluate(function, eta), expected) <= G1_TOLERANCE, (
            "gRvect.m lines 26-32", parameter_set, k)


# --------------------------------------------------------------------------
# Sign of the restoring vector (A-27 pass A finding 3)
# --------------------------------------------------------------------------
def test_restoring_sign_positive_buoyancy_lifts_the_vehicle():
    """B > W at level attitude: the block returns ``g_z = B - W > 0``
    (``gRvect.m`` line 29 with ``R(3,3) = 1``); the force on the vehicle is
    ``-g_z < 0``, i.e. upward in NED. Any position and yaw (they do not enter)."""
    p = _scaled_parameters()
    assert p["buoyancy"] > p["weight"]
    _, function = _build(p)
    rng = np.random.default_rng(SEED)
    for k in range(100):
        eta = np.concatenate([rng.uniform(-50.0, 50.0, 3), [0.0, 0.0, rng.uniform(-np.pi, np.pi)]])
        g = _evaluate(function, eta)
        assert g[2] > 0.0, (k, g)
        assert abs(g[2] - (p["buoyancy"] - p["weight"])) <= G1_TOLERANCE, (k, g)
        applied_force_z = -g[2]
        assert applied_force_z < 0.0, (k, "upward in NED is negative z")
        assert _max_diff(g[:2], 0.0) <= G1_TOLERANCE, (k, g)


# --------------------------------------------------------------------------
# G5 — published number (Fossen 2011, Section 4.1.1, p. 61; read in full)
# --------------------------------------------------------------------------
def test_G5_fossen_2011_page_61_gvect_example():
    """The book's MSS example: m = 1000, g = 9.81, W = B, r_g = [0 0 0],
    r_b = [0 0 -10], phi = 30 deg; printed g = 1e4 [0, 0, 0, 1.8324,
    9.0997, 0]. The book's script sets ``theta = 10*(180/pi)`` (a unit slip:
    572.96 rad, not 10 deg); the printed numbers are that input's, so the
    test enters theta exactly as printed. It checks the formula, not the
    angle (with 10 deg the moments would be 4.8305e4 and 1.7035e4)."""
    weight = 1000 * 9.81
    params = {"weight": weight, "buoyancy": weight,
              "center_of_gravity": np.array([0.0, 0.0, 0.0]),
              "center_of_buoyancy": np.array([0.0, 0.0, -10.0])}
    _, function = _build(params)
    eta = np.array([0.0, 0.0, 0.0, 30 * (np.pi / 180), 10 * (180 / np.pi), 0.0])
    printed = 1e4 * np.array([0.0, 0.0, 0.0, 1.8324, 9.0997, 0.0])
    half_last_digit = 0.5e-4 * 1e4
    assert _max_diff(_evaluate(function, eta), printed) <= half_last_digit


# --------------------------------------------------------------------------
# Transforms from L0 (owner E-25, ADR 0003 §4.4) and the -g statement
# --------------------------------------------------------------------------
L0_ROTATIONS = {
    "more_transformations.more_casadi_transformations.matrix_transforms.MatrixTransforms":
        ("Rzyx_explicit", "Rzyx", "Rzyx_row3"),
    "more_transformations.more_casadi_transformations.euler_ned_body_transforms.EulerNEDBodyTransforms":
        ("R_fossen", "R_bn", "R_nb"),
}


def test_L0_rotation_is_taken_from_more_casadi_transformations(monkeypatch):
    """The graph's ``R(3,:)`` comes from the CasADi L0 (``Rzyx_row3``, ``Rzyx``
    or ``R_bn``), counted by wrapping those functions while the block is built."""
    calls = []
    for path, names in L0_ROTATIONS.items():
        module_name, cls_name = path.rsplit(".", 1)
        cls = getattr(importlib.import_module(module_name), cls_name)
        for name in names:
            original = getattr(cls, name)

            def spy(*args, _original=original, _name=name, **kwargs):
                calls.append(_name)
                return _original(*args, **kwargs)
            monkeypatch.setattr(cls, name, staticmethod(spy))
    _build(_remus_like_parameters())
    assert calls, "no L0 rotation was called; R(3,:) is computed locally"


def test_L0_no_local_trigonometry_in_the_block():
    block = _contract()
    tree = ast.parse(_code_without_docstrings(block))
    trig = [ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute) and n.func.attr in {"sin", "cos", "tan"}]
    assert not trig, f"local rotation terms: {trig}"


def test_docstring_states_the_minus_g_step():
    """ADR 0003 U2: "the -g step documented" (A-27 pass A finding 3). The
    module or ``submerged_hydrostatics_casadi`` docstring must say that the
    applied force is ``-g`` (equation of motion ``... + g = tau``)."""
    block = _contract()
    text = " ".join(filter(None, [block.__doc__, block.submerged_hydrostatics_casadi.__doc__]))
    text = re.sub(r"-\s+g\b", "-g", text.replace("−", "-"))
    assert re.search(r"(?<![\w-])-g(?![\w(])", text), "the docstring does not state the -g step"


# --------------------------------------------------------------------------
# Generic by construction (owner E-16)
# --------------------------------------------------------------------------
def _code_without_docstrings(module):
    tree = ast.parse(Path(module.__file__).read_text())
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(body, list) and body and isinstance(body[0], ast.Expr) \
                and isinstance(getattr(body[0], "value", None), ast.Constant) \
                and isinstance(body[0].value.value, str):
            body.pop(0)
    return ast.unparse(tree).lower()


def test_generic_block_reads_no_vehicle_name_and_has_no_vehicle_defaults():
    block = _contract()
    code = _code_without_docstrings(block)
    assert not [n for n in VEHICLE_NAMES if n in code]
    signature = inspect.signature(block.preprocess_submerged_hydrostatics)
    defaults = [n for n, p in signature.parameters.items() if p.default is not inspect.Parameter.empty]
    assert defaults == [], f"vehicle values must come from the parameter set: {defaults}"
    # two parameter sets give two different blocks from the same code
    _, a = _build(_remus_like_parameters())
    _, b = _build(_scaled_parameters())
    eta = _random_eta()[0]
    assert _max_diff(_evaluate(a, eta), _evaluate(b, eta)) > G4_FACTOR * G1_TOLERANCE
