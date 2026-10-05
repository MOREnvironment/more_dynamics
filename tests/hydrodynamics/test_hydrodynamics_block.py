"""Gate tests for the hydrodynamics blocks: linear damping, cross-flow drag, lift/drag.

Written before the blocks exist (job A-22, 2026-10-05); the blocks are ported
by job A-24. Every block test fails today with "block not ported yet".

Rewired by job A-30 (2026-10-05): section "G1-MSS" compares the default paths
with MATLAB running current MSS (job A-26, ``*_mss_current.csv``); where the
default deviates by register (``30_checks/README.md`` D-MSS-1, Reynolds number
on the diameter) the MATLAB file is compared with the MSS flag
(``cross_flow_reynolds_length="length"``) and the default with the
transcription plus the register's stated effect. D-MSS-2 (owner E-24): Hoerner
below its table raises ``ValueError``. The test-side transcriptions stay as the
second check and are themselves checked against the new files.

Owner rulings these tests encode
--------------------------------
* E-18 / E-20 Q2 a: the **default** path equals current MSS; each template
  departure stays behind a named flag that reproduces the old source value.
* E-10 a: cross-flow drag on MSS's 20 strip midpoints with the full-precision
  ``cylinderDrag.m`` table (default); the 21-end-point grid behind a flag.
* E-20 Q3 b: the cross-flow Reynolds number uses the **diameter** by default;
  ``cross_flow_reynolds_length="length"`` reproduces ``cylinderDrag.m:77``.
* E-16: generic by construction (no vehicle name or number inside a block).

Contract the porter must provide
--------------------------------
``more_dynamics.models.hydrodynamics.linear_damping``:

* ``preprocess_submerged_linear_damping(rigid_body_mass_matrix,
  added_mass_matrix, weight, center_of_gravity, center_of_buoyancy,
  time_constants, damping_ratios, sway_damping_fade=False)``. Matrices 6x6,
  weight N, centres (3,) m from the CO, ``time_constants = (T1, T2, T6)`` s,
  ``damping_ratios = (zeta4, zeta5)``. Constants expose
  ``damping_coefficients`` (6,): the ``Dmtrx.m`` diagonal before any fade.
* ``submerged_linear_damping_casadi(constants)``: input ``["nu_r"]`` (6x1),
  outputs ``["D", "tau"]``: ``D`` (6x6, positive diagonal as ``Dmtrx.m``) with
  ``D[0,0] *= exp(-3 U_r)`` (``remus100.m:216``), ``U_r = |nu_r[0:3]|``, and
  ``D[1,1]`` faded the same way only when ``sway_damping_fade`` is True (the
  template departure, A-18); ``tau = -D @ nu_r``.
* ``preprocess_surface_linear_damping(mass_matrix, restoring_matrix,
  max_forward_thrust, max_speed, time_constants, damping_ratios,
  yaw_damping_nonlinearity)``. ``mass_matrix`` = M_RB + M_A (6x6);
  ``restoring_matrix`` = G about the centre of flotation (6x6, only G33, G44,
  G55 are used); ``max_forward_thrust`` N; ``max_speed`` m/s;
  ``time_constants = (T_sway, T_yaw)``; ``damping_ratios = (zeta3, zeta4,
  zeta5)``; ``yaw_damping_nonlinearity`` = the 10 of ``otter.m:240``.
  Constants expose ``damping_derivatives`` (6,) = ``[Xu Yv Zw Kp Mq Nr]``,
  the negative numbers of ``otter.m`` 202-207.
* ``surface_linear_damping_casadi(constants)``: input ``["nu_r"]``, outputs
  ``["D", "tau"]``: ``D = -diag(damping_derivatives)`` (positive, the same
  sign as the submerged ``D``) and ``tau`` = ``otter.m`` 235-242 (the force
  on the vehicle, quadratic yaw term included).

``more_dynamics.models.hydrodynamics.cross_flow``:

* ``preprocess_cross_flow_drag(length, beam, draft, water_density,
  drag_model, cross_flow_reynolds_length="diameter", strip_grid="midpoint")``;
  ``drag_model`` in {"cylinder", "hoerner"} (no default), flag values
  ``cross_flow_reynolds_length`` in {"diameter", "length"} and ``strip_grid``
  in {"midpoint", "endpoint"}. "diameter" is the ``beam`` argument (MSS calls
  ``crossFlowDrag(L, D, D, ...)``).
* ``cross_flow_drag_casadi(constants)``: input ``["nu_r"]``, output ``["tau"]``.

``more_dynamics.models.hydrodynamics.lift_drag``:

* ``preprocess_lift_drag(span, planform_area, parasitic_drag_coefficient,
  oswald_efficiency, water_density)``.
* ``lift_drag_casadi(constants)``: input ``["nu_r"]``, output ``["tau"]``;
  ``alpha = atan2(w_r, u_r)``, ``U_r = |nu_r[0:3]|`` (``remus100.m`` 127-128).

Map rows covered, numpy source read in full (``more_generic_models/more_generic_models/``)
---------------------------------------------------------------------------------------
* ``dynamics/plant/matrices/linear_damping.py``: ``get_damping_coeff`` 16-84,
  ``D_linear`` 88-112 (fades surge **and sway**), ``D_linear_catamaran``
  165-197, ``get_tau_damping`` 200-226. ``get_D_linear_surface_vessel`` is
  already ported (``linear_surface.py``); ``force_surge_damping`` 230-294 is
  not on either vehicle's path and is not tested here.
* ``dynamics/plant/matrices/hydrodynamics.py``: ``get_lift_drag_coeff``
  62-75, ``force_lift_drag`` 78-88, ``_hoerner`` 95-100, ``_cylinder_drag``
  103-120 (Re column rounded to 4 digits, A-2), ``cross_flow_drag`` 127-178,
  ``_get_strips`` 181-190 (21 end points).
* ``dynamics/plant/matrices/drag_models.py``: ``hoerner`` 16-64,
  ``cylinder_drag`` 68-162 (same rounded table; not on any vehicle's path).
* ``dynamics/plant/auv_spheroid/auv_spheroid.py``: ``compute_constant_values``
  157-197, ``get_D`` 235-236, ``get_tau_lift_drag`` 269-294,
  ``get_tau_cross_flow`` 296-323.
* ``dynamics/plant/asv_catamaran/asv_catamaran.py``: ``compute_constant_values``
  121-186, ``get_tau_damping`` 232-233, ``get_tau_cross_flow`` 235-244.

Conventions found (hidden assumptions, see the ledger)
------------------------------------------------------
* Two damping sign conventions: ``Dmtrx.m`` positive D subtracted
  (``- D*nu_r``); ``otter.m`` negative derivatives added (``+ tau_damp``).
  The contract outputs both as a positive ``D`` and a force ``tau``.
* The cross-flow coefficient is one number per call, set by the mid-body
  cross-flow speed ``sqrt(v_r^2 + w_r^2)``; kinematic viscosity 1e-6 m^2/s
  is inside ``Re = U * L * 1e6``; strip height = ``draft``.
* ``crossFlowDrag.m`` hard-codes rho = 1025, ``forceLiftDrag.m`` rho = 1026
  and ``coeffLiftDrag.m`` e = 0.3; the blocks take both as parameters and the
  MSS tests pass the hard-coded values.

Frozen reference: ``tests/data/hydrodynamics/`` (``SOURCE.md``).
"""

import ast
import importlib
import inspect
import re
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "hydrodynamics"
MORE = Path(__file__).resolve().parents[4]
LINEAR_DAMPING = "more_dynamics.models.hydrodynamics.linear_damping"
CROSS_FLOW = "more_dynamics.models.hydrodynamics.cross_flow"
LIFT_DRAG = "more_dynamics.models.hydrodynamics.lift_drag"

G1_TOLERANCE = 1e-9   # 30_checks/README.md, gate G1
G2_TOLERANCE = 1e-10  # 30_checks/README.md, gate G2
G4_FACTOR = 10.0      # 30_checks/README.md, gate G4
SEED = 20261005
N_RANDOM_STATES = 1000
# A-2 ledger: the rounded Re table moves the spheroid cross-flow force by 2.3e-2 N
ROUNDED_TABLE_RESIDUAL_BOUND = 0.03

VEHICLE_NAMES = ("remus", "otter", "grethe", "marie", "hugin", "lauv",
                 "mariner", "torqeedo", "cybership", "prestero")

REMUS = "source-sim/MSS/CRAFT/AUV/models/remus100.m"
OTTER = "source-sim/MSS/CRAFT/USV/models/otter.m"
DMTRX = "source-sim/MSS/LIBRARY/modeling/Dmtrx.m"
FLD = "source-sim/MSS/LIBRARY/modeling/forceLiftDrag.m"
CLD = "source-sim/MSS/LIBRARY/modeling/coeffLiftDrag.m"
CFD = "source-sim/MSS/LIBRARY/modeling/crossFlowDrag.m"
CYL = "source-sim/MSS/HYDRO/cylinderDrag.m"
HOERNER = "source-sim/MSS/HYDRO/Hoerner.m"
GRAVITY = "source-sim/MSS/INS/functions/gravity.m"
SPHEROID_M = "source-sim/MSS/LIBRARY/modeling/spheroid.m"
GEN_SPH = "more_generic_models/more_generic_models/test/plant/auv_spheroid/test_mss_matlab/test_dynamics_consitency.m"
GEN_CAT = "more_generic_models/more_generic_models/test/plant/asv_catamaran/test_mss_matlab/test_dynamics_consistency.m"

# (file relative to <more>, line) -> how the stripped line starts. Pinned by
# test_cited_lines_are_unchanged; every constant below is parsed from these.
CITED_LINES = {
    (REMUS, 97): "mu = deg2rad(63.446827);",
    (REMUS, 127): "alpha = atan2( nu_r(3), nu_r(1) );",
    (REMUS, 128): "U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );",
    (REMUS, 132): "L_auv = 1.6;",
    (REMUS, 133): "D_auv = 0.19;",
    (REMUS, 134): "S = 0.7 * L_auv * D_auv;",
    (REMUS, 136): "b = 1.0096 * D_auv/2;",
    (REMUS, 138): "r_bG = [ 0 0 0.02 ]';",
    (REMUS, 139): "r_bB = [ 0 0 0 ]';",
    (REMUS, 144): "Cd = 0.42;",
    (REMUS, 145): "CD_0 = Cd * pi * b^2 / S;",
    (REMUS, 193): "T1 = 20;",
    (REMUS, 194): "T2 = 20;",
    (REMUS, 195): "zeta4 = 0.3;",
    (REMUS, 196): "zeta5 = 0.8;",
    (REMUS, 197): "T6 = 1;",
    (REMUS, 215): "D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bG' r_bB']);",
    (REMUS, 216): "D(1,1) = D(1,1) * exp(-3 * U_r);",
    (REMUS, 218): "tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);",
    (REMUS, 219): "tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r,'cylinder');",
    (DMTRX, 30): "M = MRB + MA;",
    (DMTRX, 44): "T3 = T2;",
    (DMTRX, 45): "w4 = sqrt( W * (r_bg(3)-r_bb(3)) / M(4,4) );",
    (DMTRX, 46): "w5 = sqrt( W * (r_bg(3)-r_bb(3)) / M(5,5) );",
    (DMTRX, 48): "D = diag( [M(1,1)/T1 M(2,2)/T2 M(3,3)/T3...",
    (DMTRX, 49): "M(4,4)*2*zeta4*w4  M(5,5)*2*zeta5*w5 M(6,6)/T6 ] );",
    (FLD, 26): "rho = 1026;",
    (FLD, 28): "[CL,CD] = coeffLiftDrag(b,S,CD_0,alpha,0);",
    (FLD, 30): "F_drag = 1/2 * rho * U_r^2 * S * CD;",
    (FLD, 31): "F_lift = 1/2 * rho * U_r^2 * S * CL;",
    (FLD, 35): "cos(alpha) * (-F_drag) - sin(alpha) * (-F_lift)",
    (FLD, 37): "sin(alpha) * (-F_drag) + cos(alpha) * (-F_lift)",
    (CLD, 54): "e = 0.3;",
    (CLD, 57): "AR = b^2/S;",
    (CLD, 60): "CL_alpha = pi * AR / ( 1 + sqrt(1 + (AR/2)^2) );",
    (CLD, 61): "CL_linear = CL_alpha * alpha;",
    (CLD, 68): "CD = CD_0 + CL.^2 / (pi * e * AR);",
    (CFD, 24): "rho = 1025;",
    (CFD, 25): "nStrips = 20;",
    (CFD, 26): "dx = L/nStrips;",
    (CFD, 41): "xL = -L/2 + (i - 0.5) * dx;",
    (CFD, 46): "U_h = abs(v_r + xL * r) * (v_r + xL * r);",
    (CFD, 47): "U_v = abs(w_r + xL * q) * (w_r + xL * q);",
    (CFD, 48): "Yh = Yh - 0.5 * rho * T * Cd_2D * U_h * dx;",
    (CFD, 49): "Zh = Zh - 0.5 * rho * T * Cd_2D * U_v * dx;",
    (CFD, 50): "Mh = Mh - 0.5 * rho * T * Cd_2D * xL * U_v * dx;",
    (CFD, 51): "Nh = Nh - 0.5 * rho * T * Cd_2D * xL * U_h * dx;",
    (CYL, 76): "U_crossflow = sqrt(nu_r(2)^2+nu_r(3)^2);",
    (CYL, 77): "Re = U_crossflow * L * 1e6;",
    (CYL, 89): "if Re < 2e5",
    (HOERNER, 47): "if B/(2*T) <= 4.00309",
    (HOERNER, 48): "CY_2D = interp1(CD_DATA(:,1),CD_DATA(:,2),B/(2*T));",
    (HOERNER, 50): "CY_2D = 0.559315;",
    (GRAVITY, 11): "g = 9.7803253359 * ( 1 + 0.001931850400 * sin(mu)^2 ) /...",
    (GRAVITY, 12): "sqrt( 1 - 0.006694384442 * sin(mu)^2 );",
    (OTTER, 90): "g   = 9.81;",
    (OTTER, 91): "rho = 1025;",
    (OTTER, 92): "L = 2.0;",
    (OTTER, 94): "m = 55.0;",
    (OTTER, 99): "T_sway = 1;",
    (OTTER, 100): "T_yaw = 1;",
    (OTTER, 101): "Umax = 6 * 0.5144;",
    (OTTER, 104): "B_pont  = 0.25;",
    (OTTER, 107): "Cb_pont = 0.4;",
    (OTTER, 121): "nabla = (m+mp)/rho;",
    (OTTER, 122): "T = nabla / (2 * Cb_pont * B_pont * L);",
    (OTTER, 192): "LCF = -0.2;",
    (OTTER, 197): "w3 = sqrt( G33/M(3,3) );",
    (OTTER, 198): "w4 = sqrt( G44/M(4,4) );",
    (OTTER, 199): "w5 = sqrt( G55/M(5,5) );",
    (OTTER, 202): "Xu = -24.4 * g / Umax;",
    (OTTER, 203): "Yv = -M(2,2) / T_sway;",
    (OTTER, 204): "Zw = -2 * 0.3 * w3 * M(3,3);",
    (OTTER, 205): "Kp = -2 * 0.2 * w4 * M(4,4);",
    (OTTER, 206): "Mq = -2 * 0.4 * w5 * M(5,5);",
    (OTTER, 207): "Nr = -M(6,6) / T_yaw;",
    (OTTER, 235): "Xh = Xu * nu_r(1);",
    (OTTER, 236): "Yh = Yv * nu_r(2);",
    (OTTER, 237): "Zh = Zw * nu_r(3);",
    (OTTER, 238): "Kh = Kp * nu_r(4);",
    (OTTER, 239): "Mh = Mq * nu_r(5);",
    (OTTER, 240): "Nh = Nr * (1 + 10 * abs(nu_r(6))) * nu_r(6);",
    (OTTER, 245): "tau_crossflow = crossFlowDrag(L,B_pont,T,nu_r);",
    (GEN_SPH, 206): "rho = 1025;",
    (GEN_SPH, 267): "D(1,1) = D(1,1) * exp(-3*U_r);",
    (GEN_SPH, 268): "D(2,2) = D(2,2) * exp(-3*U_r);",
    (GEN_SPH, 467): "rho = 1025;",
    (GEN_CAT, 12): "mp = 25;",
    (GEN_CAT, 136): "mu = deg2rad(63.446827);",
    (GEN_CAT, 137): "g = gravity(mu);",
    (GEN_CAT, 187): "Xu = -24.4*g/Umax; Yv = -M(2,2)/T_sway; Zw = -2*0.3*w3*M(3,3);",
}


# --------------------------------------------------------------------------
# Helpers: pinned text, contract, data
# --------------------------------------------------------------------------
def _line(rel, number):
    path = MORE / rel
    if not path.exists():
        pytest.skip(f"cited file not found: {path}")
    line = path.read_text().splitlines()[number - 1].strip()
    if (rel, number) in CITED_LINES:
        assert line.startswith(CITED_LINES[(rel, number)]), (rel, number, line)
    return line


def _value(rel, number, env=None):
    """Right-hand side of ``name = expr;`` on a pinned line, evaluated."""
    expr = _line(rel, number).split("=", 1)[1].split(";")[0].strip()
    if expr.startswith("["):
        return np.array([float(t) for t in expr.strip("[]' ").split()])
    names = {"__builtins__": {}, "deg2rad": np.deg2rad, "pi": np.pi, "sqrt": np.sqrt}
    return float(eval(expr.replace("^", "**"), names, dict(env or {})))  # pinned arithmetic


def _number_in(rel, number, pattern):
    match = re.search(pattern, _line(rel, number))
    assert match, (rel, number, pattern)
    return float(match.group(1))


def _table(rel, name):
    """A ``NAME = [ ... ];`` numeric table of an MSS file, parsed from the text."""
    path = MORE / rel
    if not path.exists():
        pytest.skip(f"cited file not found: {path}")
    block = re.split(rf"(?m)^{name} = \[", path.read_text())[1].split("];")[0]
    rows = [[float(t) for t in ln.replace("...", "").split()] for ln in block.splitlines()]
    return np.array([r for r in rows if r])


def _gravity(mu):
    """INS/functions/gravity.m lines 11-12."""
    g0 = _number_in(GRAVITY, 11, r"g = ([0-9.]+) \*")
    k1 = _number_in(GRAVITY, 11, r"1 \+ ([0-9.]+) \*")
    e2 = _number_in(GRAVITY, 12, r"1 - ([0-9.]+) \*")
    return g0 * (1 + k1 * np.sin(mu) ** 2) / np.sqrt(1 - e2 * np.sin(mu) ** 2)


def _contract(name):
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed (owner decision E-7): {exc}")
        pytest.fail(f"block not ported yet: {exc}")


def _evaluate(function, nu_r):
    out = function(nu_r=np.asarray(nu_r, float))
    return {k: np.array(v, dtype=float) for k, v in out.items()}


def _load_csv(name):
    path = DATA_DIR / name
    header = path.read_text().splitlines()[0].split(",")
    return header, np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def _columns(header, values, prefix, n):
    return values[:, [header.index(f"{prefix}_{i:02d}") for i in range(1, n + 1)]]


LEGACY_SPHEROID = "spheroid_matlab_reference.csv"
LEGACY_CATAMARAN = "catamaran_matlab_reference.csv"
# MATLAB running current MSS (job A-26, SOURCE.md): the default-path references
MSS_SPHEROID = "spheroid_matlab_reference_mss_current.csv"
MSS_CATAMARAN = "catamaran_matlab_reference_mss_current.csv"


def _spheroid_reference(csv=LEGACY_SPHEROID):
    h, v = _load_csv(csv)
    return {
        "nu_r": _columns(h, v, "nu_r", 6),
        "alpha": v[:, h.index("alpha")],
        "U_r": v[:, h.index("U_r")],
        "tau_lift_drag": _columns(h, v, "tau_lift_drag", 6),
        "tau_crossflow": _columns(h, v, "tau_crossflow", 6),
        "M_RB": _columns(h, v, "M_RB", 36).reshape(-1, 6, 6),  # row-major
        "M_A": _columns(h, v, "M_A", 36).reshape(-1, 6, 6),
        "D": _columns(h, v, "D", 36).reshape(-1, 6, 6),
    }


def _catamaran_reference(csv=LEGACY_CATAMARAN):
    h, v = _load_csv(csv)
    return {
        "nu_r": _columns(h, v, "nu_r", 6),
        "tau_damp": _columns(h, v, "tau_damp", 6),
        "tau_crossflow": _columns(h, v, "tau_crossflow", 6),
        "M": _columns(h, v, "M_total", 36).reshape(-1, 6, 6),
        "G": _columns(h, v, "G", 36).reshape(-1, 6, 6),
        "D": _columns(h, v, "D", 36).reshape(-1, 6, 6),
    }


def _random_nu_r():
    """1000 seeded states: 700 uniform in +-3, 300 with slow cross-flow
    (|v|, |w| < 0.1 m/s, Re < 2e5 on the 1.6 m length: the sub-critical
    branch the 50 reference cases never reach, A-2 section 7)."""
    rng = np.random.default_rng(SEED)
    fast = rng.uniform(-3.0, 3.0, size=(700, 6))
    slow = rng.uniform(-3.0, 3.0, size=(300, 6))
    slow[:, 1:3] = rng.uniform(-0.1, 0.1, size=(300, 2))
    slow[:, 3:6] = rng.uniform(-0.05, 0.05, size=(300, 3))
    return np.vstack([fast, slow])


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _skew(v):
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def _hmtrx(r):
    return np.block([[np.eye(3), _skew(r).T], [np.zeros((3, 3)), np.eye(3)]])


# --------------------------------------------------------------------------
# Parameter sets (every number parsed from a pinned line)
# --------------------------------------------------------------------------
def _spheroid_geometry():
    l_auv, d_auv = _value(REMUS, 132), _value(REMUS, 133)
    s = _value(REMUS, 134, {"L_auv": l_auv, "D_auv": d_auv})
    b = _value(REMUS, 136, {"D_auv": d_auv})
    cd0 = _value(REMUS, 145, {"Cd": _value(REMUS, 144), "b": b, "S": s})
    return {"L": l_auv, "D": d_auv, "S": s, "CD_0": cd0}


def _spheroid_damping(fade, csv=LEGACY_SPHEROID):
    """Matrices and weight as the spheroid reference generator built them."""
    ref = _spheroid_reference(csv)
    assert np.ptp(ref["M_RB"], axis=0).max() == 0.0 and np.ptp(ref["M_A"], axis=0).max() == 0.0
    mrb, ma = ref["M_RB"][0], ref["M_A"][0]
    return {
        "rigid_body_mass_matrix": mrb,
        "added_mass_matrix": ma,
        "weight": mrb[0, 0] * _gravity(_value(REMUS, 97)),  # remus100.m 212
        "center_of_gravity": _value(REMUS, 138),
        "center_of_buoyancy": _value(REMUS, 139),
        "time_constants": [_value(REMUS, 193), _value(REMUS, 194), _value(REMUS, 197)],
        "damping_ratios": [_value(REMUS, 195), _value(REMUS, 196)],
        "sway_damping_fade": fade,
    }


def _scaled_damping(fade):
    """Second vehicle (test construction): masses x8, rotational inertia x32,
    weight x8, CG x2, time constants x sqrt(2)."""
    p = _spheroid_damping(fade)
    scale = np.sqrt(np.array([8.0, 8.0, 8.0, 32.0, 32.0, 32.0]))
    return {**p,
            "rigid_body_mass_matrix": np.outer(scale, scale) * p["rigid_body_mass_matrix"],
            "added_mass_matrix": np.outer(scale, scale) * p["added_mass_matrix"],
            "weight": 8.0 * p["weight"],
            "center_of_gravity": 2.0 * p["center_of_gravity"],
            "time_constants": list(np.sqrt(2.0) * np.asarray(p["time_constants"])),
            }


def _cross_flow_spheroid(reynolds, grid):
    geo = _spheroid_geometry()
    return {"length": geo["L"], "beam": geo["D"], "draft": geo["D"],   # remus100.m 219
            "water_density": _value(CFD, 24), "drag_model": "cylinder",
            "cross_flow_reynolds_length": reynolds, "strip_grid": grid}


def _catamaran_draft():
    m, rho, mp = _value(OTTER, 94), _value(OTTER, 91), _value(GEN_CAT, 12)
    nabla = _value(OTTER, 121, {"m": m, "mp": mp, "rho": rho})
    return _value(OTTER, 122, {"nabla": nabla, "Cb_pont": _value(OTTER, 107),
                               "B_pont": _value(OTTER, 104), "L": _value(OTTER, 92)})


def _cross_flow_catamaran(grid):
    return {"length": _value(OTTER, 92), "beam": _value(OTTER, 104),
            "draft": _catamaran_draft(), "water_density": _value(CFD, 24),
            "drag_model": "hoerner", "strip_grid": grid}                # otter.m 245


def _scaled_cross_flow(params):
    return {**params, "length": 2.0 * params["length"], "beam": 2.0 * params["beam"],
            "draft": 2.0 * params["draft"]}


def _lift_drag_spheroid(rho):
    geo = _spheroid_geometry()
    return {"span": geo["D"], "planform_area": geo["S"],                # remus100.m 218
            "parasitic_drag_coefficient": geo["CD_0"],
            "oswald_efficiency": _value(CLD, 54), "water_density": rho}


def _surface_damping_catamaran(gravity, csv=LEGACY_CATAMARAN):
    """Mass and restoring matrices of the catamaran reference; G moved from
    the CO to the CF with H(-r_f), r_f = [LCF 0 0] (otter.m 192-194)."""
    ref = _catamaran_reference(csv)
    m, g_co = ref["M"][0], ref["G"][0]
    h_inv = _hmtrx(-np.array([_value(OTTER, 192), 0.0, 0.0]))
    return {
        "mass_matrix": m,
        "restoring_matrix": h_inv.T @ g_co @ h_inv,
        "max_forward_thrust": _number_in(OTTER, 202, r"Xu = -([0-9.]+) \*") * gravity,
        "max_speed": _value(OTTER, 101),
        "time_constants": [_value(OTTER, 99), _value(OTTER, 100)],
        "damping_ratios": [_number_in(OTTER, n, r"-2 \* ([0-9.]+) \*") for n in (204, 205, 206)],
        "yaw_damping_nonlinearity": _number_in(OTTER, 240, r"\(1 \+ ([0-9.]+) \* abs"),
    }


def _matlab_gravity_catamaran():
    return _gravity(_value(GEN_CAT, 136))  # generator lines 136-137


# --------------------------------------------------------------------------
# MSS transcriptions (each line cited and pinned)
# --------------------------------------------------------------------------
def _mss_dmtrx(p, nu_r):
    """Dmtrx.m 30, 44-49 (submerged branch) + remus100.m 216 (surge fade only)."""
    m = p["rigid_body_mass_matrix"] + p["added_mass_matrix"]
    t1, t2, t6 = p["time_constants"]
    z4, z5 = p["damping_ratios"]
    t3 = t2
    dz = p["center_of_gravity"][2] - p["center_of_buoyancy"][2]
    w4 = np.sqrt(p["weight"] * dz / m[3, 3])
    w5 = np.sqrt(p["weight"] * dz / m[4, 4])
    d = np.diag([m[0, 0] / t1, m[1, 1] / t2, m[2, 2] / t3,
                 m[3, 3] * 2 * z4 * w4, m[4, 4] * 2 * z5 * w5, m[5, 5] / t6])
    d[0, 0] *= np.exp(-3 * np.linalg.norm(nu_r[:3]))  # remus100.m 128, 216
    return d


def _mss_otter_damping(p, nu_r):
    """otter.m 197-207 and 235-240; returns (derivatives, tau_damp)."""
    m, g_cf = p["mass_matrix"], p["restoring_matrix"]
    w3, w4, w5 = (np.sqrt(g_cf[i, i] / m[i, i]) for i in (2, 3, 4))
    z3, z4, z5 = p["damping_ratios"]
    t_sway, t_yaw = p["time_constants"]
    deriv = np.array([-p["max_forward_thrust"] / p["max_speed"], -m[1, 1] / t_sway,
                      -2 * z3 * w3 * m[2, 2], -2 * z4 * w4 * m[3, 3],
                      -2 * z5 * w5 * m[4, 4], -m[5, 5] / t_yaw])
    tau = deriv * nu_r
    tau[5] *= 1 + p["yaw_damping_nonlinearity"] * abs(nu_r[5])
    return deriv, tau


def _interp1_clamped(x, table):
    """cylinderDrag.m 80-86: interp1 inside the table, end values outside."""
    if x < table[0, 0]:
        return table[0, 1]
    if x > table[-1, 0]:
        return table[-1, 1]
    return float(np.interp(x, table[:, 0], table[:, 1]))


def _mss_cylinder_cd(length, beam, nu_r, reynolds_length):
    """cylinderDrag.m 23-107; ``reynolds_length`` replaces the L of line 77."""
    cd_data = _table(CYL, "CD_DATA")
    u_cf = np.sqrt(nu_r[1] ** 2 + nu_r[2] ** 2)            # line 76
    re_number = u_cf * reynolds_length * 1e6               # line 77
    cd = _interp1_clamped(re_number, cd_data)              # lines 80-86
    kappa_table = _table(CYL, "KAPPA_SUBCRITICAL_DATA") if re_number < 2e5 \
        else _table(CYL, "KAPPA_SUPERCRITICAL_DATA")       # lines 89-105
    return cd * _interp1_clamped(length / beam, kappa_table)


def _mss_hoerner(beam, draft):
    """Hoerner.m 25-51."""
    data = _table(HOERNER, "CD_DATA")
    ratio = beam / (2 * draft)
    return float(np.interp(ratio, data[:, 0], data[:, 1])) if ratio <= data[-1, 0] \
        else _number_in(HOERNER, 50, r"= ([0-9.]+);")


def _mss_cross_flow(p, nu_r):
    """crossFlowDrag.m 24-54 (current: 20 midpoints); ``strip_grid="endpoint"``
    gives the pre-2026-08-26 loop ``for xL = -L/2:dx:L/2`` (A-2 section 1)."""
    length, beam, draft = p["length"], p["beam"], p["draft"]
    if p["drag_model"] == "cylinder":
        re_len = beam if p.get("cross_flow_reynolds_length", "diameter") == "diameter" else length
        cd = _mss_cylinder_cd(length, beam, nu_r, re_len)
    else:
        cd = _mss_hoerner(beam, draft)
    n_strips = int(_value(CFD, 25))
    dx = length / n_strips
    if p.get("strip_grid", "midpoint") == "midpoint":
        xs = -length / 2 + (np.arange(1, n_strips + 1) - 0.5) * dx   # line 41
    else:
        xs = -length / 2 + dx * np.arange(n_strips + 1)
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    u_h = np.abs(v_r + xs * r) * (v_r + xs * r)                      # line 46
    u_v = np.abs(w_r + xs * q) * (w_r + xs * q)                      # line 47
    k = -0.5 * p["water_density"] * draft * cd * dx                 # lines 48-51
    return np.array([0.0, k * u_h.sum(), k * u_v.sum(), 0.0, k * (xs * u_v).sum(), k * (xs * u_h).sum()])


def _mss_lift_drag(p, nu_r):
    """forceLiftDrag.m 26-40 + coeffLiftDrag.m 54-68 (sigma = 0), remus100.m 127-128."""
    alpha = np.arctan2(nu_r[2], nu_r[0])
    u_r = np.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2)
    b, s, e = p["span"], p["planform_area"], p["oswald_efficiency"]
    ar = b ** 2 / s
    cl = np.pi * ar / (1 + np.sqrt(1 + (ar / 2) ** 2)) * alpha
    cd = p["parasitic_drag_coefficient"] + cl ** 2 / (np.pi * e * ar)
    f_drag = 0.5 * p["water_density"] * u_r ** 2 * s * cd
    f_lift = 0.5 * p["water_density"] * u_r ** 2 * s * cl
    return np.array([np.cos(alpha) * (-f_drag) - np.sin(alpha) * (-f_lift), 0.0,
                     np.sin(alpha) * (-f_drag) + np.cos(alpha) * (-f_lift), 0.0, 0.0, 0.0])


# --------------------------------------------------------------------------
# Builders and numpy source
# --------------------------------------------------------------------------
def _build_submerged(p):
    block = _contract(LINEAR_DAMPING)
    c = block.preprocess_submerged_linear_damping(**p)
    return c, block.submerged_linear_damping_casadi(c)


def _build_surface(p):
    block = _contract(LINEAR_DAMPING)
    c = block.preprocess_surface_linear_damping(**p)
    return c, block.surface_linear_damping_casadi(c)


def _build_cross_flow(p):
    block = _contract(CROSS_FLOW)
    return block.cross_flow_drag_casadi(block.preprocess_cross_flow_drag(**p))


def _build_lift_drag(p):
    block = _contract(LIFT_DRAG)
    return block.lift_drag_casadi(block.preprocess_lift_drag(**p))


def _source():
    try:
        from more_generic_models.dynamics.plant.asv_catamaran.asv_catamaran import ASVCatamaran
        from more_generic_models.dynamics.plant.auv_spheroid.auv_spheroid import AUVSpheroid
        from more_generic_models.dynamics.plant.matrices.hydrodynamics import HydroForces
        from more_generic_models.dynamics.plant.matrices.linear_damping import DampingMatrix
    except ImportError as exc:
        pytest.skip(
            "numpy source not importable; install it into the venv with "
            f"pip install -e <more>/more_generic_models (job A-4c) — {exc}"
        )
    return {"spheroid": AUVSpheroid, "catamaran": ASVCatamaran,
            "hydro": HydroForces, "damping": DampingMatrix}


@pytest.fixture
def full_table_source(monkeypatch):
    """The numpy source with its rounded Re column replaced, in memory only, by
    the full-precision ``cylinderDrag.m`` table (A-2 section 6, edit 1)."""
    src = _source()
    monkeypatch.setattr(src["hydro"], "_CYLINDER_CD_DATA", _table(CYL, "CD_DATA"))
    return src


# --------------------------------------------------------------------------
# Reference sanity and pinned lines (no port needed)
# --------------------------------------------------------------------------
def test_reference_shapes():
    sph, cat = _spheroid_reference(), _catamaran_reference()
    for key in ("nu_r", "tau_lift_drag", "tau_crossflow"):
        assert sph[key].shape == (50, 6), key
    for key in ("M_RB", "M_A", "D"):
        assert sph[key].shape == (50, 6, 6), key
    for key in ("nu_r", "tau_damp", "tau_crossflow"):
        assert cat[key].shape == (50, 6), key
    for key in ("M", "G", "D"):
        assert cat[key].shape == (50, 6, 6), key
    # U_r and alpha of the reference are remus100.m 127-128 of its nu_r
    np.testing.assert_allclose(sph["U_r"], np.linalg.norm(sph["nu_r"][:, :3], axis=1), atol=1e-12)
    np.testing.assert_allclose(sph["alpha"], np.arctan2(sph["nu_r"][:, 2], sph["nu_r"][:, 0]), atol=1e-12)


def test_cited_lines_are_unchanged():
    for (rel, number), text in CITED_LINES.items():
        line = _line(rel, number)
        assert line.startswith(text), (rel, number, line)
    assert _table(CYL, "CD_DATA").shape == (28, 2)
    assert _table(CYL, "KAPPA_SUBCRITICAL_DATA").shape == (7, 2)
    assert _table(CYL, "KAPPA_SUPERCRITICAL_DATA").shape == (7, 2)
    assert _table(HOERNER, "CD_DATA").shape == (20, 2)


def test_transcriptions_reproduce_the_stored_matlab_reference():
    """The test-side MSS transcriptions, run on the legacy settings, land on the
    stored MATLAB numbers; so they are fit to serve as default-path references."""
    sph, cat = _spheroid_reference(), _catamaran_reference()
    p_d = _spheroid_damping(True)
    p_cf = _cross_flow_spheroid("length", "endpoint")
    p_ld = _lift_drag_spheroid(_value(GEN_SPH, 467))
    p_cat = _surface_damping_catamaran(_matlab_gravity_catamaran())
    p_cat_cf = _cross_flow_catamaran("endpoint")
    for k in range(50):
        d = _mss_dmtrx(p_d, sph["nu_r"][k])
        d[1, 1] *= np.exp(-3 * sph["U_r"][k])   # generator line 268 (sway fade)
        assert _max_diff(d, sph["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cf, sph["nu_r"][k]), sph["tau_crossflow"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_lift_drag(p_ld, sph["nu_r"][k]), sph["tau_lift_drag"][k]) <= G1_TOLERANCE, k
        deriv, tau = _mss_otter_damping(p_cat, cat["nu_r"][k])
        assert _max_diff(np.diag(deriv), cat["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(tau, cat["tau_damp"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cat_cf, cat["nu_r"][k]), cat["tau_crossflow"][k]) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G1 — numpy source vs the stored MATLAB reference
# --------------------------------------------------------------------------
def test_G1_numpy_source_spheroid_damping_and_lift_drag():
    vehicle = _source()["spheroid"]()
    ref = _spheroid_reference()
    for k in range(50):
        assert _max_diff(vehicle.get_D(ref["U_r"][k]), ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        tau = vehicle.get_tau_lift_drag(ref["alpha"][k], ref["U_r"][k])
        assert _max_diff(tau, ref["tau_lift_drag"][k]) <= G1_TOLERANCE, ("lift_drag", k)


def test_G1_numpy_source_spheroid_cross_flow_with_full_table(full_table_source):
    vehicle = full_table_source["spheroid"]()
    ref = _spheroid_reference()
    for k in range(50):
        tau = vehicle.get_tau_cross_flow(ref["nu_r"][k], "cylinder")
        assert _max_diff(tau, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_G1_numpy_source_rounded_table_residual_is_pinned():
    """A-2: the unpatched source misses MATLAB by the rounded Re column only."""
    vehicle = _source()["spheroid"]()
    ref = _spheroid_reference()
    worst = max(_max_diff(vehicle.get_tau_cross_flow(ref["nu_r"][k], "cylinder"),
                          ref["tau_crossflow"][k]) for k in range(50))
    assert G1_TOLERANCE < worst <= ROUNDED_TABLE_RESIDUAL_BOUND, worst


def test_numpy_source_drag_models_copy_equals_hydro_forces_copy():
    """``drag_models.py::DragModels`` (map rows -> cross_flow.py) carries a second
    copy of both coefficient tables; it agrees with the ``HydroForces`` copy the
    vehicles call, so the G1/G2 tests on ``HydroForces`` cover both rows."""
    src = _source()
    from more_generic_models.dynamics.plant.matrices.drag_models import DragModels
    rng = np.random.default_rng(SEED)
    for b, t in rng.uniform(0.05, 2.0, size=(200, 2)):
        assert abs(DragModels.hoerner(b, t) - src["hydro"]._hoerner(b, t)) <= G2_TOLERANCE, (b, t)
    for nu_r in _random_nu_r():
        length, beam = rng.uniform(0.5, 6.0), rng.uniform(0.05, 0.5)
        assert abs(DragModels.cylinder_drag(length, beam, nu_r)
                   - src["hydro"]._cylinder_drag(length, beam, nu_r)) <= G2_TOLERANCE


def test_G1_numpy_source_catamaran_damping_and_cross_flow():
    vessel = _source()["catamaran"]()
    ref = _catamaran_reference()
    for k in range(50):
        assert _max_diff(vessel.get_D(), ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(vessel.get_tau_damping(ref["nu_r"][k]), ref["tau_damp"][k]) <= G1_TOLERANCE, k
        assert _max_diff(vessel.get_tau_cross_flow(ref["nu_r"][k]), ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G1 — block (legacy settings) vs the stored MATLAB reference
# --------------------------------------------------------------------------
def test_G1_block_spheroid_damping_with_sway_fade_flag():
    constants, function = _build_submerged(_spheroid_damping(True))
    ref = _spheroid_reference()
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], ref["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(out["tau"].reshape(-1), -ref["D"][k] @ ref["nu_r"][k]) <= G1_TOLERANCE, k


def test_G1_block_spheroid_cross_flow_legacy_settings():
    function = _build_cross_flow(_cross_flow_spheroid("length", "endpoint"))
    ref = _spheroid_reference()
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_G1_block_spheroid_lift_drag_template_density():
    function = _build_lift_drag(_lift_drag_spheroid(_value(GEN_SPH, 467)))
    ref = _spheroid_reference()
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_lift_drag"][k]) <= G1_TOLERANCE, k


def test_G1_block_catamaran_surface_damping_and_cross_flow():
    _, damping = _build_surface(_surface_damping_catamaran(_matlab_gravity_catamaran()))
    cross_flow = _build_cross_flow(_cross_flow_catamaran("endpoint"))
    ref = _catamaran_reference()
    for k in range(50):
        out = _evaluate(damping, ref["nu_r"][k])
        assert _max_diff(out["D"], -ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), ref["tau_damp"][k]) <= G1_TOLERANCE, ("tau", k)
        tau_cf = _evaluate(cross_flow, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(tau_cf, ref["tau_crossflow"][k]) <= G1_TOLERANCE, ("cross_flow", k)


# --------------------------------------------------------------------------
# G2 — CasADi blocks vs numpy source (reference cases + 1000 seeded states)
# --------------------------------------------------------------------------
def _g2_states(reference):
    return np.vstack([reference["nu_r"], _random_nu_r()])


def test_G2_function_signatures():
    for function, outs in (
        (_build_submerged(_spheroid_damping(False))[1], ["D", "tau"]),
        (_build_surface(_surface_damping_catamaran(_matlab_gravity_catamaran()))[1], ["D", "tau"]),
        (_build_cross_flow(_cross_flow_spheroid("diameter", "midpoint")), ["tau"]),
        (_build_lift_drag(_lift_drag_spheroid(_value(FLD, 26))), ["tau"]),
    ):
        assert function.name_in() == ["nu_r"]
        assert function.size_in(0) == (6, 1)
        assert function.name_out() == outs
        sizes = {"D": (6, 6), "tau": (6, 1)}
        for i, name in enumerate(outs):
            assert function.size_out(i) == sizes[name], name


@pytest.mark.parametrize("scaled", [False, True], ids=["remus_like", "scaled"])
@pytest.mark.parametrize("fade", [False, True], ids=["default", "sway_fade"])
def test_G2_submerged_damping_matches_source(fade, scaled):
    damping = _source()["damping"]
    p = _scaled_damping(fade) if scaled else _spheroid_damping(fade)
    constants, function = _build_submerged(p)
    coeff = damping.get_damping_coeff(p["weight"], p["center_of_gravity"], p["center_of_buoyancy"],
                                      p["rigid_body_mass_matrix"], p["added_mass_matrix"],
                                      p["time_constants"], p["damping_ratios"])
    assert _max_diff(constants.damping_coefficients, coeff) <= G2_TOLERANCE
    for k, nu_r in enumerate(_g2_states(_spheroid_reference())):
        expected = damping.D_linear(coeff, np.linalg.norm(nu_r[:3]))  # fades sway too
        if not fade:
            expected[1, 1] = coeff[1]   # default: surge-only fade, remus100.m 216
        out = _evaluate(function, nu_r)
        assert _max_diff(out["D"], expected) <= G2_TOLERANCE, k
        assert _max_diff(out["tau"].reshape(-1), -expected @ nu_r) <= G2_TOLERANCE, k


@pytest.mark.parametrize("scale", [1.0, 3.0], ids=["otter_like", "scaled"])
def test_G2_surface_damping_matches_source(scale):
    damping = _source()["damping"]
    g = _matlab_gravity_catamaran()
    p = _surface_damping_catamaran(g)
    p = {**p, "mass_matrix": scale ** 3 * p["mass_matrix"],
         "restoring_matrix": scale ** 2 * p["restoring_matrix"],
         "max_forward_thrust": scale ** 3 * p["max_forward_thrust"]}
    constants, function = _build_surface(p)
    # source D_linear_catamaran(G_CF, MRB, MA, Umax, g, zeta_345, T_126, coeff):
    # Xu = coeff[0] g / Umax, the rest -1 x the otter.m forms (asv_catamaran_params.py 25-30)
    coeff = np.array([-p["max_forward_thrust"] / g, -1.0, -1.0, -1.0, -1.0, -1.0])
    t_sway, t_yaw = p["time_constants"]
    d_src = damping.D_linear_catamaran(p["restoring_matrix"], p["mass_matrix"], np.zeros((6, 6)),
                                       p["max_speed"], g, p["damping_ratios"],
                                       [np.nan, t_sway, t_yaw], coeff)
    assert _max_diff(constants.damping_derivatives, np.diag(d_src)) <= G2_TOLERANCE
    for k, nu_r in enumerate(_g2_states(_catamaran_reference())):
        out = _evaluate(function, nu_r)
        assert _max_diff(out["D"], -d_src) <= G2_TOLERANCE, k
        # get_tau_damping hard-codes the yaw factor 10 (linear_damping.py 221)
        assert _max_diff(out["tau"].reshape(-1), damping.get_tau_damping(d_src, nu_r)) <= G2_TOLERANCE, k


@pytest.mark.parametrize("vehicle", ["spheroid", "spheroid_scaled", "catamaran", "catamaran_scaled"])
def test_G2_cross_flow_legacy_settings_match_source(vehicle, full_table_source):
    hydro = full_table_source["hydro"]
    p = (_cross_flow_spheroid("length", "endpoint") if vehicle.startswith("spheroid")
         else _cross_flow_catamaran("endpoint"))
    if vehicle.endswith("scaled"):
        p = _scaled_cross_flow(p)
    function = _build_cross_flow(p)
    ref = _spheroid_reference() if vehicle.startswith("spheroid") else _catamaran_reference()
    for k, nu_r in enumerate(_g2_states(ref)):
        expected = hydro.cross_flow_drag(p["length"], p["beam"], p["draft"], nu_r,
                                         p["drag_model"], p["water_density"])
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, expected) <= G2_TOLERANCE, (vehicle, k)


@pytest.mark.parametrize("scale", [1.0, 2.0], ids=["remus_like", "scaled"])
def test_G2_lift_drag_matches_source(scale):
    hydro = _source()["hydro"]
    p = _lift_drag_spheroid(_value(GEN_SPH, 467))
    p = {**p, "span": scale * p["span"], "planform_area": scale ** 2 * p["planform_area"]}
    function = _build_lift_drag(p)
    for k, nu_r in enumerate(_g2_states(_spheroid_reference())):
        alpha, u_r = np.arctan2(nu_r[2], nu_r[0]), np.linalg.norm(nu_r[:3])
        expected = hydro.force_lift_drag(p["span"], p["planform_area"], p["parasitic_drag_coefficient"],
                                         alpha, u_r, p["water_density"], p["oswald_efficiency"])
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, expected) <= G2_TOLERANCE, k


# --------------------------------------------------------------------------
# G4 — a broken model must fail (control first, then one perturbation each)
# --------------------------------------------------------------------------
def _worst(function, cases, expected, key):
    return max(_max_diff(_evaluate(function, nu_r)[key].reshape(np.shape(e)), e)
               for nu_r, e in zip(cases, expected))


def test_G4_submerged_damping_roll_ratio_times_zero_is_detected():
    ref = _spheroid_reference()
    p = _spheroid_damping(True)
    assert _worst(_build_submerged(p)[1], ref["nu_r"], ref["D"], "D") <= G1_TOLERANCE
    broken = {**p, "damping_ratios": [0.0, p["damping_ratios"][1]]}
    assert _worst(_build_submerged(broken)[1], ref["nu_r"], ref["D"], "D") > G4_FACTOR * G1_TOLERANCE


def test_G4_surface_damping_heave_ratio_times_zero_is_detected():
    ref = _catamaran_reference()
    p = _surface_damping_catamaran(_matlab_gravity_catamaran())
    assert _worst(_build_surface(p)[1], ref["nu_r"], ref["tau_damp"], "tau") <= G1_TOLERANCE
    broken = {**p, "damping_ratios": [0.0] + list(p["damping_ratios"][1:])}
    assert _worst(_build_surface(broken)[1], ref["nu_r"], ref["tau_damp"], "tau") > G4_FACTOR * G1_TOLERANCE


def test_G4_cross_flow_draft_plus_1_percent_is_detected():
    ref = _spheroid_reference()
    p = _cross_flow_spheroid("length", "endpoint")
    assert _worst(_build_cross_flow(p), ref["nu_r"], ref["tau_crossflow"], "tau") <= G1_TOLERANCE
    broken = {**p, "draft": 1.01 * p["draft"]}
    assert _worst(_build_cross_flow(broken), ref["nu_r"], ref["tau_crossflow"], "tau") > G4_FACTOR * G1_TOLERANCE


def test_G4_lift_drag_area_plus_1_percent_is_detected():
    ref = _spheroid_reference()
    p = _lift_drag_spheroid(_value(GEN_SPH, 467))
    assert _worst(_build_lift_drag(p), ref["nu_r"], ref["tau_lift_drag"], "tau") <= G1_TOLERANCE
    broken = {**p, "planform_area": 1.01 * p["planform_area"]}
    assert _worst(_build_lift_drag(broken), ref["nu_r"], ref["tau_lift_drag"], "tau") > G4_FACTOR * G1_TOLERANCE


# --------------------------------------------------------------------------
# Headline (E-18, E-20 Q2 a, E-10 a, E-20 Q3 b): default path vs current MSS
# --------------------------------------------------------------------------
def _states():
    return np.vstack([_spheroid_reference()["nu_r"], _random_nu_r()])


def test_MSS_submerged_damping_default_fades_surge_only():
    """Default = Dmtrx.m + remus100.m:216; the template's sway fade is the flag."""
    p = _spheroid_damping(False)
    _, function = _build_submerged(p)
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)
        assert _max_diff(out["D"], _mss_dmtrx(p, nu_r)) <= G1_TOLERANCE, ("Dmtrx.m 30-49, remus100.m 216", k)


def test_MSS_surface_damping_equals_otter():
    """otter.m 197-207, 235-240 with otter.m's own g = 9.81 (line 90)."""
    p = _surface_damping_catamaran(_value(OTTER, 90))
    constants, function = _build_surface(p)
    for k, nu_r in enumerate(np.vstack([_catamaran_reference()["nu_r"], _random_nu_r()])):
        deriv, tau = _mss_otter_damping(p, nu_r)
        out = _evaluate(function, nu_r)
        assert _max_diff(out["tau"].reshape(-1), tau) <= G1_TOLERANCE, ("otter.m 235-240", k)
    assert _max_diff(constants.damping_derivatives, deriv) <= G1_TOLERANCE


@pytest.mark.parametrize("model", ["cylinder", "hoerner"])
def test_MSS_cross_flow_length_flag_equals_current_mss(model):
    """``cross_flow_reynolds_length="length"`` + default grid = crossFlowDrag.m
    (6a2a064) and cylinderDrag.m line 77 exactly (E-20 Q3 b, comparison path)."""
    p = (_cross_flow_spheroid("length", "midpoint") if model == "cylinder"
         else _cross_flow_catamaran("midpoint"))
    if model == "hoerner":
        p = {**p, "cross_flow_reynolds_length": "length"}
    function = _build_cross_flow(p)
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, _mss_cross_flow(p, nu_r)) <= G1_TOLERANCE, ("crossFlowDrag.m 24-54", k)


def test_E20_Q3b_default_reynolds_number_uses_the_diameter():
    """Default flags (diameter, midpoint) = crossFlowDrag.m with line 77's L
    replaced by B; and the two settings really differ in the drag crisis."""
    block = _contract(CROSS_FLOW)
    signature = inspect.signature(block.preprocess_cross_flow_drag)
    assert signature.parameters["cross_flow_reynolds_length"].default == "diameter"
    assert signature.parameters["strip_grid"].default == "midpoint"
    p = _cross_flow_spheroid("diameter", "midpoint")
    defaults = {k: v for k, v in p.items() if k not in ("cross_flow_reynolds_length", "strip_grid")}
    function = _build_cross_flow(defaults)
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, _mss_cross_flow(p, nu_r)) <= G1_TOLERANCE, k
    sway = np.array([0.0, 0.3, 0.0, 0.0, 0.0, 0.0])   # A-18's example speed
    on_length = _mss_cross_flow({**p, "cross_flow_reynolds_length": "length"}, sway)
    assert _max_diff(_mss_cross_flow(p, sway), on_length) > G4_FACTOR * G1_TOLERANCE


def test_E10_endpoint_grid_flag_reproduces_the_template():
    """strip_grid="endpoint" equals the pre-2026-08-26 grid (the stored
    reference); the default grid does not."""
    ref = _spheroid_reference()
    legacy = _build_cross_flow(_cross_flow_spheroid("length", "endpoint"))
    current = _build_cross_flow(_cross_flow_spheroid("length", "midpoint"))
    worst_legacy = max(_max_diff(_evaluate(legacy, n)["tau"].reshape(-1), ref["tau_crossflow"][k])
                       for k, n in enumerate(ref["nu_r"]))
    worst_current = max(_max_diff(_evaluate(current, n)["tau"].reshape(-1), ref["tau_crossflow"][k])
                        for k, n in enumerate(ref["nu_r"]))
    assert worst_legacy <= G1_TOLERANCE
    assert worst_current > G4_FACTOR * G1_TOLERANCE


def test_MSS_lift_drag_equals_force_lift_drag():
    """forceLiftDrag.m rho = 1026 (line 26), e = 0.3 (coeffLiftDrag.m 54)."""
    p = _lift_drag_spheroid(_value(FLD, 26))
    function = _build_lift_drag(p)
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, _mss_lift_drag(p, nu_r)) <= G1_TOLERANCE, ("forceLiftDrag.m 26-40", k)


# --------------------------------------------------------------------------
# G1-MSS (job A-30): default paths vs MATLAB running current MSS (job A-26)
# --------------------------------------------------------------------------
def test_transcriptions_reproduce_the_current_mss_matlab_reference():
    """The second check: the test-side transcriptions, on current-MSS settings,
    land on the A-26 files (A-26 ledger step 4 c, now a test)."""
    sph, cat = _spheroid_reference(MSS_SPHEROID), _catamaran_reference(MSS_CATAMARAN)
    p_d = _spheroid_damping(False, MSS_SPHEROID)
    p_cf = _cross_flow_spheroid("length", "midpoint")
    p_ld = _lift_drag_spheroid(_value(FLD, 26))
    p_cat = _surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN)
    p_cat_cf = _cross_flow_catamaran("midpoint")
    for k in range(50):
        assert _max_diff(_mss_dmtrx(p_d, sph["nu_r"][k]), sph["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cf, sph["nu_r"][k]), sph["tau_crossflow"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_lift_drag(p_ld, sph["nu_r"][k]), sph["tau_lift_drag"][k]) <= G1_TOLERANCE, k
        deriv, tau = _mss_otter_damping(p_cat, cat["nu_r"][k])
        assert _max_diff(np.diag(deriv), cat["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(tau, cat["tau_damp"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cat_cf, cat["nu_r"][k]), cat["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_G1_MSS_submerged_damping_default_equals_matlab():
    """Default (surge fade only, remus100.m 216) with M_RB, M_A of the A-26 file
    (imlay61.m rho = 1026) vs MATLAB's D and -D nu_r."""
    ref = _spheroid_reference(MSS_SPHEROID)
    p = _spheroid_damping(False, MSS_SPHEROID)
    p = {k: v for k, v in p.items() if k != "sway_damping_fade"}   # the default
    _, function = _build_submerged(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), -ref["D"][k] @ ref["nu_r"][k]) <= G1_TOLERANCE, ("tau", k)


def test_G1_MSS_cross_flow_length_flag_equals_matlab():
    """D-MSS-1: MSS's behaviour (Re on the length, cylinderDrag.m:77) is the
    flag ``cross_flow_reynolds_length="length"``; with the default grid it
    equals MATLAB's ``crossFlowDrag.m`` (6a2a064)."""
    ref = _spheroid_reference(MSS_SPHEROID)
    p = {k: v for k, v in _cross_flow_spheroid("length", "midpoint").items() if k != "strip_grid"}
    function = _build_cross_flow(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_D_MSS_1_default_cross_flow_equals_transcription_and_register_effect():
    """D-MSS-1 (30_checks/README.md): the default (Re on the diameter) equals
    the transcription with line 77's L replaced by D, departs from MATLAB on
    the reference cases, and has the register's stated effect: REMUS
    Y_v|v| = -125 vs -66 N/(m/s)^2 at 0.3 m/s sideways (diameter vs length;
    the two numbers are typed from the register row, whole-number precision)."""
    ref = _spheroid_reference(MSS_SPHEROID)
    p = _cross_flow_spheroid("diameter", "midpoint")
    default = _build_cross_flow({k: v for k, v in p.items()
                                 if k not in ("cross_flow_reynolds_length", "strip_grid")})
    worst_vs_matlab = 0.0
    for k in range(50):
        out = _evaluate(default, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, _mss_cross_flow(p, ref["nu_r"][k])) <= G1_TOLERANCE, k
        worst_vs_matlab = max(worst_vs_matlab, _max_diff(out, ref["tau_crossflow"][k]))
    assert worst_vs_matlab > G4_FACTOR * G1_TOLERANCE   # the deviation is live
    v = 0.3
    sway = np.array([0.0, v, 0.0, 0.0, 0.0, 0.0])
    length_flag = _build_cross_flow({**p, "cross_flow_reynolds_length": "length"})
    y_v_diameter = _evaluate(default, sway)["tau"].reshape(-1)[1] / (v * abs(v))
    y_v_length = _evaluate(length_flag, sway)["tau"].reshape(-1)[1] / (v * abs(v))
    assert round(y_v_diameter) == -125, y_v_diameter
    assert round(y_v_length) == -66, y_v_length


def test_G1_MSS_lift_drag_equals_matlab():
    """forceLiftDrag.m rho = 1026 (line 26), e = 0.3 (coeffLiftDrag.m 54)."""
    ref = _spheroid_reference(MSS_SPHEROID)
    function = _build_lift_drag(_lift_drag_spheroid(_value(FLD, 26)))
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_lift_drag"][k]) <= G1_TOLERANCE, k


def test_G1_MSS_catamaran_surface_damping_equals_matlab():
    """otter.m g = 9.81 (line 90); M and G of the A-26 file (inertia about the
    combined CG)."""
    ref = _catamaran_reference(MSS_CATAMARAN)
    _, function = _build_surface(_surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN))
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], -ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), ref["tau_damp"][k]) <= G1_TOLERANCE, ("tau", k)


def test_G1_MSS_catamaran_cross_flow_default_equals_matlab():
    """Hoerner on the default grid (20 midpoints): no registered deviation."""
    ref = _catamaran_reference(MSS_CATAMARAN)
    p = {k: v for k, v in _cross_flow_catamaran("midpoint").items() if k != "strip_grid"}
    function = _build_cross_flow(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_D_MSS_2_hoerner_below_its_table_raises_naming_the_ratio():
    """D-MSS-2 (owner E-24): B/(2T) below Hoerner.m's first data point (MSS:
    interp1 -> NaN, lines 47-48) raises ValueError naming the ratio; the first
    data point itself is accepted."""
    block = _contract(CROSS_FLOW)
    first = _table(HOERNER, "CD_DATA")[0, 0]
    draft = 1.0
    beam = first * draft          # ratio = first / 2, below the table
    ratio = beam / (2 * draft)
    with pytest.raises(ValueError) as error:
        block.preprocess_cross_flow_drag(length=2.0, beam=beam, draft=draft,
                                         water_density=_value(CFD, 24), drag_model="hoerner")
    assert f"{ratio}" in str(error.value), str(error.value)
    block.preprocess_cross_flow_drag(length=2.0, beam=2 * first * draft, draft=draft,
                                     water_density=_value(CFD, 24), drag_model="hoerner")


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


GENERIC = {
    LINEAR_DAMPING: {"preprocess_submerged_linear_damping": {"sway_damping_fade"},
                     "preprocess_surface_linear_damping": set()},
    CROSS_FLOW: {"preprocess_cross_flow_drag": {"cross_flow_reynolds_length", "strip_grid"}},
    LIFT_DRAG: {"preprocess_lift_drag": set()},
}


@pytest.mark.parametrize("module_name", sorted(GENERIC))
def test_generic_block_reads_no_vehicle_name_and_has_no_vehicle_defaults(module_name):
    block = _contract(module_name)
    code = _code_without_docstrings(block)
    assert not [n for n in VEHICLE_NAMES if n in code]
    for function_name, flags in GENERIC[module_name].items():
        signature = inspect.signature(getattr(block, function_name))
        defaults = {n for n, p in signature.parameters.items() if p.default is not inspect.Parameter.empty}
        assert defaults == flags, (function_name, defaults)
