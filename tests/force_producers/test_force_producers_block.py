"""Gate tests for the force-producer blocks.

Written before the blocks exist (job A-22, 2026-10-05); the blocks are ported
by job A-25. Every block test fails today with "block not ported yet". One
source test fails on purpose: the outboard motor's reverse thrust (owner
decision E-20 Q5 a; see ``test_E20_Q5a_numpy_source_full_reverse_pushes_backwards``).

Changed by job A-28 (2026-10-05), owner decision E-22 Q1 a and E-23: the
Wageningen default follows MSS ``wageningen.m`` (no ``|J|``, no clamp to
[0, 1.3]); the numpy source's clamp stays reachable behind the flag
``clip_advance_ratio=True``. The tests that compare a producer with the numpy
source now run on that flag; the default is checked against MSS. Until job
A-29 ports the flag, the tests that pass it and the default-path tests fail.

Rewired by job A-30 (2026-10-05): the differential thruster's default (otter.m
``g = 9.81``) is compared with MATLAB running current MSS (job A-26,
``differential_thruster_mss_current.csv``); the legacy CSV (latitude gravity)
keeps the numpy source's G1 and the block's G1/G4 on the source's settings.
Fins and the propeller have no current-MSS MATLAB file (A-26 scope).

Producers and their references
------------------------------
=====================  ==========================================  =======================================
block                  numpy source                                reference
=====================  ==========================================  =======================================
differential thruster  ``thruster_differential.py``                stored MATLAB CSV + ``otter.m``
propeller (B-series)   ``thruster_wagenigen.py``                   stored MATLAB CSV + ``remus100.m``
K_T/K_Q polynomial     ``propeller_models/wagenigen/wagenigen.py`` MSS ``WageningData.mat`` / ``.txt``
fins (rudder, sterns)  ``fins_auv_physical.py``                    ``remus100.m`` (transcribed)
VSIM fins              ``fins_auv_actuation_vsim.py``              numpy source only (G2)
outboard motor         ``electrical_outboard_motor.py``            numpy source, reverse sign fixed (G2)
=====================  ==========================================  =======================================

Contract the porter must provide (``more_dynamics.models.force_producers``)
---------------------------------------------------------------------------
Every producer: ``preprocess_<name>(**params) -> constants`` (numpy) and
``<name>_casadi(constants) -> ca.Function`` with the command input(s) first,
then ``"nu_r"`` (6x1), and the single output ``["tau"]`` (6x1, N and N m,
BODY, about the CO). Keywords have no default except the named flags.

* ``differential_thruster``: ``preprocess_differential_thruster(
  positive_thrust_coefficients, negative_thrust_coefficients,
  thruster_positions, thruster_directions, max_forward_thrust,
  max_reverse_thrust)``; per-thruster arrays ``[left, right]``: coefficients
  (2,) in N s^2/rad^2, positions (2,3) m, directions (2,3) (normalised inside),
  ``max_forward_thrust`` / ``max_reverse_thrust`` (2,) N, positive magnitudes,
  giving ``n_max = sqrt(max_forward_thrust / k_pos)`` and ``n_min =
  -sqrt(max_reverse_thrust / k_neg)`` (``otter.m`` 136-137). Constants expose
  ``max_speed``, ``min_speed`` (2,) rad/s and ``allocation_matrix`` (6,2).
  Function inputs ``["n", "nu_r"]``, ``n`` (2x1) shaft speed in rad/s,
  saturated inside; thrust ``k_pos n|n|`` for ``n > 0``, else ``k_neg n|n|``.
* ``propeller``: ``preprocess_propeller(diameter, max_speed, thrust_deduction,
  wake_fraction, pitch_diameter_ratio, blade_area_ratio, blade_count,
  max_advance_number, roll_moment_scale, water_density, position, orientation,
  open_water_model, thrust_torque_coefficients=None,
  clip_advance_ratio=False)``; ``max_speed`` rpm;
  ``orientation`` (roll, pitch, yaw) rad of the shaft; ``open_water_model`` in
  {"linearized", "full", "bollard"} (no default);
  ``thrust_torque_coefficients`` = optional exact ``(KT_0, KQ_0, KT_max,
  KQ_max)`` overriding the polynomial values (E-16 parameter style; MSS prints
  them in ``remus100.m`` 158-162). Constants expose
  ``thrust_torque_coefficients`` (4,), the values in use. Inputs ``["n",
  "nu_r"]``, ``n`` (1x1) rpm; advance speed ``(1 - w) |nu_r[0:3]|``. The
  "full" model keeps the numpy source's own clip of ``J`` to
  ``[0, 2 max_advance_number]`` (``thruster_wagenigen.py`` 113) on both
  paths; ``clip_advance_ratio`` is passed to the Wageningen polynomial.
* ``wageningen_kt_kq``: ``preprocess_wageningen(pitch_diameter_ratio,
  blade_area_ratio, blade_count, clip_advance_ratio=False)``;
  ``wageningen_casadi(c)`` input ``["J"]`` (1x1), outputs ``["KT", "KQ"]``
  (1x1 each) — the one exception to the ``tau`` output. **Default (A-28,
  E-22 Q1 a):** MSS ``wageningen.m``, the polynomial at any ``J``, including
  ``J < 0`` and ``J > 1.3``. ``clip_advance_ratio=True``: the numpy source,
  ``J`` replaced by ``min(|J|, 1.3)`` (``wagenigen.py`` 11).
* ``fins``: ``preprocess_fins(rudder_area, stern_plane_area,
  rudder_lift_coefficient, stern_plane_lift_coefficient, rudder_position,
  stern_plane_position, max_deflection, water_density,
  convention="starboard_down_positive")``; areas m^2 (MSS: ``2 * S_fin``),
  positions = x of the surface in m, ``max_deflection`` rad. Inputs
  ``["delta", "nu_r"]``, ``delta`` = (rudder, stern plane) rad, saturated
  inside. The default convention equals ``remus100.m`` 236-252.
* ``vsim_fins``: ``preprocess_vsim_fins(max_forces, positions,
  max_deflection)``; (N,3) arrays. Inputs ``["delta", "nu_r"]``, ``delta``
  (Nx1).
* ``outboard_motor``: two functions.
  ``preprocess_outboard_motor_rpm(max_thrust, max_power, efficiency, position,
  max_speed, propeller_diameter, pitch_diameter_ratio, blade_area_ratio,
  blade_count, thrust_deduction, water_density, advance_speed_factor,
  propwash_factor, reverse_thrust_factor, reverse_torque_factor,
  power_limit=True, clip_advance_ratio=False)`` (the flag passed to the
  Wageningen polynomial; signed ``J``, so backing down with forward flow
  gives ``J < 0`` on the default path) and ``outboard_motor_rpm_casadi(c)`` with inputs
  ``["n", "delta", "nu_r"]``: ``n`` the actual propeller speed (rpm,
  saturated inside), ``delta`` the actual steering angle (rad, port-positive,
  the source's internal convention). Rpm lag, steering rate limit and battery
  are states and are not part of this force map. **Reverse fixed (E-20 Q5 a):**
  for ``n < 0`` the reverse factors scale ``K_T`` and ``K_Q`` without the
  extra ``sign(n)``, so full reverse gives a negative surge force.
  ``preprocess_outboard_motor_throttle(max_thrust, max_propulsive_power,
  position)`` and ``outboard_motor_throttle_casadi(c)`` with inputs
  ``["thrust", "delta", "nu_r"]`` (actual thrust state in N).

Map rows covered, numpy source read in full (``more_generic_models/more_generic_models/``)
---------------------------------------------------------------------------------------
* ``dynamics/propulsion/thruster/thruster_base.py`` 8-99 (mode, saturation,
  ``rpm_to_rps``): exercised through every thruster.
* ``.../thruster_differential/thruster_differential.py`` 7-326 (state
  ``update_state`` 136-143 is a lag, not tested here; ``thrust_from_speed``
  156-171, ``tau_from_thrust`` 194-196, ``_compute_B`` 272-279, limits 84-94).
* ``.../thruster_wagenigen/thruster_wagenigen.py`` 12-268 (``_compute_kt_kq``
  104-115, linearised 117-137, bollard 139-142, ``compute_open_water``
  144-175, ``compute_force`` 217-268).
* ``.../propeller_models/wagenigen/wagenigen.py`` 3-127.
* ``.../fins/fins_auv_physical/fins_auv_physical.py`` 165-322 and
  ``fins_base.py``; ``config/dataclass/fins/fins_auv_physical_params.py``.
* ``.../fins/fins_auv_actuation_model/fins_auv_actuation_vsim.py`` 5-109.
* ``.../electrical_outboard_motor/electrical_outboard_motor.py`` 6-448
  (working tree; the A-20 edits are uncommitted): throttle ``step`` 105-184,
  rpm ``_prop_thrust_torque`` 311-325 (reverse sign, lines 317-319), rpm
  ``step`` 327-433; ``config/dataclass/thruster/electrical_outboard_motor_params.py``.

Frozen reference: ``tests/data/force_producers/`` (``SOURCE.md``).
"""

import ast
import importlib
import inspect
import re
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "force_producers"
MORE = Path(__file__).resolve().parents[4]
PACKAGE = "more_dynamics.models.force_producers"

G1_TOLERANCE = 1e-9   # 30_checks/README.md, gate G1
G2_TOLERANCE = 1e-10  # 30_checks/README.md, gate G2
G4_FACTOR = 10.0      # 30_checks/README.md, gate G4
SEED = 20261005
N_RANDOM_STATES = 1000

VEHICLE_NAMES = ("remus", "otter", "grethe", "marie", "hugin", "lauv",
                 "mariner", "torqeedo", "cybership", "prestero")

REMUS = "source-sim/MSS/CRAFT/AUV/models/remus100.m"
OTTER = "source-sim/MSS/CRAFT/USV/models/otter.m"
GRAVITY = "source-sim/MSS/INS/functions/gravity.m"
WAG_MAT = "source-sim/MSS/LIBRARY/modeling/utiles/WageningData.mat"
WAG_TXT = "source-sim/MSS/LIBRARY/modeling/utiles/WageningData.txt"
WAG_M = "source-sim/MSS/LIBRARY/modeling/wageningen.m"
GEN_PROP = "more_generic_models/more_generic_models/test/thruster/thruster_wagenigen/test_mss_reference.m"
GEN_DIFF = "more_generic_models/more_generic_models/test/thruster/thruster_differential/test_diff_trhuster.m"
SRC_PROPULSION = "more_generic_models/more_generic_models/dynamics/propulsion/"
SRC_WAG = SRC_PROPULSION + "propeller_models/wagenigen/wagenigen.py"
SRC_PROP = SRC_PROPULSION + "thruster/thruster_wagenigen/thruster_wagenigen.py"
SRC_OUTBOARD = SRC_PROPULSION + "thruster/electrical_outboard_motor/electrical_outboard_motor.py"

CITED_LINES = {
    (OTTER, 90): "g   = 9.81;",
    (OTTER, 105): "y_pont  = 0.395;",
    (OTTER, 132): "l1 = -y_pont;",
    (OTTER, 133): "l2 = y_pont;",
    (OTTER, 134): "k_pos = 0.02216/2;",
    (OTTER, 135): "k_neg = 0.01289/2;",
    (OTTER, 136): "n_max =  sqrt((0.5*24.4 * g)/k_pos);",
    (OTTER, 137): "n_min = -sqrt((0.5*13.6 * g)/k_neg);",
    (OTTER, 220): "n = satlim(n, n_min, n_max);",
    (OTTER, 224): "if n(i) > 0",
    (OTTER, 225): "Thrust(i) = k_pos * n(i) * abs(n(i));",
    (OTTER, 227): "Thrust(i) = k_neg * n(i) * abs(n(i));",
    (OTTER, 232): "tau = [Thrust(1) + Thrust(2) 0 0 0 0 -l1 * Thrust(1) - l2 * Thrust(2) ]';",
    (REMUS, 99): "rho = 1026;",
    (REMUS, 110): "delta_max = deg2rad(20);",
    (REMUS, 111): "n_max = 1525;",
    (REMUS, 114): "delta_r = sat(ui(1), delta_max);",
    (REMUS, 115): "delta_s = sat(ui(2), delta_max);",
    (REMUS, 116): "n_p = sat(ui(3),n_max) / 60;",
    (REMUS, 128): "U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );",
    (REMUS, 132): "L_auv = 1.6;",
    (REMUS, 135): "a = 1.0096 * L_auv/2;",
    (REMUS, 149): "D_prop = 0.14;",
    (REMUS, 150): "t_prop = 0.1;",
    (REMUS, 151): "Va = 0.944 * U_r;",
    (REMUS, 154): "Ja_max = 0.6632;",
    (REMUS, 157): "% >> [KT_0, KQ_0] = wageningen(0,1,0.718,3)",
    (REMUS, 158): "KT_0 = 0.4566;",
    (REMUS, 159): "KQ_0 = 0.0700;",
    (REMUS, 160): "% >> [KT_max, KQ_max] = wageningen(0.6632,1,0.718,3)",
    (REMUS, 161): "KT_max = 0.1798;",
    (REMUS, 162): "KQ_max = 0.0312;",
    (REMUS, 168): "if n_p > 0",
    (REMUS, 170): "X_prop = rho * D_prop^4 * ( ...",
    (REMUS, 171): "KT_0 * abs(n_p) * n_p + (KT_max-KT_0)/Ja_max * (Va/D_prop) * abs(n_p) );",
    (REMUS, 172): "K_prop = rho * D_prop^5 * ( ...",
    (REMUS, 173): "KQ_0 * abs(n_p) * n_p + (KQ_max-KQ_0)/Ja_max * (Va/D_prop) * abs(n_p) );",
    (REMUS, 176): "X_prop = rho * D_prop^4 * KT_0 * abs(n_p) * n_p;",
    (REMUS, 177): "K_prop = rho * D_prop^5 * KQ_0 * abs(n_p) * n_p;",
    (REMUS, 180): "S_fin = 0.00665;",
    (REMUS, 183): "CL_delta_r = 0.5;",
    (REMUS, 184): "A_r = 2 * S_fin;",
    (REMUS, 185): "x_r = -a;",
    (REMUS, 188): "CL_delta_s = 0.7;",
    (REMUS, 189): "A_s = 2 * S_fin;",
    (REMUS, 190): "x_s = -a;",
    (REMUS, 232): "U_rh = sqrt( nu_r(1)^2 + nu_r(2)^2 );",
    (REMUS, 233): "U_rv = sqrt( nu_r(1)^2 + nu_r(3)^2 );",
    (REMUS, 236): "X_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r^2;",
    (REMUS, 237): "X_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s^2;",
    (REMUS, 240): "Y_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r;",
    (REMUS, 243): "Z_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s;",
    (REMUS, 247): "tau(1) = (1-t_prop) * X_prop + X_r + X_s;",
    (REMUS, 248): "tau(2) = Y_r;",
    (REMUS, 249): "tau(3) = Z_s;",
    (REMUS, 250): "tau(4) = K_prop / 10;",
    (REMUS, 251): "tau(5) = -x_s * Z_s;",
    (REMUS, 252): "tau(6) = x_r * Y_r;",
    (GRAVITY, 11): "g = 9.7803253359 * ( 1 + 0.001931850400 * sin(mu)^2 ) /...",
    (GRAVITY, 12): "sqrt( 1 - 0.006694384442 * sin(mu)^2 );",
    (GEN_DIFF, 8): "g      = gravity(deg2rad(63.446827));",
    (GEN_PROP, 83): "rho = 1025;",
    (GEN_PROP, 84): "n_max = 1525;",
    (GEN_PROP, 92): "D_prop = 0.14;",
    (GEN_PROP, 93): "t_prop = 0.1;",
    (GEN_PROP, 94): "Va = 0.944 * U;",
    (GEN_PROP, 96): "Ja_max = 0.6632;",
    (GEN_PROP, 98): "KT_0 = 0.4566;",
    (GEN_PROP, 99): "KQ_0 = 0.0700;",
    (GEN_PROP, 101): "KT_max = 0.1798;",
    (GEN_PROP, 102): "KQ_max = 0.0312;",
    (GEN_PROP, 130): "K_prop = K_prop / 10;",
    (WAG_M, 33): "load('WageningData.mat');",
    (WAG_M, 35): "KT = sum(WagCThrust_stuv.*((Ja).^WagThrust_s).*(PD.^WagThrust_t).*...",
    (WAG_M, 37): "KQ = sum(WagCTorque_stuv.*((Ja).^WagTorque_s).*(PD.^WagTorque_t).*...",
    (SRC_WAG, 11): "Ja = np.clip(np.abs(Ja), 0, 1.3)",
    (SRC_PROP, 113): "Ja = np.clip(Ja, 0.0, self.Ja_max * 2)",
    (SRC_OUTBOARD, 315): "Ja = Va / (n_rps * self.D_prop)",
    (SRC_OUTBOARD, 382): "Va = self.Va_factor * nu[0] * (1.0 - self.w_prop)",
    (SRC_OUTBOARD, 385): "T_prop = (1.0 - self.t_prop) * T_prop",
    (SRC_OUTBOARD, 391): "P_shaft_limit = self.max_power * max(self.eta, 1e-3)",
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
    names = {"__builtins__": {}, "deg2rad": np.deg2rad, "pi": np.pi, "sqrt": np.sqrt}
    return float(eval(expr.replace("^", "**"), names, dict(env or {})))  # pinned arithmetic


def _number_in(rel, number, pattern):
    match = re.search(pattern, _line(rel, number))
    assert match, (rel, number, pattern)
    return float(match.group(1))


def _gravity(mu):
    """INS/functions/gravity.m lines 11-12."""
    g0 = _number_in(GRAVITY, 11, r"g = ([0-9.]+) \*")
    k1 = _number_in(GRAVITY, 11, r"1 \+ ([0-9.]+) \*")
    e2 = _number_in(GRAVITY, 12, r"1 - ([0-9.]+) \*")
    return g0 * (1 + k1 * np.sin(mu) ** 2) / np.sqrt(1 - e2 * np.sin(mu) ** 2)


def _contract(name):
    try:
        return importlib.import_module(f"{PACKAGE}.{name}")
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed (owner decision E-7): {exc}")
        pytest.fail(f"block not ported yet: {exc}")


MODULE_OF = {"outboard_motor_rpm": "outboard_motor", "outboard_motor_throttle": "outboard_motor",
             "wageningen": "wageningen_kt_kq"}


def _build(name, params):
    """``preprocess_<name>`` and ``<name>_casadi`` from their module."""
    block = _contract(MODULE_OF.get(name, name))
    constants = getattr(block, f"preprocess_{name}")(**params)
    return constants, getattr(block, f"{name}_casadi")(constants)


def _tau(function, **inputs):
    out = function(**{k: np.asarray(v, float) for k, v in inputs.items()})
    return np.array(out["tau"], dtype=float).reshape(-1)


def _load_csv(name):
    path = DATA_DIR / name
    header = path.read_text().splitlines()[0].split(",")
    return header, np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _rng():
    return np.random.default_rng(SEED)


def _source(kind):
    try:
        if kind == "differential":
            from more_generic_models.dynamics.propulsion.thruster.thruster_differential.thruster_differential import DifferentialThruster as cls
        elif kind == "propeller":
            from more_generic_models.dynamics.propulsion.thruster.thruster_wagenigen.thruster_wagenigen import ThrusterWageningen as cls
        elif kind == "wageningen":
            from more_generic_models.dynamics.propulsion.propeller_models.wagenigen.wagenigen import wageningen as cls
        elif kind == "fins":
            from more_generic_models.dynamics.propulsion.fins.fins_auv_physical.fins_auv_physical import FinsAUVPhysical as cls
        elif kind == "vsim":
            from more_generic_models.dynamics.propulsion.fins.fins_auv_actuation_model.fins_auv_actuation_vsim import FinsActuationVSIM as cls
        elif kind == "outboard_rpm":
            from more_generic_models.dynamics.propulsion.thruster.electrical_outboard_motor.electrical_outboard_motor import ElectricOutboardMotorRPM as cls
        elif kind == "outboard_throttle":
            from more_generic_models.dynamics.propulsion.thruster.electrical_outboard_motor.electrical_outboard_motor import ElectricOutboardMotor as cls
        else:
            raise KeyError(kind)
    except ImportError as exc:
        pytest.skip(
            "numpy source not importable; install it into the venv with "
            f"pip install -e <more>/more_generic_models (job A-4c) — {exc}"
        )
    return cls


# --------------------------------------------------------------------------
# Pinned lines (no port needed)
# --------------------------------------------------------------------------
def test_cited_lines_are_unchanged():
    for (rel, number), text in CITED_LINES.items():
        line = _line(rel, number)
        assert line.startswith(text), (rel, number, line)


# ==========================================================================
# 1. Differential thruster (Otter-type twin propellers)
# ==========================================================================
def _otter_differential(gravity):
    """otter.m 105, 132-137; positions [left, right] = [0, -/+ y_pont, 0]."""
    y = _value(OTTER, 105)
    k_pos = _value(OTTER, 134)
    k_neg = _value(OTTER, 135)
    f_fwd = _number_in(OTTER, 136, r"\(\(0\.5\*([0-9.]+) \* g\)")
    f_rev = _number_in(OTTER, 137, r"\(\(0\.5\*([0-9.]+) \* g\)")
    return {
        "positive_thrust_coefficients": [k_pos, k_pos],
        "negative_thrust_coefficients": [k_neg, k_neg],
        "thruster_positions": [[0.0, _value(OTTER, 132, {"y_pont": y}), 0.0],
                               [0.0, _value(OTTER, 133, {"y_pont": y}), 0.0]],
        "thruster_directions": [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        "max_forward_thrust": [0.5 * f_fwd * gravity] * 2,
        "max_reverse_thrust": [0.5 * f_rev * gravity] * 2,
    }


def _scaled_differential():
    """Second vehicle (test construction): unequal, bigger propellers, a toe-in
    angle and an x offset."""
    p = _otter_differential(9.81)
    return {**p,
            "positive_thrust_coefficients": [4.0 * p["positive_thrust_coefficients"][0],
                                             3.5 * p["positive_thrust_coefficients"][1]],
            "negative_thrust_coefficients": [4.0 * p["negative_thrust_coefficients"][0],
                                             3.0 * p["negative_thrust_coefficients"][1]],
            "thruster_positions": [[-1.2, -0.9, 0.3], [-1.2, 0.9, 0.3]],
            "thruster_directions": [[1.0, 0.05, 0.0], [1.0, -0.05, 0.0]],
            "max_forward_thrust": [400.0, 380.0],
            "max_reverse_thrust": [200.0, 180.0]}


def _mss_otter_tau(p, n):
    """otter.m 220-232 (left = n(1), lever arm l1)."""
    k_pos, k_neg = p["positive_thrust_coefficients"][0], p["negative_thrust_coefficients"][0]
    n_max = np.sqrt(p["max_forward_thrust"][0] / k_pos)
    n_min = -np.sqrt(p["max_reverse_thrust"][0] / k_neg)
    n = np.clip(n, n_min, n_max)
    thrust = np.where(n > 0, k_pos * n * np.abs(n), k_neg * n * np.abs(n))
    l1, l2 = p["thruster_positions"][0][1], p["thruster_positions"][1][1]
    return np.array([thrust[0] + thrust[1], 0, 0, 0, 0, -l1 * thrust[0] - l2 * thrust[1]])


def _differential_csv():
    h, v = _load_csv("differential_thruster_matlab.csv")
    get = lambda *names: v[:, [h.index(n) for n in names]]
    return {"n": get("n_cmd_left", "n_cmd_right"), "thrust": get("thrust_left", "thrust_right"),
            "tau_X": v[:, h.index("tau_X")], "tau_N": v[:, h.index("tau_N")]}


def _differential_source(p):
    cls = _source("differential")
    k_pos = np.asarray(p["positive_thrust_coefficients"])
    k_neg = np.asarray(p["negative_thrust_coefficients"])
    return cls(params={
        "k_pos": k_pos, "k_neg": k_neg,
        "r_thruster_left": p["thruster_positions"][0], "r_thruster_right": p["thruster_positions"][1],
        "d_thruster_left": p["thruster_directions"][0], "d_thruster_right": p["thruster_directions"][1],
        "n_max_override": np.sqrt(np.asarray(p["max_forward_thrust"]) / k_pos),
        "n_min_override": -np.sqrt(np.asarray(p["max_reverse_thrust"]) / k_neg),
    })


def _matlab_gravity_differential():
    return _gravity(np.deg2rad(_number_in(GEN_DIFF, 8, r"deg2rad\(([0-9.]+)\)")))


def test_G1_numpy_source_differential_matches_matlab():
    source = _source("differential")()   # defaults: k = otter.m / 2, limits with gravity(mu)
    ref = _differential_csv()
    for k, n in enumerate(ref["n"]):
        _, thrust = source.thrust_from_speed(n)
        tau = source.tau_from_thrust(thrust)
        assert _max_diff(thrust, ref["thrust"][k]) <= G1_TOLERANCE, k
        assert _max_diff([tau[0], tau[5]], [ref["tau_X"][k], ref["tau_N"][k]]) <= G1_TOLERANCE, k


def test_G1_block_differential_matches_matlab():
    _, function = _build("differential_thruster", _otter_differential(_matlab_gravity_differential()))
    ref = _differential_csv()
    for k, n in enumerate(ref["n"]):
        tau = _tau(function, n=n, nu_r=np.zeros(6))
        assert _max_diff(tau[[0, 5]], [ref["tau_X"][k], ref["tau_N"][k]]) <= G1_TOLERANCE, k
        assert _max_diff(tau[1:5], 0.0) <= G1_TOLERANCE, k


def test_G2_differential_signature():
    constants, function = _build("differential_thruster", _otter_differential(9.81))
    assert function.name_in() == ["n", "nu_r"] and function.name_out() == ["tau"]
    assert function.size_in(0) == (2, 1) and function.size_in(1) == (6, 1)
    assert function.size_out(0) == (6, 1)
    for field in ("max_speed", "min_speed", "allocation_matrix"):
        assert hasattr(constants, field), field


@pytest.mark.parametrize("parameter_set", ["otter_like", "scaled"])
def test_G2_differential_matches_source(parameter_set):
    p = _otter_differential(_matlab_gravity_differential()) if parameter_set == "otter_like" \
        else _scaled_differential()
    source = _differential_source(p)
    constants, function = _build("differential_thruster", p)
    np.testing.assert_allclose(constants.allocation_matrix, source.get_B(), atol=G2_TOLERANCE, rtol=0)
    rng = _rng()
    cases = np.vstack([_differential_csv()["n"], rng.uniform(-150.0, 150.0, size=(N_RANDOM_STATES, 2))])
    nu_r = rng.uniform(-3.0, 3.0, size=(len(cases), 6))
    for k, n in enumerate(cases):
        _, thrust = source.thrust_from_speed(n)
        assert _max_diff(_tau(function, n=n, nu_r=nu_r[k]), source.tau_from_thrust(thrust)) <= G2_TOLERANCE, k


def test_G4_differential_left_coefficient_plus_1_percent_is_detected():
    p = _otter_differential(_matlab_gravity_differential())
    ref = _differential_csv()

    def worst(params):
        _, f = _build("differential_thruster", params)
        return max(_max_diff(_tau(f, n=n, nu_r=np.zeros(6))[[0, 5]], [ref["tau_X"][k], ref["tau_N"][k]])
                   for k, n in enumerate(ref["n"]))

    assert worst(p) <= G1_TOLERANCE
    k_pos = list(p["positive_thrust_coefficients"])
    k_pos[0] *= 1.01
    assert worst({**p, "positive_thrust_coefficients": k_pos}) > G4_FACTOR * G1_TOLERANCE


def test_G1_MSS_differential_default_equals_matlab():
    """otter.m's own g = 9.81 (line 90) vs MATLAB running otter.m (A-26): the
    saturation limits are the extremes of the CSV's (post-saturation) speeds,
    and tau_X, tau_N match on all 1201 rows."""
    h, v = _load_csv("differential_thruster_mss_current.csv")
    n = v[:, [h.index("n_cmd_left"), h.index("n_cmd_right")]]
    constants, function = _build("differential_thruster", _otter_differential(_value(OTTER, 90)))
    np.testing.assert_allclose(constants.max_speed, [n.max()] * 2, atol=G1_TOLERANCE, rtol=0)
    np.testing.assert_allclose(constants.min_speed, [n.min()] * 2, atol=G1_TOLERANCE, rtol=0)
    for k in range(len(n)):
        tau = _tau(function, n=n[k], nu_r=np.zeros(6))
        assert _max_diff(tau[[0, 5]], v[k, [h.index("tau_X"), h.index("tau_N")]]) <= G1_TOLERANCE, k
        assert _max_diff(tau[1:5], 0.0) <= G1_TOLERANCE, k


def test_MSS_differential_equals_otter():
    """otter.m 132-137, 220-232 with otter.m's own g = 9.81 (line 90)."""
    p = _otter_differential(_value(OTTER, 90))
    _, function = _build("differential_thruster", p)
    for k, n in enumerate(_rng().uniform(-150.0, 150.0, size=(N_RANDOM_STATES, 2))):
        assert _max_diff(_tau(function, n=n, nu_r=np.zeros(6)), _mss_otter_tau(p, n)) <= G1_TOLERANCE, (
            "otter.m 220-232", k)


# ==========================================================================
# 2. Propeller (Wageningen B-series, single screw)
# ==========================================================================
OPEN_WATER_MODES = {"linearized": "linearized_open_water", "full": "full_open_water",
                    "bollard": "bollard_pull"}


def _matlab_propeller():
    """test_mss_reference.m 83-102, 130: remus100.m's propeller with rho = 1025."""
    return {
        "diameter": _value(GEN_PROP, 92), "max_speed": _value(GEN_PROP, 84),
        "thrust_deduction": _value(GEN_PROP, 93),
        "wake_fraction": 1.0 - _number_in(GEN_PROP, 94, r"= ([0-9.]+) \* U"),
        "pitch_diameter_ratio": 1.0, "blade_area_ratio": 0.718, "blade_count": 3,  # remus100.m 157
        "max_advance_number": _value(GEN_PROP, 96),
        "roll_moment_scale": 1.0 / _number_in(GEN_PROP, 130, r"/ ([0-9.]+);"),
        "water_density": _value(GEN_PROP, 83),
        "position": [0.0, 0.0, 0.0], "orientation": [0.0, 0.0, 0.0],
        "open_water_model": "linearized",
        "thrust_torque_coefficients": [_value(GEN_PROP, n) for n in (98, 99, 101, 102)],
    }


def _remus_propeller():
    """remus100.m 99, 111, 149-162, 250 (rho = 1026); geometry of line 157."""
    pd, aeao, z = (float(t) for t in re.search(
        r"wageningen\(0,([0-9.]+),([0-9.]+),([0-9]+)\)", _line(REMUS, 157)).groups())
    return {
        "diameter": _value(REMUS, 149), "max_speed": _value(REMUS, 111),
        "thrust_deduction": _value(REMUS, 150),
        "wake_fraction": 1.0 - _number_in(REMUS, 151, r"= ([0-9.]+) \* U_r"),
        "pitch_diameter_ratio": pd, "blade_area_ratio": aeao, "blade_count": int(z),
        "max_advance_number": _value(REMUS, 154),
        "roll_moment_scale": 1.0 / _number_in(REMUS, 250, r"K_prop / ([0-9.]+);"),
        "water_density": _value(REMUS, 99),
        "position": [0.0, 0.0, 0.0], "orientation": [0.0, 0.0, 0.0],
        "open_water_model": "linearized",
        "thrust_torque_coefficients": [_value(REMUS, n) for n in (158, 159, 161, 162)],
    }


def _scaled_propeller():
    """Second vehicle (test construction): a bigger 4-blade propeller mounted
    aft, below and slightly tilted."""
    p = _remus_propeller()
    return {**p, "diameter": 0.3, "max_speed": 800.0, "thrust_deduction": 0.15,
            "wake_fraction": 0.1, "pitch_diameter_ratio": 0.9, "blade_area_ratio": 0.55,
            "blade_count": 4, "max_advance_number": 0.8, "roll_moment_scale": 0.2,
            "water_density": 1025.0, "position": [-0.8, 0.0, 0.05],
            "orientation": [0.0, 0.1, 0.05], "thrust_torque_coefficients": None}


def _mss_remus_propeller(p, rpm, nu_r):
    """remus100.m 116, 128, 151, 168-178, 247 (propeller part), 250."""
    kt0, kq0, ktm, kqm = p["thrust_torque_coefficients"]
    n_p = np.clip(rpm, -p["max_speed"], p["max_speed"]) / 60.0
    va = (1.0 - p["wake_fraction"]) * np.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2)
    rho, d, j_max = p["water_density"], p["diameter"], p["max_advance_number"]
    if n_p > 0:
        x = rho * d ** 4 * (kt0 * abs(n_p) * n_p + (ktm - kt0) / j_max * (va / d) * abs(n_p))
        k = rho * d ** 5 * (kq0 * abs(n_p) * n_p + (kqm - kq0) / j_max * (va / d) * abs(n_p))
    else:
        x = rho * d ** 4 * kt0 * abs(n_p) * n_p
        k = rho * d ** 5 * kq0 * abs(n_p) * n_p
    return np.array([(1 - p["thrust_deduction"]) * x, 0, 0, p["roll_moment_scale"] * k, 0, 0])


def _propeller_source(p):
    cls = _source("propeller")
    # the mode goes in the config: ThrusterWageningen.__init__ (line 76) re-reads it from there
    src = cls(params={
        "n_max": p["max_speed"], "D_prop": p["diameter"], "t_prop": p["thrust_deduction"],
        "rho": p["water_density"], "PD": p["pitch_diameter_ratio"], "AEAO": p["blade_area_ratio"],
        "z": p["blade_count"], "Ja_max": p["max_advance_number"], "scale_roll": p["roll_moment_scale"],
        "w_factor": p["wake_fraction"], "mode": OPEN_WATER_MODES[p["open_water_model"]],
        "position": tuple(p["position"]), "orientation_rpy": tuple(p["orientation"])})
    assert src.mode == OPEN_WATER_MODES[p["open_water_model"]]
    if p["thrust_torque_coefficients"] is not None:   # exact printed values, as the block's override
        src.KT_0, src.KQ_0, src.KT_max, src.KQ_max = p["thrust_torque_coefficients"]
    return src


def _propeller_csv():
    h, v = _load_csv("propeller_matlab.csv")
    return {"rpm": v[:, h.index("RPM")], "U": v[:, h.index("U")],
            "tau": v[:, [h.index(f"tau{i}") for i in range(1, 7)]]}


def test_propeller_reference_is_the_linearised_remus_model():
    """The stored CSV is remus100.m 168-178 with test_mss_reference.m's rho."""
    p, ref = _matlab_propeller(), _propeller_csv()
    for k in range(len(ref["rpm"])):
        expected = _mss_remus_propeller(p, ref["rpm"][k], [ref["U"][k], 0, 0, 0, 0, 0])
        assert _max_diff(expected, ref["tau"][k]) <= G1_TOLERANCE, k


def test_G1_numpy_source_propeller_matches_matlab():
    source = _propeller_source(_matlab_propeller())
    ref = _propeller_csv()
    for k in range(len(ref["rpm"])):
        tau = source.compute_force(ref["rpm"][k], np.array([ref["U"][k], 0, 0, 0, 0, 0]))
        assert _max_diff(tau, ref["tau"][k]) <= G1_TOLERANCE, k


def test_G1_block_propeller_matches_matlab():
    _, function = _build("propeller", _matlab_propeller())
    ref = _propeller_csv()
    for k in range(len(ref["rpm"])):
        tau = _tau(function, n=[ref["rpm"][k]], nu_r=[ref["U"][k], 0, 0, 0, 0, 0])
        assert _max_diff(tau, ref["tau"][k]) <= G1_TOLERANCE, k


def test_G2_propeller_signature():
    constants, function = _build("propeller", _remus_propeller())
    assert function.name_in() == ["n", "nu_r"] and function.name_out() == ["tau"]
    assert function.size_in(0) == (1, 1) and function.size_in(1) == (6, 1)
    assert function.size_out(0) == (6, 1)
    np.testing.assert_allclose(constants.thrust_torque_coefficients,
                               _remus_propeller()["thrust_torque_coefficients"], atol=0, rtol=0)


@pytest.mark.parametrize("mode", sorted(OPEN_WATER_MODES))
@pytest.mark.parametrize("parameter_set", ["remus_like", "scaled"])
def test_G2_propeller_matches_source(parameter_set, mode):
    """On the flag path (A-28): the numpy source clamps J in the polynomial."""
    base = _remus_propeller() if parameter_set == "remus_like" else _scaled_propeller()
    p = {**base, "open_water_model": mode, "thrust_torque_coefficients": None}
    source = _propeller_source(p)
    constants, function = _build("propeller", {**p, "clip_advance_ratio": True})
    np.testing.assert_allclose(constants.thrust_torque_coefficients,
                               [source.KT_0, source.KQ_0, source.KT_max, source.KQ_max],
                               atol=G2_TOLERANCE, rtol=0)
    rng = _rng()
    ref = _propeller_csv()
    rpm = np.concatenate([ref["rpm"], rng.uniform(-2000.0, 2000.0, N_RANDOM_STATES)])
    nu_r = np.zeros((len(rpm), 6))
    nu_r[: len(ref["U"]), 0] = ref["U"]
    nu_r[len(ref["U"]):] = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    for k in range(len(rpm)):
        expected = source.compute_force(rpm[k], nu_r[k])
        assert _max_diff(_tau(function, n=[rpm[k]], nu_r=nu_r[k]), expected) <= G2_TOLERANCE, (mode, k)


def test_MSS_propeller_full_model_default_beyond_J_1_3_uses_unclamped_polynomial():
    """A-28 (E-22 Q1 a). Test construction: the scaled propeller on the shaft
    axis at the CO (position and orientation zero), "full" model, J = 1.5,
    inside the model's own clip 2 x 0.8 = 1.6 (thruster_wagenigen.py 113) and
    beyond 1.3 (wagenigen.py 11). Default: K_T, K_Q = MSS wageningen.m at J;
    flag: the numpy source (clamped)."""
    _line(SRC_PROP, 113)
    _line(SRC_WAG, 11)
    p = {**_scaled_propeller(), "position": [0.0, 0.0, 0.0], "orientation": [0.0, 0.0, 0.0],
         "open_water_model": "full", "thrust_torque_coefficients": None}
    j, rpm = 1.5, 400.0
    assert 1.3 < j < 2.0 * p["max_advance_number"]
    n_rps = rpm / 60.0
    nu_r = np.array([j * n_rps * p["diameter"] / (1.0 - p["wake_fraction"]), 0, 0, 0, 0, 0])
    kt, kq = _mss_wageningen(j, p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"])
    rho, d = p["water_density"], p["diameter"]
    expected = np.array([(1.0 - p["thrust_deduction"]) * rho * d ** 4 * kt * n_rps ** 2, 0, 0,
                         p["roll_moment_scale"] * rho * d ** 5 * kq * n_rps ** 2, 0, 0])
    default = _tau(_build("propeller", p)[1], n=[rpm], nu_r=nu_r)
    assert _max_diff(default, expected) <= G1_TOLERANCE, ("wageningen.m at J = 1.5", default, expected)
    flagged = _tau(_build("propeller", {**p, "clip_advance_ratio": True})[1], n=[rpm], nu_r=nu_r)
    assert _max_diff(flagged, _propeller_source(p).compute_force(rpm, nu_r)) <= G2_TOLERANCE
    assert _max_diff(default, flagged) > G4_FACTOR * G1_TOLERANCE


def test_G4_propeller_diameter_plus_1_percent_is_detected():
    p, ref = _matlab_propeller(), _propeller_csv()

    def worst(params):
        _, f = _build("propeller", params)
        return max(_max_diff(_tau(f, n=[ref["rpm"][k]], nu_r=[ref["U"][k], 0, 0, 0, 0, 0]), ref["tau"][k])
                   for k in range(len(ref["rpm"])))

    assert worst(p) <= G1_TOLERANCE
    assert worst({**p, "diameter": 1.01 * p["diameter"]}) > G4_FACTOR * G1_TOLERANCE


def test_MSS_propeller_equals_remus100():
    """remus100.m 99 (rho 1026), 111, 116, 149-178, 247, 250 with the printed
    coefficients (158-162) as the exact override."""
    p = _remus_propeller()
    _, function = _build("propeller", p)
    rng = _rng()
    for k in range(N_RANDOM_STATES):
        rpm, nu_r = rng.uniform(-2000.0, 2000.0), rng.uniform(-3.0, 3.0, 6)
        expected = _mss_remus_propeller(p, rpm, nu_r)
        assert _max_diff(_tau(function, n=[rpm], nu_r=nu_r), expected) <= G1_TOLERANCE, ("remus100.m", k)


def test_G5_polynomial_coefficients_round_to_remus100_printed_values():
    """remus100.m 157-162 print KT/KQ at J = 0 and J = Ja_max from
    ``wageningen(J,1,0.718,3)`` to 4 decimals; the polynomial must land on them."""
    wageningen = _source("wageningen")
    p = _remus_propeller()
    printed = p["thrust_torque_coefficients"]
    kt0, kq0 = wageningen(0.0, p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"])
    ktm, kqm = wageningen(p["max_advance_number"], p["pitch_diameter_ratio"],
                          p["blade_area_ratio"], p["blade_count"])
    np.testing.assert_array_equal(np.round([kt0, kq0, ktm, kqm], 4), printed,
                                  err_msg="remus100.m lines 157-162 (printed precision)")


def test_G5_block_polynomial_coefficients_round_to_remus100_printed_values():
    p = {**_remus_propeller(), "thrust_torque_coefficients": None}
    constants, _ = _build("propeller", p)
    np.testing.assert_array_equal(np.round(constants.thrust_torque_coefficients, 4),
                                  _remus_propeller()["thrust_torque_coefficients"],
                                  err_msg="remus100.m lines 157-162 (printed precision)")


# ==========================================================================
# 3. Wageningen K_T / K_Q polynomial
# ==========================================================================
def _mss_wageningen_tables():
    """WageningData.mat, loaded as wageningen.m line 33 does."""
    import scipy.io
    path = MORE / WAG_MAT
    if not path.exists():
        pytest.skip(f"MSS checkout not found at {path}")
    _line(WAG_M, 33)
    d = scipy.io.loadmat(path)
    thrust = np.hstack([d["WagCThrust_stuv"], d["WagThrust_s"], d["WagThrust_t"],
                        d["WagThrust_u"], d["WagThrust_v"]])
    torque = np.hstack([d["WagCTorque_stuv"], d["WagTorque_s"], d["WagTorque_t"],
                        d["WagTorque_u"], d["WagTorque_v"]])
    return thrust, torque


def _mss_wageningen(j, pd, aeao, z):
    """wageningen.m: KT = sum(C .* J.^s .* PD.^t .* AEAO.^u .* z.^v), same for KQ."""
    thrust, torque = _mss_wageningen_tables()
    terms = lambda t: t[:, 0] * j ** t[:, 1] * pd ** t[:, 2] * aeao ** t[:, 3] * z ** t[:, 4]
    return float(np.sum(terms(thrust))), float(np.sum(terms(torque)))


WAGENINGEN_SETS = {"remus_like": (1.0, 0.718, 3), "four_blade": (0.83, 0.55, 4)}


def test_G5_mss_table_equals_its_printed_text():
    """WageningData.mat equals WageningData.txt ("Re = 2 x 10^6" table). The
    paper both transcribe (Barnitsas, Ray and Kinley 1981, as MSS spells it) is
    not on disk: the published-table check proper is owed (needs access)."""
    path = MORE / WAG_TXT
    if not path.exists():
        pytest.skip(f"MSS checkout not found at {path}")
    text = path.read_text()
    rows = lambda s: np.array([[float(x) for x in ln.split()] for ln in s.splitlines()
                               if re.match(r"\s*[-+][0-9]", ln)])
    thrust_txt, torque_txt = (rows(part) for part in text.split("TORQUE:"))
    thrust, torque = _mss_wageningen_tables()
    assert thrust.shape == (39, 5) and torque.shape == (47, 5)
    np.testing.assert_array_equal(thrust, thrust_txt)
    np.testing.assert_array_equal(torque, torque_txt)


@pytest.mark.parametrize("name", sorted(WAGENINGEN_SETS))
def test_G1_numpy_source_wageningen_matches_mss_table(name):
    wageningen = _source("wageningen")
    pd, aeao, z = WAGENINGEN_SETS[name]
    for j in np.linspace(0.0, 1.3, 131):
        assert _max_diff(wageningen(j, pd, aeao, z), _mss_wageningen(j, pd, aeao, z)) <= G1_TOLERANCE, j


def _build_wageningen(pd, aeao, z, **flag):
    return _build("wageningen", {"pitch_diameter_ratio": pd, "blade_area_ratio": aeao, "blade_count": z,
                                 **flag})


def _wageningen_grid():
    """A-28: J on -0.5 ... 2.0, past both ends of the source's [0, 1.3]."""
    return np.concatenate([np.linspace(-0.5, 2.0, 251), _rng().uniform(-0.5, 2.0, N_RANDOM_STATES)])


def _kt_kq(function, j):
    out = function(J=float(j))
    return float(out["KT"]), float(out["KQ"])


def test_G2_wageningen_signature():
    _, function = _build_wageningen(*WAGENINGEN_SETS["remus_like"])
    assert function.name_in() == ["J"] and function.name_out() == ["KT", "KQ"]
    assert function.size_in(0) == (1, 1)
    assert function.size_out(0) == (1, 1) and function.size_out(1) == (1, 1)


@pytest.mark.parametrize("name", sorted(WAGENINGEN_SETS))
def test_G2_wageningen_flag_matches_source_including_clamp(name):
    """``clip_advance_ratio=True`` is the numpy source (wagenigen.py 11)."""
    _line(SRC_WAG, 11)
    wageningen = _source("wageningen")
    pd, aeao, z = WAGENINGEN_SETS[name]
    _, function = _build_wageningen(pd, aeao, z, clip_advance_ratio=True)
    for j in _wageningen_grid():
        assert _max_diff(_kt_kq(function, j), wageningen(j, pd, aeao, z)) <= G2_TOLERANCE, j


@pytest.mark.parametrize("name", sorted(WAGENINGEN_SETS))
def test_G4_wageningen_flag_is_not_dead(name):
    """Default and flag agree on 0 <= J <= 1.3 and differ outside it."""
    pd, aeao, z = WAGENINGEN_SETS[name]
    _, default = _build_wageningen(pd, aeao, z)
    _, flagged = _build_wageningen(pd, aeao, z, clip_advance_ratio=True)
    grid = _wageningen_grid()
    outside = grid[(grid < -1e-3) | (grid > 1.3 + 1e-3)]
    inside = grid[(grid >= 0.0) & (grid <= 1.3)]
    assert len(outside) > 0 and len(inside) > 0
    for j in inside:
        assert _max_diff(_kt_kq(default, j), _kt_kq(flagged, j)) <= G2_TOLERANCE, j
    for j in outside:
        assert _max_diff(_kt_kq(default, j), _kt_kq(flagged, j)) > G4_FACTOR * G1_TOLERANCE, j


def test_G4_wageningen_pitch_ratio_plus_1_percent_is_detected():
    pd, aeao, z = WAGENINGEN_SETS["remus_like"]
    grid = np.linspace(0.0, 1.3, 27)
    worst = lambda f: max(_max_diff(_kt_kq(f, j), _mss_wageningen(j, pd, aeao, z)) for j in grid)
    assert worst(_build_wageningen(pd, aeao, z)[1]) <= G1_TOLERANCE
    assert worst(_build_wageningen(1.01 * pd, aeao, z)[1]) > G4_FACTOR * G1_TOLERANCE


@pytest.mark.parametrize("name", sorted(WAGENINGEN_SETS))
def test_MSS_wageningen_default_equals_mss_at_any_J(name):
    """A-28 (E-22 Q1 a): the default is wageningen.m 35-38 at J itself, also
    for J < 0 and J > 1.3 (MSS has no |J| and no clamp)."""
    _line(WAG_M, 35)
    _line(WAG_M, 37)
    pd, aeao, z = WAGENINGEN_SETS[name]
    _, function = _build_wageningen(pd, aeao, z)
    for j in _wageningen_grid():
        assert _max_diff(_kt_kq(function, j), _mss_wageningen(j, pd, aeao, z)) <= G1_TOLERANCE, j


# ==========================================================================
# 4. Fins: rudder and stern planes (remus100.m)
# ==========================================================================
CONVENTIONS = ("starboard_down_positive", "starboard_up_positive", "port_down_positive", "port_up_positive")


def _remus_fins():
    """remus100.m 99, 110, 135, 180-190."""
    a = _value(REMUS, 135, {"L_auv": _value(REMUS, 132)})
    s_fin = _value(REMUS, 180)
    return {"rudder_area": _value(REMUS, 184, {"S_fin": s_fin}),
            "stern_plane_area": _value(REMUS, 189, {"S_fin": s_fin}),
            "rudder_lift_coefficient": _value(REMUS, 183),
            "stern_plane_lift_coefficient": _value(REMUS, 188),
            "rudder_position": _value(REMUS, 185, {"a": a}),
            "stern_plane_position": _value(REMUS, 190, {"a": a}),
            "max_deflection": _value(REMUS, 110),
            "water_density": _value(REMUS, 99)}


def _template_fins():
    """The numpy source's defaults (fins_auv_physical_params.py): rho 1025,
    x = -0.8, 25 deg (the template departures of A-18)."""
    src = _source("fins")()
    return {"rudder_area": src.Ar, "stern_plane_area": src.Ae,
            "rudder_lift_coefficient": src.CL_delta_r, "stern_plane_lift_coefficient": src.CL_delta_e,
            "rudder_position": src.x_r, "stern_plane_position": src.x_e,
            "max_deflection": src.delta_max, "water_density": src.rho}


def _mss_fins(p, delta, nu_r):
    """remus100.m 114-115, 232-252 (fin part of tau)."""
    dr = np.clip(delta[0], -p["max_deflection"], p["max_deflection"])
    ds = np.clip(delta[1], -p["max_deflection"], p["max_deflection"])
    u_rh = np.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2)
    u_rv = np.sqrt(nu_r[0] ** 2 + nu_r[2] ** 2)
    rho = p["water_density"]
    x_r = -0.5 * rho * u_rh ** 2 * p["rudder_area"] * p["rudder_lift_coefficient"] * dr ** 2
    x_s = -0.5 * rho * u_rv ** 2 * p["stern_plane_area"] * p["stern_plane_lift_coefficient"] * ds ** 2
    y_r = -0.5 * rho * u_rh ** 2 * p["rudder_area"] * p["rudder_lift_coefficient"] * dr
    z_s = -0.5 * rho * u_rv ** 2 * p["stern_plane_area"] * p["stern_plane_lift_coefficient"] * ds
    return np.array([x_r + x_s, y_r, z_s, 0.0, -p["stern_plane_position"] * z_s, p["rudder_position"] * y_r])


def _fin_cases():
    rng = _rng()
    delta = rng.uniform(-0.6, 0.6, size=(N_RANDOM_STATES, 2))   # beyond both limits
    nu_r = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    return delta, nu_r


def test_G2_fins_signature():
    _, function = _build("fins", _remus_fins())
    assert function.name_in() == ["delta", "nu_r"] and function.name_out() == ["tau"]
    assert function.size_in(0) == (2, 1) and function.size_in(1) == (6, 1)


@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("parameter_set", ["template", "scaled"])
def test_G2_fins_match_source(parameter_set, convention):
    p = _template_fins()
    if parameter_set == "scaled":
        p = {**p, "rudder_area": 4 * p["rudder_area"], "stern_plane_area": 3 * p["stern_plane_area"],
             "rudder_position": 2 * p["rudder_position"], "stern_plane_position": 1.8 * p["stern_plane_position"],
             "max_deflection": 0.3}
    cls = _source("fins")
    source = cls(S_fin=0.5 * p["rudder_area"], CL_delta_r=p["rudder_lift_coefficient"],
                 CL_delta_e=p["stern_plane_lift_coefficient"], x_r=p["rudder_position"],
                 x_e=p["stern_plane_position"], delta_max=np.rad2deg(p["max_deflection"]),
                 rho=p["water_density"], convention=convention)
    source.Ae = p["stern_plane_area"]   # the source ties both areas to S_fin; the block does not
    _, function = _build("fins", {**p, "convention": convention})
    for k, (delta, nu_r) in enumerate(zip(*_fin_cases())):
        expected = source.compute_force(nu_r, np.zeros(6), np.zeros(3), delta)
        assert _max_diff(_tau(function, delta=delta, nu_r=nu_r), expected) <= G2_TOLERANCE, (convention, k)


def test_G4_fins_rudder_area_plus_1_percent_is_detected():
    p = _remus_fins()
    delta, nu_r = _fin_cases()
    worst = lambda params: max(_max_diff(_tau(_build("fins", params)[1], delta=d, nu_r=n), _mss_fins(p, d, n))
                               for d, n in zip(delta[:100], nu_r[:100]))
    assert worst(p) <= G1_TOLERANCE
    assert worst({**p, "rudder_area": 1.01 * p["rudder_area"]}) > G4_FACTOR * G1_TOLERANCE


def test_MSS_fins_default_convention_equals_remus100():
    p = _remus_fins()
    _, function = _build("fins", p)    # default convention
    for k, (delta, nu_r) in enumerate(zip(*_fin_cases())):
        assert _max_diff(_tau(function, delta=delta, nu_r=nu_r), _mss_fins(p, delta, nu_r)) <= G1_TOLERANCE, (
            "remus100.m 114-115, 232-252", k)


def test_departure_template_fin_values_are_parameters_not_code():
    """A-18: the template's 25 deg limit (and rho 1025, x = -0.8) differ from
    remus100.m; the block reproduces both from parameters alone."""
    p_mss, p_tmpl = _remus_fins(), _template_fins()
    assert p_tmpl["max_deflection"] != p_mss["max_deflection"]
    delta, nu_r = np.array([0.4, -0.4]), np.array([1.5, 0.1, -0.1, 0, 0, 0])
    tau_mss = _tau(_build("fins", p_mss)[1], delta=delta, nu_r=nu_r)
    tau_tmpl = _tau(_build("fins", p_tmpl)[1], delta=delta, nu_r=nu_r)
    assert _max_diff(tau_mss, _mss_fins(p_mss, delta, nu_r)) <= G1_TOLERANCE
    assert _max_diff(tau_tmpl, _mss_fins(p_tmpl, delta, nu_r)) <= G1_TOLERANCE
    assert _max_diff(tau_mss, tau_tmpl) > G4_FACTOR * G1_TOLERANCE


# ==========================================================================
# 5. VSIM fins (numpy source only)
# ==========================================================================
def _vsim_parameters(scaled):
    src = _source("vsim")()
    p = {"max_forces": src.max_force, "positions": src.position, "max_deflection": src.max_act}
    if scaled:
        p = {"max_forces": 2.5 * p["max_forces"], "positions": 1.7 * p["positions"],
             "max_deflection": 0.35}
    return p


def _vsim_source(p):
    return _source("vsim")(max_force=np.asarray(p["max_forces"]), position=np.asarray(p["positions"]),
                           max_act=p["max_deflection"])


def _vsim_cases(n_fins):
    rng = _rng()
    delta = rng.uniform(-0.6, 0.6, size=(N_RANDOM_STATES, n_fins))
    nu_r = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    nu_r[:20, :3] = rng.uniform(-1e-7, 1e-7, size=(20, 3))   # the source's zero-speed branch
    return delta, nu_r


def test_G2_vsim_fins_signature():
    p = _vsim_parameters(False)
    _, function = _build("vsim_fins", p)
    assert function.name_in() == ["delta", "nu_r"] and function.name_out() == ["tau"]
    assert function.size_in(0) == (len(p["max_forces"]), 1)


@pytest.mark.parametrize("scaled", [False, True], ids=["lauv_like", "scaled"])
def test_G2_vsim_fins_match_source(scaled):
    p = _vsim_parameters(scaled)
    source = _vsim_source(p)
    _, function = _build("vsim_fins", p)
    for k, (delta, nu_r) in enumerate(zip(*_vsim_cases(len(p["max_forces"])))):
        assert _max_diff(_tau(function, delta=delta, nu_r=nu_r), source.compute_force(nu_r, delta)) <= G2_TOLERANCE, k


def test_G4_vsim_fins_force_plus_1_percent_is_detected():
    p = _vsim_parameters(False)
    source = _vsim_source(p)
    delta, nu_r = _vsim_cases(len(p["max_forces"]))
    worst = lambda params: max(_max_diff(_tau(_build("vsim_fins", params)[1], delta=d, nu_r=n),
                                         source.compute_force(n, d)) for d, n in zip(delta[:100], nu_r[:100]))
    assert worst(p) <= G2_TOLERANCE
    assert worst({**p, "max_forces": 1.01 * np.asarray(p["max_forces"])}) > G4_FACTOR * G2_TOLERANCE


# ==========================================================================
# 6. Outboard motor (E-20 Q5 a: reverse sign fixed)
# ==========================================================================
def _outboard_rpm_parameters(scaled):
    """The numpy source's defaults (electrical_outboard_motor_params.py, A-20
    values); ``scaled`` is a test construction."""
    c = _source("outboard_rpm")().config
    p = {"max_thrust": c.max_thrust, "max_power": c.max_power, "efficiency": c.efficiency,
         "position": np.asarray(c.r_body, float), "max_speed": c.n_max_rpm,
         "propeller_diameter": c.D_prop, "pitch_diameter_ratio": c.PD,
         "blade_area_ratio": c.AEAO, "blade_count": c.z, "thrust_deduction": c.t_prop,
         "water_density": c.rho, "advance_speed_factor": c.Va_factor, "propwash_factor": c.w_prop,
         "reverse_thrust_factor": c.reverse_KT_factor, "reverse_torque_factor": c.reverse_KQ_factor}
    if scaled:
        p = {**p, "max_thrust": 1500.0, "max_power": 6000.0, "position": np.array([-3.5, 0.4, 0.7]),
             "max_speed": 2000.0, "propeller_diameter": 0.4, "pitch_diameter_ratio": 1.1,
             "blade_count": 4, "reverse_thrust_factor": 0.4, "reverse_torque_factor": 0.7}
    return p


def _outboard_rpm_source(p, power_limit=True):
    cls = _source("outboard_rpm")
    return cls(max_thrust=p["max_thrust"], max_power=p["max_power"], efficiency=p["efficiency"],
               r_body=p["position"], n_max_rpm=p["max_speed"], D_prop=p["propeller_diameter"],
               PD=p["pitch_diameter_ratio"], AEAO=p["blade_area_ratio"], z=p["blade_count"],
               t_prop=p["thrust_deduction"], rho=p["water_density"], Va_factor=p["advance_speed_factor"],
               w_prop=p["propwash_factor"], reverse_KT_factor=p["reverse_thrust_factor"],
               reverse_KQ_factor=p["reverse_torque_factor"], tau_rpm=0.0,
               disable_power_limit=not power_limit)


def _source_force_at(source, rpm, delta, nu_r):
    """One source step that leaves the states where the block reads them: rpm
    lag off (tau_rpm = 0), steering command equal to the current angle (no
    rate), battery full."""
    source.delta = float(delta)
    u_steer = delta / source.delta_max * source._steer_input_sign
    return source.step(rpm, u_steer, np.asarray(nu_r, float), 0.1)


def _source_reverse_sign_is_wrong(p):
    tau = _source_force_at(_outboard_rpm_source(p), -p["max_speed"], 0.0, np.zeros(6))
    return tau[0] > 0.0


def _fixed_reference(p, rpm, delta, nu_r, power_limit=True):
    """The fixed force. In the current source a reverse K_T, K_Q carry an extra
    sign(n) = -1 (electrical_outboard_motor.py 317-319). Thrust and torque are
    odd in that factor; the power cap uses |Q| and the thrust cap is
    symmetric; the vectoring and moment are linear in thrust. So for n < 0 the
    fixed tau is exactly -1 x the current source's tau. Once A-25 fixes the
    source, the probe returns False and the source is used as is."""
    tau = _source_force_at(_outboard_rpm_source(p, power_limit), rpm, delta, nu_r)
    flip = -1.0 if (rpm < 0 and _source_reverse_sign_is_wrong(p)) else 1.0
    return flip * tau


def _outboard_cases(p):
    rng = _rng()
    rpm = rng.uniform(-1.2 * p["max_speed"], 1.2 * p["max_speed"], N_RANDOM_STATES)
    delta = rng.uniform(-0.6, 0.6, N_RANDOM_STATES)   # within the source's 35 deg
    nu_r = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    return rpm, delta, nu_r


def test_E20_Q5a_numpy_source_full_reverse_pushes_backwards():
    """FAILS TODAY ON PURPOSE (owner decision E-20 Q5 a, A-20 finding 2): the
    current source gives +166 N of surge at -1300 rpm. Job A-25 fixes the sign
    in both repositories; this test then passes."""
    p = _outboard_rpm_parameters(False)
    tau = _source_force_at(_outboard_rpm_source(p), -p["max_speed"], 0.0, np.zeros(6))
    assert tau[0] < 0.0, f"numpy source: full reverse gives surge force {tau[0]:+.1f} N (sign bug, E-20 Q5 a)"


def test_E20_Q5a_block_full_reverse_pushes_backwards():
    p = _outboard_rpm_parameters(False)
    _, function = _build("outboard_motor_rpm", p)
    tau = _tau(function, n=[-p["max_speed"]], delta=[0.0], nu_r=np.zeros(6))
    assert tau[0] < 0.0, tau
    forward = _tau(function, n=[p["max_speed"]], delta=[0.0], nu_r=np.zeros(6))
    assert forward[0] > 0.0, forward


def test_G2_outboard_signatures():
    _, rpm_fn = _build("outboard_motor_rpm", _outboard_rpm_parameters(False))
    assert rpm_fn.name_in() == ["n", "delta", "nu_r"] and rpm_fn.name_out() == ["tau"]
    _, thr_fn = _build("outboard_motor_throttle", _outboard_throttle_parameters(False))
    assert thr_fn.name_in() == ["thrust", "delta", "nu_r"] and thr_fn.name_out() == ["tau"]
    for f in (rpm_fn, thr_fn):
        assert f.size_in(0) == (1, 1) and f.size_in(1) == (1, 1) and f.size_in(2) == (6, 1)


@pytest.mark.parametrize("power_limit", [True, False], ids=["power_limit", "no_power_limit"])
@pytest.mark.parametrize("scaled", [False, True], ids=["grethe_like", "scaled"])
def test_G2_outboard_rpm_matches_source_with_fixed_reverse(scaled, power_limit):
    """On the flag path (A-28): the numpy source clamps J in the polynomial."""
    p = _outboard_rpm_parameters(scaled)
    _, function = _build("outboard_motor_rpm", {**p, "power_limit": power_limit, "clip_advance_ratio": True})
    for k, (rpm, delta, nu_r) in enumerate(zip(*_outboard_cases(p))):
        expected = _fixed_reference(p, rpm, delta, nu_r, power_limit)
        assert _max_diff(_tau(function, n=[rpm], delta=[delta], nu_r=nu_r), expected) <= G2_TOLERANCE, k


def test_G4_outboard_propeller_diameter_plus_1_percent_is_detected():
    p = _outboard_rpm_parameters(False)
    cases = list(zip(*_outboard_cases(p)))[:100]
    expected = [_fixed_reference(p, *c) for c in cases]
    worst = lambda params: max(_max_diff(_tau(_build("outboard_motor_rpm", params)[1], n=[r], delta=[d], nu_r=n), e)
                               for (r, d, n), e in zip(cases, expected))
    flagged = {**p, "clip_advance_ratio": True}   # A-28: the source's clamp
    assert worst(flagged) <= G2_TOLERANCE
    assert worst({**flagged, "propeller_diameter": 1.01 * p["propeller_diameter"]}) > G4_FACTOR * G2_TOLERANCE


def _outboard_rpm_force(p, rpm, delta, nu_r, kt, kq, power_limit=True):
    """electrical_outboard_motor.py 315-321, 380-400 and the steered wrench
    (403-414) for given K_T, K_Q, with n > 0 (no reverse factors)."""
    for number in (382, 385, 391):
        _line(SRC_OUTBOARD, number)
    n_rps = np.clip(rpm, -p["max_speed"], p["max_speed"]) / 60.0
    rho, d = p["water_density"], p["propeller_diameter"]
    thrust = (1.0 - p["thrust_deduction"]) * rho * d ** 4 * kt * abs(n_rps) * n_rps
    torque = rho * d ** 5 * kq * abs(n_rps) * n_rps
    shaft_power, limit = abs(2.0 * np.pi * n_rps * torque), p["max_power"] * max(p["efficiency"], 1e-3)
    if power_limit and shaft_power > limit:
        thrust *= limit / shaft_power
    thrust = np.clip(thrust, -p["max_thrust"], p["max_thrust"])
    force = thrust * np.array([np.cos(delta), np.sin(delta), 0.0])
    return np.concatenate([force, np.cross(p["position"], force)])


@pytest.mark.parametrize("scaled", [False, True], ids=["grethe_like", "scaled"])
def test_MSS_outboard_rpm_default_beyond_J_1_3_uses_unclamped_polynomial(scaled):
    """A-28 (E-22 Q1 a). Test construction: J = 1.6 at 300 rpm (surge speed
    from J = Va_factor (1 - w_prop) u / (n D), lines 315 and 382), steering
    0.1 rad. Default: K_T, K_Q = MSS wageningen.m at J; flag: the numpy source."""
    _line(SRC_OUTBOARD, 315)
    p = _outboard_rpm_parameters(scaled)
    j, rpm, delta = 1.6, 300.0, 0.1
    n_rps = rpm / 60.0
    u = j * n_rps * p["propeller_diameter"] / (p["advance_speed_factor"] * (1.0 - p["propwash_factor"]))
    nu_r = np.array([u, 0, 0, 0, 0, 0])
    kt, kq = _mss_wageningen(j, p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"])
    expected = _outboard_rpm_force(p, rpm, delta, nu_r, kt, kq)
    assert abs(expected[0]) < p["max_thrust"]   # the thrust cap does not hide the polynomial
    default = _tau(_build("outboard_motor_rpm", p)[1], n=[rpm], delta=[delta], nu_r=nu_r)
    assert _max_diff(default, expected) <= G1_TOLERANCE, ("wageningen.m at J = 1.6", default, expected)
    flagged = _tau(_build("outboard_motor_rpm", {**p, "clip_advance_ratio": True})[1],
                   n=[rpm], delta=[delta], nu_r=nu_r)
    assert _max_diff(flagged, _fixed_reference(p, rpm, delta, nu_r)) <= G2_TOLERANCE
    assert _max_diff(default, flagged) > G4_FACTOR * G1_TOLERANCE


def _outboard_throttle_parameters(scaled):
    c = _source("outboard_throttle")().config
    p = {"max_thrust": c.max_thrust, "max_propulsive_power": c.max_propulsive_power,
         "position": np.asarray(c.r_body, float)}
    if scaled:
        p = {"max_thrust": 1500.0, "max_propulsive_power": 4000.0, "position": np.array([-3.5, 0.4, 0.7])}
    return p


@pytest.mark.parametrize("scaled", [False, True], ids=["grethe_like", "scaled"])
def test_G2_outboard_throttle_matches_source(scaled):
    p = _outboard_throttle_parameters(scaled)
    cls = _source("outboard_throttle")
    _, function = _build("outboard_motor_throttle", p)
    rng = _rng()
    thrust = rng.uniform(-p["max_thrust"], p["max_thrust"], N_RANDOM_STATES)
    delta = rng.uniform(-0.6, 0.6, N_RANDOM_STATES)
    nu_r = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    nu_r[:50, :2] = rng.uniform(-0.05, 0.05, size=(50, 2))   # below the source's 0.1 m/s speed floor
    for k in range(N_RANDOM_STATES):
        src = cls(max_thrust=p["max_thrust"], max_propulsive_power=p["max_propulsive_power"],
                  r_body=p["position"])
        src.thrust, src.delta = float(thrust[k]), float(delta[k])   # states the block reads
        u_steer = delta[k] / src.delta_max * src._steer_input_sign
        expected = src.step(thrust[k] / p["max_thrust"], u_steer, nu_r[k], 0.1)
        got = _tau(function, thrust=[thrust[k]], delta=[delta[k]], nu_r=nu_r[k])
        assert _max_diff(got, expected) <= G2_TOLERANCE, k


# ==========================================================================
# Generic by construction (owner E-16)
# ==========================================================================
def _code_without_docstrings(module):
    tree = ast.parse(Path(module.__file__).read_text())
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(body, list) and body and isinstance(body[0], ast.Expr) \
                and isinstance(getattr(body[0], "value", None), ast.Constant) \
                and isinstance(body[0].value.value, str):
            body.pop(0)
    return ast.unparse(tree).lower()


GENERIC = {  # module -> {preprocess function -> keywords allowed a default}
    "differential_thruster": {"preprocess_differential_thruster": set()},
    "propeller": {"preprocess_propeller": {"thrust_torque_coefficients", "clip_advance_ratio"}},
    "wageningen_kt_kq": {"preprocess_wageningen": {"clip_advance_ratio"}},
    "fins": {"preprocess_fins": {"convention"}},
    "vsim_fins": {"preprocess_vsim_fins": set()},
    "outboard_motor": {"preprocess_outboard_motor_rpm": {"power_limit", "clip_advance_ratio"},
                       "preprocess_outboard_motor_throttle": set()},
}


@pytest.mark.parametrize("module_name", sorted(GENERIC))
def test_generic_block_reads_no_vehicle_name_and_has_no_vehicle_defaults(module_name):
    block = _contract(module_name)
    code = _code_without_docstrings(block)
    assert not [n for n in VEHICLE_NAMES if n in code]
    for function_name, allowed in GENERIC[module_name].items():
        signature = inspect.signature(getattr(block, function_name))
        defaults = {n for n, p in signature.parameters.items() if p.default is not inspect.Parameter.empty}
        assert defaults == allowed, (function_name, defaults)
