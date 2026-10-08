"""Gate tests for the force-producer blocks.

Written 2026-10-05, before the blocks existed. Gates (the test names carry
them): G1 against MATLAB running MSS (frozen CSVs), 1e-9 absolute; G2 against
an independent transcription of the cited MSS lines, 1e-10; G4 a perturbed model is detected; G5 a printed number reproduced.

The outboard motor's reverse thrust (owner's decision of 2026-10-05): full
reverse must push backwards (``test_outboard_block_full_reverse_pushes_backwards``).

Wageningen default (owner's decision of 2026-10-05; MSS is the reference):
MSS ``wageningen.m`` (no ``|J|``, no clamp to [0, 1.3]); the template clamp
stays reachable behind the flag ``clip_advance_ratio=True``; the default is
checked against MSS.

The differential thruster's default (otter.m ``g = 9.81``) is compared with
MATLAB running current MSS (2026-10-05,
``differential_thruster_mss_current.csv``); the template CSV (latitude gravity)
keeps the block's G1/G4 on the template settings.

Added 2026-10-06, section 7 at the end of this file:

* MATLAB R2026a running current MSS backs the fins and the propeller
  (``remus100_actuators_mss_current.csv``: ``remus100.m``'s own ``X_prop``,
  ``K_prop``, ``X_r``, ``X_s``, ``Y_r``, ``Z_s`` and ``tau`` on 1036 seeded
  cases) and the Wageningen polynomial (``wageningen_kt_kq_mss_current.csv``:
  ``wageningen.m`` on three geometries, J = -0.5 ... 2.0). Equality needs the
  parameter set to supply ``remus100.m``'s values: fin limit 20 deg, rho 1026,
  fin x = -a, and the printed K_T/K_Q (tests ``test_*_depart_from_matlab*``).
* Outboard: the reverse-sign fix is pinned; throttle mode gets its G4; the
  default path's reverse branch (signed J, reverse factors) gets a test.
* Transforms: no local rotation or cross product in
  ``models/force_producers``; rotations and skews come from
  ``more_transformations`` (numpy) and ``more_casadi_transformations``
  (graph).
* The allocation (thrust from the requested wrench, ``B_prop``) belongs to
  the ``more_control`` library and is tested there.

Producers and their references
------------------------------
=====================  ====================================================================
block                  reference in this repository
=====================  ====================================================================
differential thruster  MATLAB CSVs (MSS ``otter.m``; template) + ``otter.m`` transcription
propeller (B-series)   MATLAB CSVs (``remus100.m``; template) + ``remus100.m`` transcription
K_T/K_Q polynomial     MATLAB ``wageningen.m`` CSV + ``WageningData.mat`` / ``.txt``
fins (rudder, sterns)  MATLAB ``remus100.m`` CSV + transcription
VSIM fins              signature and physical tests only (no independent reference: gap)
outboard motor         MATLAB K_T/K_Q on its geometry + the force law written in the test
=====================  ====================================================================

Contract of the block (``more_dynamics.models.force_producers``)
----------------------------------------------------------------
All-CasADi blocks with named parameters. Every producer: a declaration
``<name>_parameters(...) -> tuple[Parameter, ...]`` (name, shape, SI unit,
meaning, range) and a builder ``<name>_casadi(*, <selectors>) ->
ca.Function`` with the command input(s) first, then ``"nu_r"`` (6x1), then
one input per declared parameter, and the output ``"tau"`` (6x1, N and N m,
BODY, about the CO) first. The tests check the numbers (``check_values`` or
the block's ``check_<name>_values``) and freeze them in (``freeze``, the
path a plugin carries), so every function below takes the commands and
``nu_r`` only. Only the named selectors have defaults.

* ``differential_thruster``: parameters ``positive_thrust_coefficients``,
  ``negative_thrust_coefficients`` (2x1, N s^2/rad^2), ``thruster_positions``,
  ``thruster_directions`` (2x3, rows ``[left; right]``; directions normalised
  inside, a zero row refused by ``check_differential_thruster_values``),
  ``max_forward_thrust`` / ``max_reverse_thrust`` (2x1 N, positive),
  giving ``n_max = sqrt(max_forward_thrust / k_pos)`` and ``n_min =
  -sqrt(max_reverse_thrust / k_neg)`` (``otter.m`` 136-137). Outputs
  ``tau``, ``allocation_matrix`` (6x2), ``max_speed``, ``min_speed`` (2x1
  rad/s). Inputs ``["n", "nu_r"]``, ``n`` (2x1) shaft speed in rad/s,
  saturated inside; thrust ``k_pos n|n|`` for ``n > 0``, else ``k_neg n|n|``.
* ``propeller``: ``propeller_casadi(*, open_water_model,
  open_water_coefficients="polynomial", clip_advance_ratio=False)``;
  parameters ``propeller_diameter``, ``max_shaft_speed`` (rpm), ``thrust_deduction``,
  ``wake_fraction``, ``pitch_diameter_ratio``, ``blade_area_ratio``,
  ``blade_count``, ``max_advance_number``, ``roll_moment_scale``,
  ``water_density``, ``position``, ``orientation`` (roll, pitch, yaw of the
  shaft, rad) and, with ``open_water_coefficients="given"``,
  ``thrust_torque_coefficients`` = ``(KT_0, KQ_0, KT_max, KQ_max)`` (MSS
  prints them in ``remus100.m`` 157-161; the tests pass ``None`` for the
  polynomial values). ``open_water_model`` in {"linearized", "full",
  "bollard"} (no default). Outputs ``tau``, ``shaft_axis``,
  ``coefficients_in_use`` (read here as ``thrust_torque_coefficients``).
  Inputs ``["n", "nu_r"]``, ``n`` (1x1) rpm; advance speed ``(1 - w)
  |nu_r[0:3]|``. The "full" model clips ``J`` to ``[0, 2
  max_advance_number]`` on both paths; ``clip_advance_ratio`` is passed to
  the Wageningen polynomial.
* ``wageningen_kt_kq``: ``wageningen_casadi(*, clip_advance_ratio=False)``,
  parameters ``pitch_diameter_ratio``, ``blade_area_ratio``,
  ``blade_count``; input ``["J"]`` (1x1), outputs ``["KT", "KQ"]`` (1x1
  each) — the one exception to the ``tau`` output. **Default:** MSS
  ``wageningen.m``, the polynomial at any ``J``, including ``J < 0`` and
  ``J > 1.3``. ``clip_advance_ratio=True``: the template clamp, ``J``
  replaced by ``min(|J|, 1.3)``.
* ``fins``: ``fins_casadi(*, convention="starboard_down_positive")``;
  parameters ``rudder_area``, ``stern_plane_area`` (m^2; MSS: ``2 * S_fin``),
  ``rudder_lift_coefficient``, ``stern_plane_lift_coefficient``,
  ``rudder_position``, ``stern_plane_position`` (x of the surface, m),
  ``max_deflection`` (rad), ``water_density``. Inputs ``["delta", "nu_r"]``,
  ``delta`` = (rudder, stern plane) rad, saturated inside. The default
  convention equals ``remus100.m`` 238-254.
* ``vsim_fins``: ``vsim_fins_casadi(*, fin_count)``; parameters
  ``max_forces``, ``positions`` (fin_count x 3), ``max_deflection``. Inputs
  ``["delta", "nu_r"]``, ``delta`` (fin_count x 1).
* ``outboard_motor``: two blocks. ``outboard_motor_rpm_casadi(*,
  power_limit=True, clip_advance_ratio=False)`` with parameters
  ``max_thrust``, ``max_power``, ``efficiency``, ``position``, ``max_speed``,
  ``propeller_diameter``, ``pitch_diameter_ratio``, ``blade_area_ratio``,
  ``blade_count``, ``thrust_deduction``, ``water_density``,
  ``advance_speed_factor``, ``propwash_factor``, ``reverse_thrust_factor``,
  ``reverse_torque_factor`` (signed ``J``, so backing down with forward flow
  gives ``J < 0`` on the default path; ranges: ``efficiency`` in ``(0, 1]``,
  ``thrust_deduction`` and ``propwash_factor`` in ``[0, 1]``, the reverse
  factors ``> 0``, refused outside by ``check_values``) and inputs ``["n", "delta", "nu_r"]``:
  ``n`` the actual propeller speed (rpm, saturated inside), ``delta`` the
  actual steering angle (rad, port-positive). Rpm lag, steering rate limit
  and battery are states and are not part of this force map. **Reverse
  fixed:** for ``n < 0`` the reverse factors scale ``K_T`` and ``K_Q``
  without the extra ``sign(n)``, so full reverse gives a negative surge
  force. ``outboard_motor_throttle_casadi()`` with parameters
  ``max_thrust``, ``max_propulsive_power``, ``position`` and inputs
  ``["thrust", "delta", "nu_r"]`` (actual thrust state in N).

Frozen reference: ``tests/data/force_producers/`` (``SOURCE.md``).
"""

import ast
import importlib
import inspect
import json
import os
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "force_producers"
# Nothing relative to one machine (owner, 2026-10-06).
# Files outside this repo are found through environment variables, default
# unset: MSS_DIR (the MSS checkout root). A cited path below starts with the
# key of EXTERNAL_ROOTS (the rest is relative to that root) or with "tests/"
# (a file of this repository).
EXTERNAL_ROOTS = {"source-sim/MSS/": "MSS_DIR"}
TESTS_ROOT = Path(__file__).resolve().parents[1]   # "tests/..." paths: this repository
PACKAGE = "more_dynamics.models.force_producers"

G1_TOLERANCE = 1e-9   # G1: block vs MATLAB running MSS, absolute
G2_TOLERANCE = 1e-10  # G2: block vs a transcription, absolute
G4_FACTOR = 10.0      # G4: a perturbed model must differ by more than this x G1
SEED = 20261005
N_RANDOM_STATES = 1000

VEHICLE_NAMES = ("remus", "otter", "grethe", "marie", "hugin", "lauv",
                 "mariner", "torqeedo", "cybership", "prestero")

REMUS = "source-sim/MSS/CRAFT/AUV/models/remus100.m"
OTTER = "source-sim/MSS/CRAFT/USV/models/otter.m"
SIMOTTER = "source-sim/MSS/CRAFT/USV/SIMotter.m"
GRAVITY = "source-sim/MSS/INS/functions/gravity.m"
WAG_MAT = "source-sim/MSS/LIBRARY/modeling/utiles/WageningData.mat"
WAG_TXT = "source-sim/MSS/LIBRARY/modeling/utiles/WageningData.txt"
WAG_M = "source-sim/MSS/LIBRARY/modeling/wageningen.m"
# The template generators of the two template references, in this repository
# ("tests/..." = this tests/ folder; always present, no variable needed).
GEN_PROP = "tests/data/force_producers/generate_propeller_template.m"
GEN_DIFF = "tests/data/force_producers/generate_differential_thruster_template.m"

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
    (REMUS, 98): "rho = 1026;",
    (REMUS, 109): "delta_max = deg2rad(20);",
    (REMUS, 110): "n_max = 1525;",
    (REMUS, 113): "delta_r = sat(ui(1), delta_max);",
    (REMUS, 114): "delta_s = sat(ui(2), delta_max);",
    (REMUS, 115): "n_p = sat(ui(3),n_max) / 60;",
    (REMUS, 127): "U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );",
    (REMUS, 131): "L_auv = 1.6;",
    (REMUS, 134): "a = 1.0096 * L_auv/2;",
    (REMUS, 148): "D_prop = 0.14;",
    (REMUS, 149): "t_prop = 0.1;",
    (REMUS, 150): "Va = 0.944 * U_r;",
    (REMUS, 153): "Ja_max = 0.6632;",
    (REMUS, 156): "% >> [KT_0, KQ_0] = wageningen(0,1,0.718,3)",
    (REMUS, 157): "KT_0 = 0.4566;",
    (REMUS, 158): "KQ_0 = 0.0700;",
    (REMUS, 159): "% >> [KT_max, KQ_max] = wageningen(0.6632,1,0.718,3)",
    (REMUS, 160): "KT_max = 0.1798;",
    (REMUS, 161): "KQ_max = 0.0312;",
    (REMUS, 167): "if n_p > 0",
    (REMUS, 169): "X_prop = rho * D_prop^4 * ( ...",
    (REMUS, 170): "KT_0 * abs(n_p) * n_p + (KT_max-KT_0)/Ja_max * (Va/D_prop) * abs(n_p) );",
    (REMUS, 171): "K_prop = rho * D_prop^5 * ( ...",
    (REMUS, 172): "KQ_0 * abs(n_p) * n_p + (KQ_max-KQ_0)/Ja_max * (Va/D_prop) * abs(n_p) );",
    (REMUS, 175): "X_prop = rho * D_prop^4 * KT_0 * abs(n_p) * n_p;",
    (REMUS, 176): "K_prop = rho * D_prop^5 * KQ_0 * abs(n_p) * n_p;",
    (REMUS, 179): "S_fin = 0.00665;",
    (REMUS, 182): "CL_delta_r = 0.5;",
    (REMUS, 183): "A_r = 2 * S_fin;",
    (REMUS, 184): "x_r = -a;",
    (REMUS, 187): "CL_delta_s = 0.7;",
    (REMUS, 188): "A_s = 2 * S_fin;",
    (REMUS, 189): "x_s = -a;",
    (REMUS, 234): "U_rh = sqrt( nu_r(1)^2 + nu_r(2)^2 );",
    (REMUS, 235): "U_rv = sqrt( nu_r(1)^2 + nu_r(3)^2 );",
    (REMUS, 238): "X_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r^2;",
    (REMUS, 239): "X_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s^2;",
    (REMUS, 242): "Y_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r;",
    (REMUS, 245): "Z_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s;",
    (REMUS, 249): "tau(1) = (1-t_prop) * X_prop + X_r + X_s;",
    (REMUS, 250): "tau(2) = Y_r;",
    (REMUS, 251): "tau(3) = Z_s;",
    (REMUS, 252): "tau(4) = K_prop / 10;",
    (REMUS, 253): "tau(5) = -x_s * Z_s;",
    (REMUS, 254): "tau(6) = x_r * Y_r;",
    (GRAVITY, 11): "g = 9.7803253359 * ( 1 + 0.001931850400 * sin(mu)^2 ) /...",
    (GRAVITY, 12): "sqrt( 1 - 0.006694384442 * sin(mu)^2 );",
    (GEN_DIFF, 25): "g      = gravity(deg2rad(63.446827));",
    (GEN_PROP, 97): "rho = 1025;",
    (GEN_PROP, 98): "n_max = 1525;",
    (GEN_PROP, 106): "D_prop = 0.14;",
    (GEN_PROP, 107): "t_prop = 0.1;",
    (GEN_PROP, 108): "Va = 0.944 * U;",
    (GEN_PROP, 110): "Ja_max = 0.6632;",
    (GEN_PROP, 112): "KT_0 = 0.4566;",
    (GEN_PROP, 113): "KQ_0 = 0.0700;",
    (GEN_PROP, 115): "KT_max = 0.1798;",
    (GEN_PROP, 116): "KQ_max = 0.0312;",
    (GEN_PROP, 144): "K_prop = K_prop / 10;",
    (WAG_M, 33): "load('WageningData.mat');",
    (WAG_M, 35): "KT = sum(WagCThrust_stuv.*((Ja).^WagThrust_s).*(PD.^WagThrust_t).*...",
    (WAG_M, 37): "KQ = sum(WagCTorque_stuv.*((Ja).^WagTorque_s).*(PD.^WagTorque_t).*...",
    # added 2026-10-06
    (OTTER, 212): "B_prop = k_pos * [...",
    (OTTER, 213): "1 1",
    (OTTER, 214): "y_pont -y_pont ];",
    (SIMOTTER, 182): "u = Binv * [tau_X; tau_N];",
    (SIMOTTER, 183): "n_c = sign(u) .* sqrt(abs(u));",
    (WAG_M, 24): "%   Barnitsas, M.M., Ray, D. and Kinley, P. (1981).",
    (WAG_M, 25): "%   KT, KQ and Efficiency Curves for the Wageningen B-Series Propellers",
    (WAG_M, 26): "%   http://deepblue.lib.umich.edu/handle/2027.42/3557",
}


# --------------------------------------------------------------------------
# Helpers: pinned text, contract, data
# --------------------------------------------------------------------------
def _external(rel):
    """Path of a cited file outside this repo, or (None, reason) when the
    environment variable for its root is unset or the file is missing."""
    if rel.startswith("tests/"):
        return TESTS_ROOT / rel[len("tests/"):], None
    for prefix, variable in EXTERNAL_ROOTS.items():
        if rel.startswith(prefix):
            root = os.environ.get(variable)
            if not root:
                return None, f"{variable} is not set (nothing relative to one machine); needed for {rel}"
            path = Path(root) / rel[len(prefix):]
            return (path, None) if path.exists() else (None, f"{variable}: cited file not found: {path}")
    raise KeyError(rel)


def _external_or_skip(rel):
    path, reason = _external(rel)
    if path is None:
        pytest.skip(reason)
    return path


def _line(rel, number):
    """The cited line, checked against its pinned text. Without the external
    checkout the pinned text itself is used (frozen in CITED_LINES), so the
    gates still run; only the unchanged-line check is skipped."""
    path, reason = _external(rel)
    if path is None:
        if (rel, number) in CITED_LINES:
            return CITED_LINES[(rel, number)]
        pytest.skip(reason)
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


# modules under PACKAGE: the thruster and propulsor blocks sit in their own sub-folders
MODULE_PATH = {"differential_thruster": "thrusters.differential_thruster", "outboard_motor": "thrusters.outboard_motor",
               "propeller": "propulsor.propeller", "wageningen_kt_kq": "propulsor.wageningen_kt_kq"}


def _contract(name):
    try:
        return importlib.import_module(f"{PACKAGE}.{MODULE_PATH.get(name, name)}")
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed: {exc}")
        pytest.fail(f"block not ported yet: {exc}")


MODULE_OF = {"outboard_motor_rpm": "outboard_motor", "outboard_motor_throttle": "outboard_motor",
             "wageningen": "wageningen_kt_kq"}


# selectors of each builder (the rest of a parameter set is numbers)
SELECTORS = {
    "differential_thruster": (), "propeller": ("open_water_model", "clip_advance_ratio"),
    "wageningen": ("clip_advance_ratio",), "fins": ("convention",), "vsim_fins": (),
    "outboard_motor_rpm": ("power_limit", "clip_advance_ratio"), "outboard_motor_throttle": (),
}
OUTPUT_FIELDS = {"coefficients_in_use": "thrust_torque_coefficients"}   # output -> constants attribute


def _as_field(value):
    array = np.array(value, dtype=float)
    if array.size == 1:
        return float(array.ravel()[0])
    return array.ravel() if 1 in array.shape else array


def _build(name, params):
    """The block ``<name>_casadi(<selectors>)`` with its numbers checked and
    frozen in (``freeze``): a function of the commands and ``nu_r``; and a
    namespace with the checked numbers and the named outputs at zero
    commands (the fields the assertions read)."""
    from more_transformations.more_casadi_transformations import check_values, freeze

    block = _contract(MODULE_OF.get(name, name))
    numbers = dict(params)
    selectors = {k: numbers.pop(k) for k in SELECTORS[name] if k in numbers}
    declaration = {}
    if name == "propeller":
        given = numbers.pop("thrust_torque_coefficients", None)
        declaration["open_water_coefficients"] = "polynomial" if given is None else "given"
        if given is not None:
            numbers["thrust_torque_coefficients"] = given
    if name == "vsim_fins":
        declaration["fin_count"] = len(numbers["max_forces"])
    declared = getattr(block, f"{name}_parameters")(**declaration)
    checker = getattr(block, f"check_{name}_values", None)
    values = checker(numbers) if checker else check_values(declared, numbers)
    function = freeze(getattr(block, f"{name}_casadi")(**selectors, **declaration), declared, values)
    rest = function(**{n: np.zeros(function.size_in(n)) for n in function.name_in()})
    fields = {k: _as_field(v) for k, v in values.items()}
    fields.update({OUTPUT_FIELDS.get(k, k): _as_field(v) for k, v in rest.items() if k not in ("tau", "KT", "KQ")})
    return SimpleNamespace(**fields), function


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


# --------------------------------------------------------------------------
# Pinned lines (no port needed)
# --------------------------------------------------------------------------
def test_cited_lines_are_unchanged():
    for (rel, number), text in CITED_LINES.items():
        path, reason = _external(rel)
        if path is None:
            pytest.skip(reason)
        line = path.read_text().splitlines()[number - 1].strip()
        assert line.startswith(text), (rel, number, line)


# ==========================================================================
# 1. Differential thruster (Otter-type twin propellers)
# ==========================================================================
def _otter_differential(gravity):
    """otter.m 104, 131-136; positions [left, right] = [0, -/+ y_pont, 0]."""
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
    """otter.m 219-231 (left = n(1), lever arm l1)."""
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


def _matlab_gravity_differential():
    return _gravity(np.deg2rad(_number_in(GEN_DIFF, 25, r"deg2rad\(([0-9.]+)\)")))


def test_G1_block_differential_matches_matlab():
    _, function = _build("differential_thruster", _otter_differential(_matlab_gravity_differential()))
    ref = _differential_csv()
    for k, n in enumerate(ref["n"]):
        tau = _tau(function, n=n, nu_r=np.zeros(6))
        assert _max_diff(tau[[0, 5]], [ref["tau_X"][k], ref["tau_N"][k]]) <= G1_TOLERANCE, k
        assert _max_diff(tau[1:5], 0.0) <= G1_TOLERANCE, k


def test_G2_differential_signature():
    constants, function = _build("differential_thruster", _otter_differential(9.81))
    assert function.name_in() == ["n", "nu_r"] and function.name_out() == ["tau", "allocation_matrix", "max_speed", "min_speed"]
    assert function.size_in(0) == (2, 1) and function.size_in(1) == (6, 1)
    assert function.size_out(0) == (6, 1)
    for field in ("max_speed", "min_speed", "allocation_matrix"):
        assert hasattr(constants, field), field


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
    """otter.m's own g = 9.81 (line 89) vs MATLAB running otter.m (2026-10-05): the
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
    """otter.m 131-136, 219-231 with otter.m's own g = 9.81 (line 89)."""
    p = _otter_differential(_value(OTTER, 90))
    _, function = _build("differential_thruster", p)
    for k, n in enumerate(_rng().uniform(-150.0, 150.0, size=(N_RANDOM_STATES, 2))):
        assert _max_diff(_tau(function, n=n, nu_r=np.zeros(6)), _mss_otter_tau(p, n)) <= G1_TOLERANCE, (
            "otter.m 219-231", k)


# ==========================================================================
# 2. Propeller (Wageningen B-series, single screw)
# ==========================================================================
OPEN_WATER_MODES = {"linearized": "linearized_open_water", "full": "full_open_water",
                    "bollard": "bollard_pull"}


def _matlab_propeller():
    """test_mss_reference.m 83-102, 130: remus100.m's propeller with rho = 1025."""
    return {
        "propeller_diameter": _value(GEN_PROP, 106), "max_shaft_speed": _value(GEN_PROP, 98),
        "thrust_deduction": _value(GEN_PROP, 107),
        "wake_fraction": 1.0 - _number_in(GEN_PROP, 108, r"= ([0-9.]+) \* U"),
        "pitch_diameter_ratio": 1.0, "blade_area_ratio": 0.718, "blade_count": 3,  # remus100.m 156
        "max_advance_number": _value(GEN_PROP, 110),
        "roll_moment_scale": 1.0 / _number_in(GEN_PROP, 144, r"/ ([0-9.]+);"),
        "water_density": _value(GEN_PROP, 97),
        "position": [0.0, 0.0, 0.0], "orientation": [0.0, 0.0, 0.0],
        "open_water_model": "linearized",
        "thrust_torque_coefficients": [_value(GEN_PROP, n) for n in (112, 113, 115, 116)],
    }


def _remus_propeller():
    """remus100.m 98, 110, 148-161, 252 (rho = 1026); geometry of line 156."""
    pd, aeao, z = (float(t) for t in re.search(
        r"wageningen\(0,([0-9.]+),([0-9.]+),([0-9]+)\)", _line(REMUS, 156)).groups())
    return {
        "propeller_diameter": _value(REMUS, 148), "max_shaft_speed": _value(REMUS, 110),
        "thrust_deduction": _value(REMUS, 149),
        "wake_fraction": 1.0 - _number_in(REMUS, 150, r"= ([0-9.]+) \* U_r"),
        "pitch_diameter_ratio": pd, "blade_area_ratio": aeao, "blade_count": int(z),
        "max_advance_number": _value(REMUS, 153),
        "roll_moment_scale": 1.0 / _number_in(REMUS, 252, r"K_prop / ([0-9.]+);"),
        "water_density": _value(REMUS, 98),
        "position": [0.0, 0.0, 0.0], "orientation": [0.0, 0.0, 0.0],
        "open_water_model": "linearized",
        "thrust_torque_coefficients": [_value(REMUS, n) for n in (157, 158, 160, 161)],
    }


def _scaled_propeller():
    """Second vehicle (test construction): a bigger 4-blade propeller mounted
    aft, below and slightly tilted."""
    p = _remus_propeller()
    return {**p, "propeller_diameter": 0.3, "max_shaft_speed": 800.0, "thrust_deduction": 0.15,
            "wake_fraction": 0.1, "pitch_diameter_ratio": 0.9, "blade_area_ratio": 0.55,
            "blade_count": 4, "max_advance_number": 0.8, "roll_moment_scale": 0.2,
            "water_density": 1025.0, "position": [-0.8, 0.0, 0.05],
            "orientation": [0.0, 0.1, 0.05], "thrust_torque_coefficients": None}


def _mss_remus_propeller(p, rpm, nu_r):
    """remus100.m 115, 127, 150, 167-177, 249 (propeller part), 250."""
    kt0, kq0, ktm, kqm = p["thrust_torque_coefficients"]
    n_p = np.clip(rpm, -p["max_shaft_speed"], p["max_shaft_speed"]) / 60.0
    va = (1.0 - p["wake_fraction"]) * np.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2)
    rho, d, j_max = p["water_density"], p["propeller_diameter"], p["max_advance_number"]
    if n_p > 0:
        x = rho * d ** 4 * (kt0 * abs(n_p) * n_p + (ktm - kt0) / j_max * (va / d) * abs(n_p))
        k = rho * d ** 5 * (kq0 * abs(n_p) * n_p + (kqm - kq0) / j_max * (va / d) * abs(n_p))
    else:
        x = rho * d ** 4 * kt0 * abs(n_p) * n_p
        k = rho * d ** 5 * kq0 * abs(n_p) * n_p
    return np.array([(1 - p["thrust_deduction"]) * x, 0, 0, p["roll_moment_scale"] * k, 0, 0])


def _propeller_csv():
    h, v = _load_csv("propeller_matlab.csv")
    return {"rpm": v[:, h.index("RPM")], "U": v[:, h.index("U")],
            "tau": v[:, [h.index(f"tau{i}") for i in range(1, 7)]]}


def test_propeller_reference_is_the_linearised_remus_model():
    """The stored CSV is remus100.m 167-177 with test_mss_reference.m's rho."""
    p, ref = _matlab_propeller(), _propeller_csv()
    for k in range(len(ref["rpm"])):
        expected = _mss_remus_propeller(p, ref["rpm"][k], [ref["U"][k], 0, 0, 0, 0, 0])
        assert _max_diff(expected, ref["tau"][k]) <= G1_TOLERANCE, k


def test_G1_block_propeller_matches_matlab():
    _, function = _build("propeller", _matlab_propeller())
    ref = _propeller_csv()
    for k in range(len(ref["rpm"])):
        tau = _tau(function, n=[ref["rpm"][k]], nu_r=[ref["U"][k], 0, 0, 0, 0, 0])
        assert _max_diff(tau, ref["tau"][k]) <= G1_TOLERANCE, k


def test_G2_propeller_signature():
    constants, function = _build("propeller", _remus_propeller())
    assert function.name_in() == ["n", "nu_r"] and function.name_out() == ["tau", "shaft_axis", "coefficients_in_use"]
    assert function.size_in(0) == (1, 1) and function.size_in(1) == (6, 1)
    assert function.size_out(0) == (6, 1)
    np.testing.assert_allclose(constants.thrust_torque_coefficients,
                               _remus_propeller()["thrust_torque_coefficients"], atol=0, rtol=0)


def test_MSS_propeller_full_model_default_beyond_J_1_3_uses_unclamped_polynomial():
    """Default = MSS (owner, 2026-10-05). Test construction: the scaled propeller on the shaft
    axis at the CO (position and orientation zero), "full" model, J = 1.5,
    inside the "full" model's clip of J to 2 x 0.8 = 1.6 and beyond 1.3, the end
    of the regression's fitted range. Default: K_T, K_Q = MSS wageningen.m at
    J; the ``clip_advance_ratio`` flag changes the force."""
    p = {**_scaled_propeller(), "position": [0.0, 0.0, 0.0], "orientation": [0.0, 0.0, 0.0],
         "open_water_model": "full", "thrust_torque_coefficients": None}
    j, rpm = 1.5, 400.0
    assert 1.3 < j < 2.0 * p["max_advance_number"]
    n_rps = rpm / 60.0
    nu_r = np.array([j * n_rps * p["propeller_diameter"] / (1.0 - p["wake_fraction"]), 0, 0, 0, 0, 0])
    kt, kq = _mss_wageningen(j, p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"])
    rho, d = p["water_density"], p["propeller_diameter"]
    expected = np.array([(1.0 - p["thrust_deduction"]) * rho * d ** 4 * kt * n_rps ** 2, 0, 0,
                         p["roll_moment_scale"] * rho * d ** 5 * kq * n_rps ** 2, 0, 0])
    default = _tau(_build("propeller", p)[1], n=[rpm], nu_r=nu_r)
    assert _max_diff(default, expected) <= G1_TOLERANCE, ("wageningen.m at J = 1.5", default, expected)
    flagged = _tau(_build("propeller", {**p, "clip_advance_ratio": True})[1], n=[rpm], nu_r=nu_r)
    assert _max_diff(default, flagged) > G4_FACTOR * G1_TOLERANCE


def test_G4_propeller_diameter_plus_1_percent_is_detected():
    p, ref = _matlab_propeller(), _propeller_csv()

    def worst(params):
        _, f = _build("propeller", params)
        return max(_max_diff(_tau(f, n=[ref["rpm"][k]], nu_r=[ref["U"][k], 0, 0, 0, 0, 0]), ref["tau"][k])
                   for k in range(len(ref["rpm"])))

    assert worst(p) <= G1_TOLERANCE
    assert worst({**p, "propeller_diameter": 1.01 * p["propeller_diameter"]}) > G4_FACTOR * G1_TOLERANCE


def test_MSS_propeller_equals_remus100():
    """remus100.m 98 (rho 1026), 111, 116, 149-178, 247, 250 with the printed
    coefficients (158-162) as the exact override."""
    p = _remus_propeller()
    _, function = _build("propeller", p)
    rng = _rng()
    for k in range(N_RANDOM_STATES):
        rpm, nu_r = rng.uniform(-2000.0, 2000.0), rng.uniform(-3.0, 3.0, 6)
        expected = _mss_remus_propeller(p, rpm, nu_r)
        assert _max_diff(_tau(function, n=[rpm], nu_r=nu_r), expected) <= G1_TOLERANCE, ("remus100.m", k)


def test_G5_polynomial_coefficients_round_to_remus100_printed_values():
    """remus100.m 156-161 print KT/KQ at J = 0 and J = Ja_max from
    ``wageningen(J,1,0.718,3)`` to 4 decimals; the polynomial must land on them
    (the WageningData.mat transcription, itself G1 against MATLAB)."""
    wageningen = _mss_wageningen
    p = _remus_propeller()
    printed = p["thrust_torque_coefficients"]
    kt0, kq0 = wageningen(0.0, p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"])
    ktm, kqm = wageningen(p["max_advance_number"], p["pitch_diameter_ratio"],
                          p["blade_area_ratio"], p["blade_count"])
    np.testing.assert_array_equal(np.round([kt0, kq0, ktm, kqm], 4), printed,
                                  err_msg="remus100.m lines 156-161 (printed precision)")


def test_G5_block_polynomial_coefficients_round_to_remus100_printed_values():
    p = {**_remus_propeller(), "thrust_torque_coefficients": None}
    constants, _ = _build("propeller", p)
    np.testing.assert_array_equal(np.round(constants.thrust_torque_coefficients, 4),
                                  _remus_propeller()["thrust_torque_coefficients"],
                                  err_msg="remus100.m lines 156-161 (printed precision)")


# ==========================================================================
# 3. Wageningen K_T / K_Q polynomial
# ==========================================================================
def _mss_wageningen_tables():
    """WageningData.mat, loaded as wageningen.m line 33 does."""
    import scipy.io
    path = _external_or_skip(WAG_MAT)
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
    path = _external_or_skip(WAG_TXT)
    text = path.read_text()
    rows = lambda s: np.array([[float(x) for x in ln.split()] for ln in s.splitlines()
                               if re.match(r"\s*[-+][0-9]", ln)])
    thrust_txt, torque_txt = (rows(part) for part in text.split("TORQUE:"))
    thrust, torque = _mss_wageningen_tables()
    assert thrust.shape == (39, 5) and torque.shape == (47, 5)
    np.testing.assert_array_equal(thrust, thrust_txt)
    np.testing.assert_array_equal(torque, torque_txt)


def _build_wageningen(pd, aeao, z, **flag):
    return _build("wageningen", {"pitch_diameter_ratio": pd, "blade_area_ratio": aeao, "blade_count": z,
                                 **flag})


def _wageningen_grid():
    """J on -0.5 ... 2.0, past both ends of the source's [0, 1.3]."""
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
    """The default is wageningen.m 35-38 at J itself, also
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
    """remus100.m 98, 109, 134, 179-189."""
    a = _value(REMUS, 134, {"L_auv": _value(REMUS, 131)})
    s_fin = _value(REMUS, 179)
    return {"rudder_area": _value(REMUS, 183, {"S_fin": s_fin}),
            "stern_plane_area": _value(REMUS, 188, {"S_fin": s_fin}),
            "rudder_lift_coefficient": _value(REMUS, 182),
            "stern_plane_lift_coefficient": _value(REMUS, 187),
            "rudder_position": _value(REMUS, 184, {"a": a}),
            "stern_plane_position": _value(REMUS, 189, {"a": a}),
            "max_deflection": _value(REMUS, 109),
            "water_density": _value(REMUS, 98)}


def _mss_fins(p, delta, nu_r):
    """remus100.m 113-114, 234-254 (fin part of tau)."""
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
            "remus100.m 113-114, 234-254", k)


# ==========================================================================
# 5. VSIM fins (signature and physical tests)
# ==========================================================================
def _vsim_parameters(scaled):
    p = _parameter_set("vsim_fins")
    if scaled:
        p = {"max_forces": 2.5 * p["max_forces"], "positions": 1.7 * p["positions"],
             "max_deflection": 0.35}
    return p


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


# ==========================================================================
# 6. Outboard motor (reverse sign fixed, owner 2026-10-05)
# ==========================================================================
def _parameter_set(name):
    """A frozen parameter set of tests/data/force_producers/parameter_sets.json (SOURCE.md)."""
    p = json.loads((DATA_DIR / "parameter_sets.json").read_text())[name]
    return {k: (np.asarray(v, float) if isinstance(v, list) else v) for k, v in p.items()}


def _outboard_rpm_parameters(scaled):
    """The outboard parameter set (``parameter_sets.json``); ``scaled`` is a test construction."""
    p = _parameter_set("outboard_rpm")
    if scaled:
        p = {**p, "max_thrust": 1500.0, "max_power": 6000.0, "position": np.array([-3.5, 0.4, 0.7]),
             "max_speed": 2000.0, "propeller_diameter": 0.4, "pitch_diameter_ratio": 1.1,
             "blade_count": 4, "reverse_thrust_factor": 0.4, "reverse_torque_factor": 0.7}
    return p


def _outboard_cases(p):
    rng = _rng()
    rpm = rng.uniform(-1.2 * p["max_speed"], 1.2 * p["max_speed"], N_RANDOM_STATES)
    delta = rng.uniform(-0.6, 0.6, N_RANDOM_STATES)   # within the source's 35 deg
    nu_r = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    return rpm, delta, nu_r


def test_outboard_block_full_reverse_pushes_backwards():
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


def _outboard_rpm_force(p, rpm, delta, nu_r, kt, kq, power_limit=True):
    """The outboard force law for given K_T, K_Q (n > 0, no reverse factors):
    T = (1 - t) rho D^4 K_T |n| n, Q = rho D^5 K_Q |n| n, T scaled down when
    the shaft power |2 pi n Q| exceeds P_max * efficiency, capped at T_max,
    steered by delta: F = T [cos delta, sin delta, 0], M = r x F."""
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
    """Default = MSS (owner, 2026-10-05). Test construction: J = 1.6 at 300 rpm (surge speed
    from J = Va_factor (1 - w_prop) u / (n D), lines 315 and 382), steering
    0.1 rad. Default: K_T, K_Q = MSS wageningen.m at J; the ``clip_advance_ratio``
    flag changes the force."""
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
    assert _max_diff(default, flagged) > G4_FACTOR * G1_TOLERANCE


# Physical ranges of the outboard's dimensionless numbers (refused by check_values, naming them)
OUTBOARD_RANGE_REFUSALS = [
    ("efficiency", -1.0, "parameter 'efficiency' [1] must be > 0 and <= 1 in every entry, got -1.0"),
    ("efficiency", 0.0, "parameter 'efficiency' [1] must be > 0 and <= 1 in every entry, got 0.0"),
    ("efficiency", 1.5, "parameter 'efficiency' [1] must be > 0 and <= 1 in every entry, got 1.5"),
    ("thrust_deduction", 2.0, "parameter 'thrust_deduction' [1] must be >= 0 and <= 1 in every entry, got 2.0"),
    ("thrust_deduction", -0.1, "parameter 'thrust_deduction' [1] must be >= 0 and <= 1 in every entry, got -0.1"),
    ("propwash_factor", 2.0, "parameter 'propwash_factor' [1] must be >= 0 and <= 1 in every entry, got 2.0"),
    ("propwash_factor", -0.1, "parameter 'propwash_factor' [1] must be >= 0 and <= 1 in every entry, got -0.1"),
    ("reverse_thrust_factor", -1.0, "parameter 'reverse_thrust_factor' [1] must be > 0 in every entry, got -1.0"),
    ("reverse_thrust_factor", 0.0, "parameter 'reverse_thrust_factor' [1] must be > 0 in every entry, got 0.0"),
    ("reverse_torque_factor", -1.0, "parameter 'reverse_torque_factor' [1] must be > 0 in every entry, got -1.0"),
    ("reverse_torque_factor", 0.0, "parameter 'reverse_torque_factor' [1] must be > 0 in every entry, got 0.0"),
]


@pytest.mark.parametrize("name, value, message", OUTBOARD_RANGE_REFUSALS,
                         ids=[f"{n}={v}" for n, v, _ in OUTBOARD_RANGE_REFUSALS])
def test_outboard_rpm_ratio_outside_its_physical_range_is_refused(name, value, message):
    """efficiency in (0, 1], thrust_deduction and propwash_factor in [0, 1],
    the reverse factors > 0 (module docstring): a value outside is refused
    with the parameter's name, unit and range."""
    from more_transformations.more_casadi_transformations import check_values

    block = _contract("outboard_motor")
    p = {k: v for k, v in _outboard_rpm_parameters(False).items() if k not in SELECTORS["outboard_motor_rpm"]}
    with pytest.raises(ValueError) as error:
        check_values(block.outboard_motor_rpm_parameters(), {**p, name: value})
    assert str(error.value) == message


@pytest.mark.parametrize("name, value", [("efficiency", 1.0), ("efficiency", 1e-4), ("thrust_deduction", 0.0),
                                         ("thrust_deduction", 1.0), ("propwash_factor", 0.0),
                                         ("propwash_factor", 1.0), ("reverse_thrust_factor", 2.0),
                                         ("reverse_torque_factor", 2.0)])
def test_outboard_rpm_ratio_at_the_ends_of_its_range_is_accepted(name, value):
    """The closed ends (efficiency 1, t and w 0 and 1) and values inside
    build a finite force; efficiency 1e-4 reaches the efficiency floor of the
    power cap (module docstring), the force stays finite."""
    p = _outboard_rpm_parameters(False)
    _, function = _build("outboard_motor_rpm", {**p, name: value})
    for rpm in (-p["max_speed"], 0.5 * p["max_speed"], p["max_speed"]):
        assert np.all(np.isfinite(_tau(function, n=[rpm], delta=[0.1], nu_r=np.array([1.0, 0.1, 0, 0, 0, 0.1]))))


def test_propeller_fractions_have_the_outboard_range():
    """A name is one physical quantity: the propeller's thrust_deduction and
    wake_fraction carry the [0, 1] range of the outboard's thrust_deduction
    and propwash_factor; max_advance_number > 0 and roll_moment_scale state
    their provenance."""
    propeller = _contract(MODULE_OF.get("propeller", "propeller"))
    outboard = _contract("outboard_motor")
    by_name = {d.name: d for d in propeller.propeller_parameters()}
    rpm = {d.name: d for d in outboard.outboard_motor_rpm_parameters()}
    for name, twin in (("thrust_deduction", "thrust_deduction"), ("wake_fraction", "propwash_factor")):
        a, b = by_name[name], rpm[twin]
        assert (a.minimum, a.maximum, a.minimum_exclusive, a.maximum_exclusive) == \
               (b.minimum, b.maximum, b.minimum_exclusive, b.maximum_exclusive) == (0.0, 1.0, False, False)
    j_max = by_name["max_advance_number"]
    assert (j_max.minimum, j_max.minimum_exclusive, j_max.maximum) == (0.0, True, None)
    assert "not a physics value" in by_name["roll_moment_scale"].meaning
    vsim = _contract(MODULE_OF.get("vsim_fins", "vsim_fins"))
    assert "not a physics value" in vsim.vsim_fins_parameters(fin_count=2)[0].meaning


def _outboard_throttle_parameters(scaled):
    p = _parameter_set("outboard_throttle")
    if scaled:
        p = {"max_thrust": 1500.0, "max_propulsive_power": 4000.0, "position": np.array([-3.5, 0.4, 0.7])}
    return p


# ==========================================================================
# Generic by construction (no vehicle name or number)
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


GENERIC = {  # module -> {builder or declaration -> selectors allowed a default}
    "differential_thruster": {"differential_thruster_casadi": set(), "differential_thruster_parameters": set()},
    "propeller": {"propeller_casadi": {"open_water_coefficients", "clip_advance_ratio"},
                  "propeller_parameters": {"open_water_coefficients"}},
    "wageningen_kt_kq": {"wageningen_casadi": {"clip_advance_ratio"}, "wageningen_parameters": set()},
    "fins": {"fins_casadi": {"convention"}, "fins_parameters": set()},
    "vsim_fins": {"vsim_fins_casadi": set(), "vsim_fins_parameters": set()},
    "outboard_motor": {"outboard_motor_rpm_casadi": {"power_limit", "clip_advance_ratio"},
                       "outboard_motor_throttle_casadi": set(), "outboard_motor_rpm_parameters": set(),
                       "outboard_motor_throttle_parameters": set()},
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


# ==========================================================================
# 7. Added 2026-10-06: current-MSS MATLAB for fins, propeller and the
#    Wageningen polynomial; outboard pins; transform imports
# ==========================================================================
def _actuators_csv():
    """remus100_actuators_mss_current.csv: MATLAB R2026a running remus100.m
    (MSS ac77394) on 1036 seeded cases; columns are remus100.m's own
    workspace variables (SOURCE.md, 2026-10-06)."""
    h, v = _load_csv("remus100_actuators_mss_current.csv")
    get = lambda *names: v[:, [h.index(n) for n in names]]
    one = lambda name: v[:, h.index(name)]
    return {"ui": get("ui1", "ui2", "ui3"),
            "nu_r": get(*[f"nu_r_{i:02d}" for i in range(1, 7)]),
            "n_p": one("n_p"), "delta": get("delta_r", "delta_s"),
            "X_prop": one("X_prop"), "K_prop": one("K_prop"), "X_r": one("X_r"),
            "X_s": one("X_s"), "Y_r": one("Y_r"), "Z_s": one("Z_s"),
            "tau": get(*[f"tau_{i:02d}" for i in range(1, 7)])}


def _matlab_propeller_part(ref):
    """remus100.m 249 (propeller term) and 250, from the MATLAB columns."""
    t_prop = _value(REMUS, 149)
    scale = _number_in(REMUS, 252, r"K_prop / ([0-9.]+);")
    part = np.zeros((len(ref["X_prop"]), 6))
    part[:, 0] = (1.0 - t_prop) * ref["X_prop"]
    part[:, 3] = ref["K_prop"] / scale
    return part


def _matlab_fins_part(ref):
    """remus100.m 249 (fin terms), 248, 249, 251, 252, from the MATLAB columns."""
    p = _remus_fins()
    part = np.zeros((len(ref["X_r"]), 6))
    part[:, 0] = ref["X_r"] + ref["X_s"]
    part[:, 1] = ref["Y_r"]
    part[:, 2] = ref["Z_s"]
    part[:, 4] = -p["stern_plane_position"] * ref["Z_s"]
    part[:, 5] = p["rudder_position"] * ref["Y_r"]
    return part


def _block_parts(propeller_params, fin_params, ref):
    _, propeller = _build("propeller", propeller_params)
    _, fins = _build("fins", fin_params)
    prop = np.array([_tau(propeller, n=[u[2]], nu_r=nu) for u, nu in zip(ref["ui"], ref["nu_r"])])
    fin = np.array([_tau(fins, delta=u[:2], nu_r=nu) for u, nu in zip(ref["ui"], ref["nu_r"])])
    return prop, fin


def test_remus100_actuator_reference_is_consistent():
    """The CSV holds 1036 cases, reaches beyond both saturation limits, and its
    tau is the sum of its own propeller and fin terms (lines 249-254)."""
    ref = _actuators_csv()
    assert ref["tau"].shape == (1036, 6)
    assert np.any(np.abs(ref["ui"][:, 2]) > _value(REMUS, 110))
    assert np.any(np.abs(ref["ui"][:, :2]) > _value(REMUS, 109))
    total = _matlab_propeller_part(ref) + _matlab_fins_part(ref)
    assert _max_diff(total, ref["tau"]) <= G1_TOLERANCE


def test_transcriptions_reproduce_remus100_actuator_matlab():
    """The test-side transcriptions (_mss_remus_propeller, _mss_fins) equal MATLAB."""
    ref = _actuators_csv()
    pp, fp = _remus_propeller(), _remus_fins()
    for k, (u, nu) in enumerate(zip(ref["ui"], ref["nu_r"])):
        expected = _mss_remus_propeller(pp, u[2], nu) + _mss_fins(fp, u[:2], nu)
        assert _max_diff(expected, ref["tau"][k]) <= G1_TOLERANCE, k


def test_G1_MSS_propeller_block_equals_matlab():
    """Block "linearized" on remus100.m's parameters (rho 1026, printed K_T/K_Q
    as the exact override) vs MATLAB's (1 - t_prop) X_prop and K_prop / 10;
    the saturated speed equals MATLAB's n_p."""
    ref = _actuators_csv()
    p = _remus_propeller()
    constants, _ = _build("propeller", p)
    np.testing.assert_allclose(np.clip(ref["ui"][:, 2], -constants.max_shaft_speed, constants.max_shaft_speed) / 60.0,
                               ref["n_p"], atol=G1_TOLERANCE, rtol=0)
    prop, _ = _block_parts(p, _remus_fins(), ref)
    expected = _matlab_propeller_part(ref)
    for k in range(len(prop)):
        assert _max_diff(prop[k], expected[k]) <= G1_TOLERANCE, k


def test_G1_MSS_fins_block_equals_matlab():
    """Block fins on remus100.m's parameters (20 deg, rho 1026, x = -a),
    default convention, vs MATLAB's X_r, X_s, Y_r, Z_s; the saturated angles
    equal MATLAB's delta_r, delta_s."""
    ref = _actuators_csv()
    p = _remus_fins()
    np.testing.assert_allclose(np.clip(ref["ui"][:, :2], -p["max_deflection"], p["max_deflection"]),
                               ref["delta"], atol=G1_TOLERANCE, rtol=0)
    _, fin = _block_parts(_remus_propeller(), p, ref)
    expected = _matlab_fins_part(ref)
    for k in range(len(fin)):
        assert _max_diff(fin[k], expected[k]) <= G1_TOLERANCE, k


def test_G1_MSS_propeller_plus_fins_equals_matlab_tau():
    """No transcription at all: block propeller + block fins = remus100.m's tau."""
    ref = _actuators_csv()
    prop, fin = _block_parts(_remus_propeller(), _remus_fins(), ref)
    for k in range(len(prop)):
        assert _max_diff(prop[k] + fin[k], ref["tau"][k]) <= G1_TOLERANCE, k


def _worst_against_matlab(propeller_params, fin_params):
    ref = _actuators_csv()
    prop, fin = _block_parts(propeller_params, fin_params, ref)
    return (_max_diff(prop, _matlab_propeller_part(ref)), _max_diff(fin, _matlab_fins_part(ref)))


def test_G4_MSS_matlab_detects_propeller_and_fin_perturbations():
    pp, fp = _remus_propeller(), _remus_fins()
    assert max(_worst_against_matlab(pp, fp)) <= G1_TOLERANCE
    assert _worst_against_matlab({**pp, "propeller_diameter": 1.01 * pp["propeller_diameter"]}, fp)[0] > G4_FACTOR * G1_TOLERANCE
    fp_bad = {**fp, "stern_plane_lift_coefficient": 1.01 * fp["stern_plane_lift_coefficient"]}
    assert _worst_against_matlab(pp, fp_bad)[1] > G4_FACTOR * G1_TOLERANCE


def test_polynomial_coefficients_depart_from_matlab_by_the_printed_rounding_only():
    """Without the override the block computes K_T/K_Q at J = 0 and
    J = Ja_max from the polynomial; remus100.m 157-161 prints them to 4
    decimals. The block then departs from MATLAB by more than G1 but never by
    more than a 5e-5 error in each coefficient can cause."""
    ref = _actuators_csv()
    p = _remus_propeller()
    prop, _ = _block_parts({**p, "thrust_torque_coefficients": None}, _remus_fins(), ref)
    diff = np.abs(prop - _matlab_propeller_part(ref))
    half_unit = 0.5e-4   # half a unit in the 4th decimal of the printed values
    rho, d, j_max, t = p["water_density"], p["propeller_diameter"], p["max_advance_number"], p["thrust_deduction"]
    n_p, va = ref["n_p"], (1.0 - p["wake_fraction"]) * np.linalg.norm(ref["nu_r"][:, :3], axis=1)
    forward = n_p > 0
    shape = half_unit * n_p ** 2 + forward * (2.0 * half_unit / j_max) * (va / d) * np.abs(n_p)
    assert np.all(diff[:, 0] <= (1.0 - t) * rho * d ** 4 * shape + G1_TOLERANCE)
    assert np.all(diff[:, 3] <= p["roll_moment_scale"] * rho * d ** 5 * shape + G1_TOLERANCE)
    assert diff.max() > G4_FACTOR * G1_TOLERANCE


def _wageningen_csv():
    h, v = _load_csv("wageningen_kt_kq_mss_current.csv")
    return {name: v[:, h.index(name)] for name in h}


def _wageningen_csv_sets():
    """remus_like / four_blade as WAGENINGEN_SETS; outboard = the
    outboard parameter set's geometry (parameter_sets.json)."""
    p = _outboard_rpm_parameters(False)
    return {**WAGENINGEN_SETS,
            "outboard": (p["pitch_diameter_ratio"], p["blade_area_ratio"], int(p["blade_count"]))}


def _wageningen_rows(name):
    w = _wageningen_csv()
    pd, aeao, z = _wageningen_csv_sets()[name]
    rows = (w["PD"] == pd) & (w["AEAO"] == aeao) & (w["z"] == z)
    assert rows.sum() == 1252, (name, rows.sum())
    return w["J"][rows], w["KT"][rows], w["KQ"][rows]


@pytest.mark.parametrize("name", ["four_blade", "outboard", "remus_like"])
def test_G1_MSS_wageningen_block_default_equals_matlab(name):
    """MATLAB running wageningen.m, J = -0.5 ... 2.0 (the default: any J)."""
    pd, aeao, z = _wageningen_csv_sets()[name]
    _, function = _build_wageningen(pd, aeao, z)
    for j, kt, kq in zip(*_wageningen_rows(name)):
        assert _max_diff(_kt_kq(function, j), (kt, kq)) <= G1_TOLERANCE, j


@pytest.mark.parametrize("name", ["four_blade", "outboard", "remus_like"])
def test_mss_table_transcription_reproduces_matlab_wageningen(name):
    """_mss_wageningen (WageningData.mat, used by the Wageningen tests) equals
    MATLAB at every J."""
    pd, aeao, z = _wageningen_csv_sets()[name]
    for j, kt, kq in zip(*_wageningen_rows(name)):
        assert _max_diff(_mss_wageningen(j, pd, aeao, z), (kt, kq)) <= G1_TOLERANCE, j


def test_G4_MSS_wageningen_blade_area_plus_1_percent_is_detected_by_matlab():
    pd, aeao, z = _wageningen_csv_sets()["remus_like"]
    rows = list(zip(*_wageningen_rows("remus_like")))
    worst = lambda f: max(_max_diff(_kt_kq(f, j), (kt, kq)) for j, kt, kq in rows)
    assert worst(_build_wageningen(pd, aeao, z)[1]) <= G1_TOLERANCE
    assert worst(_build_wageningen(pd, 1.01 * aeao, z)[1]) > G4_FACTOR * G1_TOLERANCE


@pytest.mark.parametrize("rpm", [300.0, -300.0], ids=["ahead", "astern"])
def test_G1_MSS_outboard_rpm_default_with_matlab_kt_kq(rpm):
    """Default path (MSS polynomial, signed J) of the outboard rpm block on the
    outboard geometry's MATLAB rows: the surge speed is set so that
    J = Va / (n D) (lines 315, 382) hits each MATLAB J, steering 0.1 rad;
    expected = source lines 315-325, 380-414 with MATLAB's K_T, K_Q, times the
    reverse factors for n < 0 (lines 317-321). Covers the default path's
    reverse branch."""
    p = _outboard_rpm_parameters(False)
    _, function = _build("outboard_motor_rpm", p)
    n_rps, delta = rpm / 60.0, 0.1
    gain = p["advance_speed_factor"] * (1.0 - p["propwash_factor"])
    reverse = (p["reverse_thrust_factor"], p["reverse_torque_factor"]) if rpm < 0 else (1.0, 1.0)
    for j, kt, kq in zip(*_wageningen_rows("outboard")):
        nu_r = np.array([j * n_rps * p["propeller_diameter"] / gain, 0, 0, 0, 0, 0])
        expected = _outboard_rpm_force(p, rpm, delta, nu_r, reverse[0] * kt, reverse[1] * kq)
        got = _tau(function, n=[rpm], delta=[delta], nu_r=nu_r)
        assert _max_diff(got, expected) <= G1_TOLERANCE, (j, got, expected)


# --- Transforms only from more_transformations (owner, 2026-10-05) ---------
TRANSFORM_FUNCTION_NAMES = {"rotation_zyx", "rzyx", "rzyx_explicit", "rx", "ry", "rz", "r_fossen", "r_bn", "r_nb",
                     "skew", "smtrx", "h_matrix", "hmtrx", "t_euler", "j_body_to_eta", "gravity"}
TRANSFORMS_NUMPY = "more_transformations"
TRANSFORMS_CASADI = "more_transformations.more_casadi_transformations"


def _force_producer_modules():
    package = Path(importlib.import_module(PACKAGE).__file__).resolve().parent
    files = [path for folder in (package, package / "thrusters", package / "propulsor")
             for path in folder.glob("*.py")]
    return {path.stem: ast.parse(path.read_text()) for path in sorted(files)}


def test_transforms_no_local_rotation_or_skew_in_force_producers():
    """No module of models/force_producers defines a rotation, skew,
    H, T, J or gravity of its own."""
    local = {(name, node.name) for name, tree in _force_producer_modules().items()
             for node in ast.walk(tree)
             if isinstance(node, ast.FunctionDef) and node.name.lower() in TRANSFORM_FUNCTION_NAMES}
    assert not local, sorted(local)


def test_transforms_no_local_cross_product_in_force_producers():
    """Every skew from the library: a moment r x F is S(r) F with the
    library skew (MatrixTransforms.skew, numpy in preprocess, CasADi in the graph),
    not np.cross / ca.cross."""
    calls = {(name, node.lineno) for name, tree in _force_producer_modules().items()
             for node in ast.walk(tree)
             if isinstance(node, ast.Call) and getattr(node.func, "attr", getattr(node.func, "id", "")) == "cross"}
    assert not calls, sorted(calls)


def test_transforms_force_producers_import_only_the_casadi_transforms():
    """All-CasADi blocks (no numpy in a model module): the shaft rotation and
    every skew come from ``more_transformations.more_casadi_transformations``;
    the numpy ``more_transformations`` modules are not imported."""
    modules = set()
    for tree in _force_producer_modules().values():
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                modules.add(node.module)
            elif isinstance(node, ast.Import):
                modules.update(alias.name for alias in node.names)
    assert not any(m == TRANSFORMS_NUMPY or (m.startswith(TRANSFORMS_NUMPY + ".") and not m.startswith(TRANSFORMS_CASADI)) for m in modules), modules
    assert any(m == TRANSFORMS_CASADI or m.startswith(TRANSFORMS_CASADI + ".") for m in modules), modules


def test_transforms_shaft_axis_and_moment_arms_equal_more_transformations():
    """Numerical twin of the three transform tests: the shaft axis is
    Rzyx(orientation)[:, 0] (equal to the numpy library's), and every moment
    arm equals the CasADi library skew, S(r) F."""
    from more_transformations.matrix_transforms import MatrixTransforms as numpy_transforms
    from more_transformations.more_casadi_transformations.matrix_transforms import MatrixTransforms as casadi_transforms
    rng = _rng()
    base = _scaled_propeller()
    for _ in range(200):
        orientation = rng.uniform(-np.pi, np.pi, 3)
        constants, _ = _build("propeller", {**base, "orientation": orientation})
        np.testing.assert_allclose(constants.shaft_axis, numpy_transforms.Rzyx(orientation)[:, 0], atol=G2_TOLERANCE, rtol=0)

    p = _outboard_rpm_parameters(True)
    skew = np.array(casadi_transforms.skew(p["position"]), dtype=float)
    _, function = _build("outboard_motor_rpm", p)
    for rpm, delta, nu_r in list(zip(*_outboard_cases(p)))[:200]:
        tau = _tau(function, n=[rpm], delta=[delta], nu_r=nu_r)
        np.testing.assert_allclose(tau[3:], skew @ tau[:3], atol=G2_TOLERANCE, rtol=0)

    d = _scaled_differential()
    constants, _ = _build("differential_thruster", d)
    for i in range(2):
        direction = np.asarray(d["thruster_directions"][i]) / np.linalg.norm(d["thruster_directions"][i])
        arm = np.array(casadi_transforms.skew(d["thruster_positions"][i]), dtype=float) @ direction
        np.testing.assert_allclose(constants.allocation_matrix[3:, i], arm, atol=G2_TOLERANCE, rtol=0)


# --- Source functions not ported, justified by a test ----------------------


@pytest.mark.parametrize("convention", CONVENTIONS)
def test_physics_fins_only_drag_and_moments_are_r_cross_f(convention):
    """Physical sign, no model needed: a deflected fin never pushes the
    vehicle forward (tau_X <= 0), produces no roll moment, and its moments
    are r x F at the fin positions (rudder at x_r: N = x_r Y; stern planes at
    x_s: M = -x_s Z), computed with the numpy library skew. Holds for all four
    conventions, so the convention only chooses the sign of the deflection."""
    from more_transformations.matrix_transforms import MatrixTransforms as numpy_transforms
    p = _remus_fins()
    _, function = _build("fins", {**p, "convention": convention})
    r_rudder = np.array([p["rudder_position"], 0.0, 0.0])
    r_stern = np.array([p["stern_plane_position"], 0.0, 0.0])
    for k, (delta, nu_r) in enumerate(zip(*_fin_cases())):
        tau = _tau(function, delta=delta, nu_r=nu_r)
        assert tau[0] <= 0.0, (k, tau)
        moment = numpy_transforms.skew(r_rudder) @ np.array([0.0, tau[1], 0.0]) \
            + numpy_transforms.skew(r_stern) @ np.array([0.0, 0.0, tau[2]])
        assert _max_diff(tau[3:], moment) <= G2_TOLERANCE, (k, tau)


# ==========================================================================
# 8. The all-CasADi blocks: committed blocks, frozen path, gradients, no numpy
# ==========================================================================
COMMITTED_REVISION = "5d67caa"  # last revision of the blocks with numpy pre-processing
COMMITTED_FILES = ("__init__.py", "_common.py", "differential_thruster.py", "fins.py", "outboard_motor.py",
                   "propeller.py", "vsim_fins.py", "wageningen_kt_kq.py")
COMMITTED_SEED = 20261007
N_COMMITTED_STATES = 300
GRADIENT_TOLERANCE = 1e-8


def _committed_package(tmp_path):
    """The force-producer package as committed at ``COMMITTED_REVISION``,
    read with ``git show`` into a temporary package (the old modules are
    imported, no values are frozen in a file). Skips when git or the revision
    is not available (a source archive without history)."""
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[2]
    package = tmp_path / "force_producers_committed"
    package.mkdir()
    for name in COMMITTED_FILES:
        try:
            shown = subprocess.run(
                ["git", "-C", str(root), "show",
                 f"{COMMITTED_REVISION}:more_dynamics/models/force_producers/{name}"],
                capture_output=True, text=True, check=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            pytest.skip(f"revision {COMMITTED_REVISION} not readable with git here: {exc}")
        (package / name).write_text(shown.stdout)
    sys.path.insert(0, str(tmp_path))
    try:
        return importlib.import_module("force_producers_committed")
    finally:
        sys.path.remove(str(tmp_path))


def _committed_cases():
    """(label, block name, parameter set, command generator) for every
    producer and selector choice; commands reach beyond the saturation
    limits."""
    cases = []
    for label, p in (("otter", _otter_differential(9.81)), ("scaled", _scaled_differential())):
        cases.append((f"differential {label}", "differential_thruster", p,
                      lambda r: {"n": r.uniform(-150.0, 150.0, 2)}))
    rpm = lambda r: {"n": [r.uniform(-2000.0, 2000.0)]}
    for model in OPEN_WATER_MODES:
        for clip in (False, True):
            for label, p in (("remus", _remus_propeller()), ("scaled", _scaled_propeller()),
                             ("remus_polynomial", {**_remus_propeller(), "thrust_torque_coefficients": None})):
                cases.append((f"propeller {label} {model} clip={clip}", "propeller",
                              {**p, "open_water_model": model, "clip_advance_ratio": clip}, rpm))
    for convention in CONVENTIONS:
        cases.append((f"fins {convention}", "fins", {**_remus_fins(), "convention": convention},
                      lambda r: {"delta": r.uniform(-0.6, 0.6, 2)}))
    for scaled in (False, True):
        p = _vsim_parameters(scaled)
        cases.append((f"vsim scaled={scaled}", "vsim_fins", p,
                      lambda r, k=len(p["max_forces"]): {"delta": r.uniform(-0.6, 0.6, k)}))
        for power_limit in (True, False):
            for clip in (False, True):
                q = {**_outboard_rpm_parameters(scaled), "power_limit": power_limit, "clip_advance_ratio": clip}
                cases.append((f"outboard rpm scaled={scaled} power={power_limit} clip={clip}",
                              "outboard_motor_rpm", q,
                              lambda r, m=q["max_speed"]: {"n": [r.uniform(-1.2 * m, 1.2 * m)],
                                                           "delta": [r.uniform(-0.6, 0.6)]}))
        q = _outboard_throttle_parameters(scaled)
        cases.append((f"outboard throttle scaled={scaled}", "outboard_motor_throttle", q,
                      lambda r, m=q["max_thrust"]: {"thrust": [r.uniform(-1.2 * m, 1.2 * m)],
                                                    "delta": [r.uniform(-0.6, 0.6)]}))
    return cases


def test_G2_blocks_equal_the_committed_blocks(tmp_path):
    """Every producer and selector choice (differential 2 sets; propeller 3
    sets x 3 models x clip; fins 4 conventions; VSIM 2 sets; outboard rpm
    2 sets x power cap x clip; throttle 2 sets) and the Wageningen polynomial
    (2 geometries x clip): equal to the blocks committed at
    ``COMMITTED_REVISION`` (numpy pre-processing) on 300 seeded commands and
    states, 1e-10."""
    old = _committed_package(tmp_path)
    cases = _committed_cases()
    assert len(cases) == 2 + 18 + 4 + 2 + 8 + 2
    for label, name, p, commands in cases:
        _, function = _build(name, p)
        old_p = p
        if name == "propeller":
            # propeller_diameter / max_shaft_speed are renamed only from this round
            # (A-50, shared names): the committed revision still calls them
            # diameter / max_speed.
            old_p = {**{k: v for k, v in p.items() if k not in ("propeller_diameter", "max_shaft_speed")},
                    "diameter": p["propeller_diameter"], "max_speed": p["max_shaft_speed"]}
        old_function = getattr(old, f"{name}_casadi")(getattr(old, f"preprocess_{name}")(**old_p))
        rng = np.random.default_rng(COMMITTED_SEED)
        for k in range(N_COMMITTED_STATES):
            inputs = {**commands(rng), "nu_r": rng.uniform(-3.0, 3.0, 6)}
            assert _max_diff(_tau(function, **inputs), _tau(old_function, **inputs)) <= G2_TOLERANCE, (label, k)
    for name, (pd, aeao, z) in WAGENINGEN_SETS.items():
        for clip in (False, True):
            _, function = _build_wageningen(pd, aeao, z, clip_advance_ratio=clip)
            old_function = old.wageningen_casadi(old.preprocess_wageningen(pd, aeao, z, clip))
            for j in _wageningen_grid()[:N_COMMITTED_STATES]:
                assert _max_diff(_kt_kq(function, j), _kt_kq(old_function, j)) <= G2_TOLERANCE, (name, clip, j)


@pytest.mark.parametrize("name", ["differential_thruster", "propeller", "fins", "vsim_fins",
                                  "outboard_motor_rpm", "outboard_motor_throttle"])
def test_frozen_block_equals_the_block_called_with_numbers(name):
    """``freeze`` (what a plugin carries: a function of the commands and
    ``nu_r``) against the block called with the numbers, 300 seeded inputs;
    entry by entry within 1e-14 relative to max(1, |value|)."""
    from more_transformations.more_casadi_transformations import check_values

    label, _, p, commands = next(c for c in _committed_cases() if c[1] == name)
    _, frozen = _build(name, p)
    block = _contract(MODULE_OF.get(name, name))
    numbers = {k: v for k, v in p.items() if k not in SELECTORS[name] and k != "thrust_torque_coefficients"}
    selectors = {k: p[k] for k in SELECTORS[name] if k in p}
    if name == "vsim_fins":
        selectors["fin_count"] = len(p["max_forces"])
    if name == "propeller":
        given = p.get("thrust_torque_coefficients")
        selectors["open_water_coefficients"] = "polynomial" if given is None else "given"
        if given is not None:
            numbers["thrust_torque_coefficients"] = given
    unfrozen = getattr(block, f"{name}_casadi")(**selectors)
    declared = [n for n in unfrozen.name_in() if n not in frozen.name_in()]
    values = {n: np.asarray(numbers[n], float) for n in declared}
    rng = np.random.default_rng(COMMITTED_SEED)
    for k in range(N_COMMITTED_STATES):
        inputs = {**{n: np.asarray(v, float) for n, v in commands(rng).items()}, "nu_r": rng.uniform(-3.0, 3.0, 6)}
        a, b = frozen(**inputs), unfrozen(**inputs, **values)
        for output in frozen.name_out():
            want = np.array(b[output], dtype=float)
            error = np.abs(np.array(a[output], dtype=float) - want) / np.maximum(1.0, np.abs(want))
            assert error.max() <= 1e-14, (label, output, k)


# block -> (parameter, central-difference step, commands); tau is linear in
# each (rudder area; density; thrust coefficient with the speed inside its limits), so
# the difference is exact up to rounding
GRADIENT_CASES = {
    "fins": (lambda: _remus_fins(), "rudder_area", 1e-4, lambda r: {"delta": r.uniform(-0.3, 0.3, 2)}),
    "propeller": (lambda: _remus_propeller(), "water_density", 1.0, lambda r: {"n": [r.uniform(-1500.0, 1500.0)]}),
    "differential_thruster": (lambda: _otter_differential(9.81), "positive_thrust_coefficients", 1e-4,
                              lambda r: {"n": r.uniform(10.0, 60.0, 2)}),
}


@pytest.mark.parametrize("name", sorted(GRADIENT_CASES))
def test_gradient_of_tau_with_respect_to_a_parameter(name):
    """Identification path: the block called with one declared parameter left
    as a symbol; d tau / d parameter from CasADi equals a central difference
    of the block called with numbers, 20 seeded inputs, 1e-8."""
    import casadi as ca

    make, parameter, step, commands = GRADIENT_CASES[name]
    p = make()
    block = _contract(name)
    selectors = {k: p[k] for k in SELECTORS[name] if k in p}
    numbers = {k: ca.DM(np.asarray(v, float)) for k, v in p.items() if k not in SELECTORS[name]}
    if name == "propeller":
        selectors["open_water_coefficients"] = "given"
    function = getattr(block, f"{name}_casadi")(**selectors)
    base = np.array(numbers[parameter], dtype=float).reshape(-1)
    symbol = ca.SX.sym(parameter, base.size)
    rng = np.random.default_rng(COMMITTED_SEED)

    def tau_of(x, inputs):
        return function(**inputs, **{**numbers, parameter: x})["tau"]

    found = 0.0
    for k in range(20):
        inputs = {**{n: np.asarray(v, float) for n, v in commands(rng).items()}, "nu_r": rng.uniform(-2.0, 2.0, 6)}
        gradient = ca.Function("gradient", [symbol], [ca.jacobian(tau_of(symbol, inputs), symbol)])
        exact = np.array(gradient(base), dtype=float)
        for i in range(base.size):
            shift = np.zeros(base.size)
            shift[i] = step
            upper = np.array(tau_of(base + shift, inputs), dtype=float).ravel()
            lower = np.array(tau_of(base - shift, inputs), dtype=float).ravel()
            assert _max_diff(exact[:, i], (upper - lower) / (2.0 * step)) <= GRADIENT_TOLERANCE, (name, k, i)
        found = max(found, np.abs(exact).max())
    assert found > 1e-2, name  # the gradient is not trivially zero


def test_blocks_import_no_numpy():
    """No module of models/force_producers imports numpy or scipy, or the
    numpy ``more_transformations`` modules (AST scan of every file)."""
    for name, tree in _force_producer_modules().items():
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                modules = [node.module]
            for module in modules:
                assert module.split(".")[0] not in ("numpy", "scipy"), (name, module)
                if module.startswith("more_transformations"):
                    assert module.startswith(TRANSFORMS_CASADI), (name, module)
