"""The contract of the fin skeleton and its parts, as the gate tests see it.

Written 2026-10-08, before the parts exist. Every test module of this folder
imports its module paths, signatures, tolerances, reference loaders and
builders from here, so a module that moves is one line in ``PART_MODULES``.

A fin is one skeleton, ``LiftingFin``, with five typed sockets (servo,
inflow, flow angle, interference, section). Each socket holds one part; a
part is one CasADi function. Alternative physics are separate part modules,
named by what they compute, never by where they were read; the source is in
the citation. Several force producers (fins, one-fin leaves) are summed by
one generic ``ForceProducerSet`` behind a command map.

Model layer (``more_dynamics.models.force_producers``)
------------------------------------------------------
* A part form ``<form>`` lives in its own module (``PART_MODULES``) and
  exposes ``<form>_parameters() -> tuple[Parameter, ...]`` and
  ``<form>_casadi() -> ca.Function``. The function's inputs are the
  coupling inputs of its type (``COUPLINGS``, in that order), then one
  input per declared parameter, in declaration order; its outputs are the
  type's outputs, in that order, and nothing else. A part never declares
  ``water_density``.
* ``fin_parts.lifting_fin``: ``lifting_fin_parameters()`` (``fin_position``
  3x1 m, ``chord_axis`` 3x1, ``lift_axis`` 3x1, ``fin_area`` m^2 > 0; no
  ``water_density``), ``check_lifting_fin_values(values)`` (the checked
  numbers; refuses a chord or lift axis that is not a unit vector and axes
  that are not perpendicular) and ``lifting_fin_casadi(parts)`` with
  ``parts`` a dict slot -> part function (frozen or with parameters open).
  It refuses (``ValueError`` naming the slot) a missing or unknown slot, a
  part without a coupling name of its type, a coupling of the wrong size,
  and a part with an input named ``water_density``. Its function: inputs
  ``command`` (1), ``state`` (the servo's), ``nu_r`` (6, BODY FRD, relative
  to the water), ``water_density`` (1), the skeleton's parameters, then
  every still-open part parameter as ``"<slot>.<name>"``; outputs ``tau``
  (6, N and N m about the CO, BODY), ``state_dot``, ``deflection``,
  ``angle_of_attack``.
* ``force_producer_set``: ``force_producer_set_parameters(producers,
  command_count)`` (``command_map``, sum of the producers' command sizes x
  ``command_count``) and ``force_producer_set_casadi(producers,
  command_count)``; every producer has inputs starting ``command, state,
  nu_r, water_density`` and outputs starting ``tau, state_dot`` (else
  ``ValueError`` naming the producer's index). The set's inputs:
  ``command``, ``state`` (stacked in list order), ``nu_r``,
  ``water_density`` (one, fed to every producer), ``command_map``, then the
  producers' open inputs as ``"producers[i].<name>"``; outputs ``tau`` (the
  sum) and ``state_dot`` (stacked).
* ``fin_parts.proportional_force_fin``: one fin of DUNE's simulator law (force
  proportional to the deflection and to the speed squared along a fixed
  direction, DUNE's moment arms), the per-fin leaf of the existing
  ``vsim_fins`` block. ``proportional_force_fin_parameters()``
  (``max_force`` 3x1 N, ``fin_position`` 3x1 m, ``max_deflection`` rad > 0)
  and ``proportional_force_fin_casadi()`` with the producer interface; it
  takes ``water_density`` and does not use it (DUNE's force gain holds its
  own density).

The vehicle feeds ``nu_r`` and ``water_density`` (owner's answers of
2026-10-08: the velocity input by convention after the commands; the
density fed as an input, declared by no part).

Author:    Enio Krizman
Date:      2026-10-08
"""

import hashlib
import importlib
import json
import os
from pathlib import Path

import numpy as np
import pytest

PACKAGE = "more_dynamics.models.force_producers"
DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "force_producers" / "fin_parts"
PARENT_DATA_DIR = DATA_DIR.parent

G1_TOLERANCE = 1e-9    # against MATLAB running MSS, absolute (as tests/force_producers)
G2_TOLERANCE = 1e-10   # against a transcription of a cited source or another block, absolute
G4_FACTOR = 10.0       # a perturbed model must differ by more than this x G1
SEED = 20261008
N_RANDOM = 500

# Servo settings: one consolidated module (``servo.servo``) with three
# switches (``servo.servo.DYNAMICS``, ``.ANGLE_LIMIT``); a "form" name here
# is a fixed combination of them, not a separate module. There is no
# rate-only servo: no module, no "DUNE" or "1:1" claim anywhere.
SERVO_SETTINGS = {
    "ideal": {"dynamics": "none", "rate_limit": False, "angle_limit": "on_command"},
    "no_limits": {"dynamics": "none", "rate_limit": False, "angle_limit": "none"},
    "lag_rate_angle": {"dynamics": "first_order_lag", "rate_limit": True, "angle_limit": "on_command"},
    "output_saturated": {"dynamics": "first_order_lag", "rate_limit": True, "angle_limit": "on_output"},
}

# (slot, form) -> module under PACKAGE. Forms of the waiting list have no module yet.
PART_MODULES = {
    ("servo", "ideal"): "fin.servo.servo",
    ("servo", "no_limits"): "fin.servo.servo",
    ("servo", "lag_rate_angle"): "fin.servo.servo",
    ("servo", "output_saturated"): "fin.servo.servo",
    ("inflow", "translational"): "fin.inflow.translational",
    ("inflow", "rigid_point"): "fin.inflow.rigid_point",
    ("flow_angle", "none"): "fin.flow_angle.none",
    ("flow_angle", "small_angle"): "fin.flow_angle.small_angle",
    ("interference", "none"): "fin.interference.none",
    ("section", "quadratic_drag"): "fin.section.quadratic_drag",
    ("section", "linear_section"): "fin.section.linear_section",
}
# Second-fidelity forms whose sources are not yet read in full (Pitts 1957,
# Bhattacharyya / lifting-line and stall references): their tests are written
# and skipped with this reason.
WAITING = {
    ("section", "lifting_line"): "fin.section.lifting_line",
    ("interference", "slender_body"): "fin.interference.slender_body",
}
WAITING_REASON = ("waiting on the second-fidelity source reading (Pitts 1957 slender-body "
                  "factors, the lifting-line induced drag): form and parameter names not fixed")

SKELETON_MODULE = "fin.lifting_fin"
SET_MODULE = "force_producer_set"
LEAF_MODULE = "fin_parts.proportional_force_fin"

SLOTS = ("servo", "inflow", "flow_angle", "interference", "section")
# type -> (coupling inputs {name: rows, None = any}, outputs {name: rows, None = any})
COUPLINGS = {
    "servo": ({"command": 1, "servo_state": None}, {"deflection": 1, "servo_state_dot": None}),
    "inflow": ({"nu_r": 6, "fin_position": 3}, {"fin_velocity": 3}),
    "flow_angle": ({"fin_velocity": 3, "chord_axis": 3, "lift_axis": 3}, {"flow_angle": 1, "speed_squared": 1}),
    "interference": ({}, {"deflection_factor": 1, "flow_angle_factor": 1}),
    "section": ({"angle_of_attack": 1}, {"lift_coefficient": 1, "drag_coefficient": 1}),
}
PRODUCER_INPUTS = ["command", "state", "nu_r", "water_density"]
PRODUCER_OUTPUTS = ["tau", "state_dot"]
FIN_OUTPUTS = ["tau", "state_dot", "deflection", "angle_of_attack"]

# Declarations as the tests expect them: name -> (shape, unit, minimum, minimum_exclusive)
DECLARED = {
    ("servo", "ideal"): {"max_deflection": ((1, 1), "rad", 0.0, True)},
    ("servo", "no_limits"): {},
    ("servo", "lag_rate_angle"): {"max_deflection": ((1, 1), "rad", 0.0, True),
                                  "max_rate": ((1, 1), "rad/s", 0.0, True),
                                  "time_constant": ((1, 1), "s", 0.0, True)},
    ("servo", "output_saturated"): {"max_deflection": ((1, 1), "rad", 0.0, True),
                                    "max_rate": ((1, 1), "rad/s", 0.0, True),
                                    "time_constant": ((1, 1), "s", 0.0, True)},
    ("inflow", "translational"): {},
    ("inflow", "rigid_point"): {},
    ("flow_angle", "none"): {},
    ("flow_angle", "small_angle"): {},
    ("interference", "none"): {},
    ("section", "quadratic_drag"): {"lift_slope": ((1, 1), "1/rad", 0.0, True)},
    ("section", "linear_section"): {"lift_slope": ((1, 1), "1/rad", 0.0, True),
                                    "zero_lift_drag": ((1, 1), "1", 0.0, False)},
}
SKELETON_DECLARED = {"fin_position": ((3, 1), "m", None, False),
                     "chord_axis": ((3, 1), "1", None, False),
                     "lift_axis": ((3, 1), "1", None, False),
                     "fin_area": ((1, 1), "m^2", 0.0, True)}
LEAF_DECLARED = {"max_force": ((3, 1), "N", None, False),
                 "fin_position": ((3, 1), "m", None, False),
                 "max_deflection": ((1, 1), "rad", 0.0, True)}

# Unit axes of BODY (FRD)
E_X = [1.0, 0.0, 0.0]
E_Y = [0.0, 1.0, 0.0]
E_Z = [0.0, 0.0, 1.0]
NEG_Y = [0.0, -1.0, 0.0]
NEG_Z = [0.0, 0.0, -1.0]


# --------------------------------------------------------------------------
# Modules and functions
# --------------------------------------------------------------------------
def module(relative):
    """Import ``PACKAGE.<relative>``; a missing module fails the test with
    the name of the part it waits for."""
    try:
        return importlib.import_module(f"{PACKAGE}.{relative}")
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed: {exc}")
        pytest.fail(f"not ported yet: {PACKAGE}.{relative} ({exc})")


def part_module(slot, form):
    return module(PART_MODULES[(slot, form)])


def part_declared(slot, form):
    if slot == "servo":
        return part_module(slot, form).servo_parameters(**SERVO_SETTINGS[form])
    return getattr(part_module(slot, form), f"{form}_parameters")()


def part_block(slot, form):
    """The part function with every own parameter open."""
    if slot == "servo":
        return part_module(slot, form).servo_casadi(**SERVO_SETTINGS[form])
    return getattr(part_module(slot, form), f"{form}_casadi")()


def part(slot, form, values=None):
    """The part function with its numbers checked and frozen in (``values``
    None: parameters left open)."""
    from more_transformations.more_casadi_transformations import freeze

    block = part_block(slot, form)
    if values is None:
        return block
    return freeze(block, part_declared(slot, form), values)


def skeleton():
    return module(SKELETON_MODULE)


def fin(parts, geometry):
    """A ``LiftingFin`` from part functions and its geometry, the geometry
    checked and frozen in: a function of ``command, state, nu_r,
    water_density`` and any still-open part parameter."""
    from more_transformations.more_casadi_transformations import freeze

    s = skeleton()
    values = s.check_lifting_fin_values(geometry)
    return freeze(s.lifting_fin_casadi(parts), s.lifting_fin_parameters(), values)


def fin_from_forms(forms, part_values, geometry):
    """``forms`` slot -> form name; ``part_values`` slot -> numbers (absent:
    the part has none, or its parameters stay open when the value is None)."""
    parts = {}
    for slot in SLOTS:
        form = forms[slot]
        values = part_values.get(slot)
        parts[slot] = part(slot, form, values if part_declared(slot, form) else None)
    return fin(parts, geometry)


def geometry(position, lift_axis, area, chord_axis=E_X):
    return {"fin_position": list(position), "chord_axis": list(chord_axis),
            "lift_axis": list(lift_axis), "fin_area": float(area)}


def producer_set(producers, command_map):
    """A ``ForceProducerSet`` with its command map frozen in."""
    from more_transformations.more_casadi_transformations import freeze

    s = module(SET_MODULE)
    command_map = np.asarray(command_map, dtype=float)
    k = command_map.shape[1]
    block = s.force_producer_set_casadi(producers, k)
    return freeze(block, s.force_producer_set_parameters(producers, k), {"command_map": command_map.tolist()})


def leaf(values=None):
    from more_transformations.more_casadi_transformations import freeze

    m = module(LEAF_MODULE)
    block = m.proportional_force_fin_casadi()
    if values is None:
        return block
    return freeze(block, m.proportional_force_fin_parameters(), values)


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------
def call(function, **inputs):
    """Numbers out of a CasADi function, every output as a flat array."""
    out = function(**{k: np.asarray(v, dtype=float) for k, v in inputs.items()})
    return {k: np.array(v, dtype=float).reshape(-1) for k, v in out.items()}


def tau(function, command, nu_r, water_density, state=None):
    n = function.size1_in("state")
    state = np.zeros(n) if state is None else state
    return call(function, command=command, state=state, nu_r=nu_r, water_density=water_density)["tau"]


def max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def rng():
    return np.random.default_rng(SEED)


def random_nu_r(generator, n, speed=3.0, rate=1.0):
    return np.hstack([generator.uniform(-speed, speed, (n, 3)), generator.uniform(-rate, rate, (n, 3))])


# --------------------------------------------------------------------------
# Reference data
# --------------------------------------------------------------------------
MSS_FINS_CSV = "remus100_fins_mss_cc07579.csv"
PRESTERO_JSON = "prestero_2001_remus_fins.json"


def mss_fins():
    """MATLAB running MSS ``remus100.m`` (cc07579), propeller at rest: one
    dict of columns (SOURCE.md of this folder)."""
    path = DATA_DIR / MSS_FINS_CSV
    header = path.read_text().splitlines()[0].split(",")
    values = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {name: values[:, i] for i, name in enumerate(header)}


def mss_columns(ref, prefix, count):
    return np.column_stack([ref[f"{prefix}_{i:02d}"] for i in range(1, count + 1)])


def mss_constant(ref, name):
    column = ref[name]
    assert np.all(column == column[0]), f"{name} is not constant in {MSS_FINS_CSV}"
    return float(column[0])


def prestero():
    return json.loads((DATA_DIR / PRESTERO_JSON).read_text())


def prestero_value(name):
    return float(prestero()["parameters"][name]["value"])


def parameter_set(name):
    """A parameter set of ``tests/data/force_producers/parameter_sets.json``."""
    return json.loads((PARENT_DATA_DIR / "parameter_sets.json").read_text())[name]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def external(variable, relative):
    """A file outside this repository, found only through an environment
    variable (default unset); the test skips with the variable's name."""
    root = os.environ.get(variable)
    if not root:
        pytest.skip(f"{variable} is not set (nothing relative to one machine); needed for {relative}")
    path = Path(root) / relative
    if not path.exists():
        pytest.skip(f"{variable}: cited file not found: {path}")
    return path
