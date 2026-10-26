"""The contract of the propulsor skeleton and its parts, as the gate tests see it.

Written 2026-10-10, before the parts exist. Every test module of this folder
imports its module paths, signatures, tolerances, reference loaders and
builders from here, so a module that moves is one line in ``FORMS``.

A propulsor is one skeleton, ``Propulsor``, with six typed sockets (drive,
steering, mounting, inflow, law, housing). Each socket holds one part; a part
is one CasADi function. Alternative physics are separate part modules, named
by what they compute, never by where they were read or which craft uses
them; features of one device that switch on or off (a servo's lag, rate limit
and angle limit) are build-time selectors of one part. Several propulsors
(two fixed thrusters, a rudder and a propeller) are several force producers
of the vehicle's own list, summed by the vehicle: no producer set.

Model layer (``more_dynamics.models``: ``propulsor``, ``servo``)
-----------------------------------------------------------------
* A part form lives in its own module (``FORMS``) and exposes
  ``<stem>_parameters(**kwargs) -> tuple[Parameter, ...]`` and
  ``<stem>_casadi(**kwargs) -> ca.Function``; both take the same keywords
  (the variant's ``kwargs`` in ``FORMS``: unit and coefficient selectors).
  The function's inputs are the coupling inputs of its type (``COUPLINGS``,
  in that order), then one input per declared parameter, in declaration
  order; its outputs are the type's outputs, in that order, and nothing
  else. No part declares ``water_density`` as a parameter; the ``law`` and
  ``housing`` types take it as a coupling input.
* ``propulsor.propulsor``: ``propulsor_parameters()`` (only
  ``shaft_torque_scale``: the share of the shaft torque that reacts on the
  body), ``check_propulsor_values(values)`` and ``propulsor_casadi(parts,
  units)`` with ``parts`` a dict slot -> part function (frozen or with
  parameters open) and ``units`` a dict ``{"drive": unit of the drive's
  output, "law": unit the law accepts}`` (a proposal not yet accepted: the unit
  travels with the pair, a mismatch is refused). It
  refuses (``ValueError`` naming the slot) a missing or unknown slot, a part
  without a coupling name of its type, a coupling of the wrong size, a part
  of a slot other than ``law`` and ``housing`` with an input named
  ``water_density``, a drive whose shaft output depends on the load torque
  (load-dependent speed must be a state), and a drive / law pair whose units
  differ (naming both slots and both units). Its function: inputs
  ``command`` (the drive's command then the steering's), ``state`` (the
  drive's state then the steering's), ``nu_r`` (6, BODY FRD, relative to the
  water), ``water_density`` (1), ``shaft_torque_scale``, then every
  still-open part parameter as ``"<slot>.<name>"``; outputs ``tau`` (6, N and
  N m about the CO, BODY), ``state_dot``, ``thrust`` (the thrust along the
  axis after deduction, power cap and limit, N), ``steering_angle`` (rad).
  The steps, in this order: steering -> mounting -> inflow -> drive (speed
  without the load) -> law -> drive again with the load ``torque_factor *
  torque`` (state derivative, thrust scale, thrust limit) -> ``thrust =
  sat(thrust_scale * thrust_factor * T, +-thrust_limit)`` -> ``tau =
  [thrust a; r x thrust a] + [0; shaft_torque_scale * torque_factor * Q a] +
  housing`` (Fossen 2011, eq. 12.226, p. 400).
* Part modules of a drive expose ``DRIVE_OUTPUT_UNIT`` ("rev/s" for a shaft
  speed, "1" for a normalised effort) and of a law ``LAW_INPUT_UNIT``.
* ``servo.servo`` is the steering part of the servo kind (one device with
  switches); the fixed-thruster and nozzle forms of the steering
  socket are in ``propulsor.steering``.

The vehicle feeds ``nu_r`` and ``water_density`` (nothing hardcoded: the density comes from
the vehicle today and from its site component later). The command is ``[drive command; steering
command]`` and the state ``[drive state; steering state]``, the order of the
shared repository's jet nozzle (``desired_thrust``, ``desired_normalized_
nozzle_angle``; states ``thrust``, ``nozzle_angle``), so his payload and ours
line up entry for entry.

Proposed, not yet accepted (marked wherever the tests rely on them): the
``drive_output`` unit check; that the waterjet laws give no negative thrust
for a reversed shaft (reverse only through a reversing bucket, which is not
contracted); the momentum-flux laws as the default jet law (the tests build
them as options).

Author:    Enio Krizman
Date:      2026-10-10
"""

import csv
import hashlib
import importlib
import json
import re
from pathlib import Path
from types import SimpleNamespace

import casadi as ca
import numpy as np
import pytest

PACKAGE = "more_dynamics.models"
HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
FORCE_DATA_DIR = HERE.parent / "data"                                  # tests/force_producers/data
VEHICLE_DATA_DIR = HERE.parents[1] / "vehicle_models" / "data"         # tests/vehicle_models/data

G1_TOLERANCE = 1e-9    # against MATLAB running MSS, absolute (as tests/force_producers)
G2_TOLERANCE = 1e-10   # against a transcription of a cited source or another block, absolute
EXACT_TOLERANCE = 1e-12  # forms that must reproduce another implementation exactly
G4_FACTOR = 10.0       # a perturbed model must differ by more than this x G1
SEED = 20261010
N_RANDOM = 300

SLOTS = ("drive", "steering", "mounting", "inflow", "law", "housing")
SKELETON_MODULE = "propulsor.propulsor"

# Servo settings of the steering socket (module ``servo.servo``, a name here is a fixed combination of
# the three switches, not a separate module.
SERVO_SETTINGS = {
    "ideal": {"dynamics": "none", "rate_limit": False, "angle_limit": "on_command"},
    "no_limits": {"dynamics": "none", "rate_limit": False, "angle_limit": "none"},
    "lag_rate_angle": {"dynamics": "first_order_lag", "rate_limit": True, "angle_limit": "on_command"},
    "output_saturated": {"dynamics": "first_order_lag", "rate_limit": True, "angle_limit": "on_output"},
}

# Sign convention of the azimuth mounting: the thrust axis is [cos d, sin d, 0] in BODY (FRD), so a positive steering
# angle pushes the stern to starboard (+y). The outboard block's docstring calls the angle "port-positive" and its
# formula says otherwise (which one is wrong is a decision still to be taken): this one
# constant is the line to change, with the steering test, if the formula turns out to be the wrong one.
STEERING_SIGN = +1.0


def _form(module, stem=None, **kwargs):
    return SimpleNamespace(module=module, stem=stem or module.rsplit(".", 1)[1], kwargs=kwargs)


# (slot, variant) -> where it lives and with which keywords it is built. A variant is a form plus its selectors.
FORMS = {
    ("drive", "speed_command_rpm"): _form("propulsor.drive.speed_command", command_unit="rpm"),
    ("drive", "speed_command_rad_s"): _form("propulsor.drive.speed_command", command_unit="rad/s"),
    ("drive", "power_limited_speed_command_rpm"): _form("propulsor.drive.power_limited_speed_command",
                                                       command_unit="rpm"),
    ("drive", "normalized_lag"): _form("propulsor.drive.normalized_lag"),
    ("drive", "speed_command_normalized"): _form("propulsor.drive.speed_command_normalized"),
    ("drive", "thrust_command"): _form("propulsor.drive.thrust_command"),
    ("drive", "thrust_command_lag"): _form("propulsor.drive.thrust_command_lag"),
    ("steering", "fixed"): _form("propulsor.steering.fixed"),
    ("steering", "servo_ideal"): _form("servo.servo", "servo", **SERVO_SETTINGS["ideal"]),
    ("steering", "servo_no_limits"): _form("servo.servo", "servo", **SERVO_SETTINGS["no_limits"]),
    ("steering", "servo_lag_rate_angle"): _form("servo.servo", "servo", **SERVO_SETTINGS["lag_rate_angle"]),
    ("steering", "servo_output_saturated"): _form("servo.servo", "servo", **SERVO_SETTINGS["output_saturated"]),
    ("steering", "jet_nozzle_clamp"): _form("propulsor.steering.jet_nozzle_clamp"),
    ("mounting", "fixed_euler"): _form("propulsor.mounting.fixed_euler"),
    ("mounting", "fixed_direction"): _form("propulsor.mounting.fixed_direction"),
    ("mounting", "azimuth"): _form("propulsor.mounting.azimuth"),
    ("mounting", "jet_nozzle_direction"): _form("propulsor.mounting.jet_nozzle_direction"),
    ("mounting", "steered_rotation"): _form("propulsor.mounting.steered_rotation"),
    ("inflow", "none"): _form("propulsor.inflow.none"),
    ("inflow", "speed_magnitude"): _form("propulsor.inflow.speed_magnitude"),
    ("inflow", "surge_component"): _form("propulsor.inflow.surge_component"),
    ("inflow", "horizontal_speed_floored"): _form("propulsor.inflow.horizontal_speed_floored"),
    ("law", "open_water_bollard"): _form("propulsor.law.open_water_bollard", open_water_coefficients="polynomial"),
    ("law", "open_water_linearized"): _form("propulsor.law.open_water_linearized", open_water_coefficients="polynomial"),
    ("law", "open_water_linearized_given"): _form("propulsor.law.open_water_linearized", open_water_coefficients="given"),
    ("law", "open_water_polynomial"): _form("propulsor.law.open_water_polynomial", open_water_coefficients="polynomial"),
    ("law", "open_water_signed_advance"): _form("propulsor.law.open_water_signed_advance"),
    ("law", "quadratic_asymmetric"): _form("propulsor.law.quadratic_asymmetric"),
    ("law", "normalized_linear"): _form("propulsor.law.normalized_linear"),
    ("law", "jet_bollard_quadratic"): _form("propulsor.law.jet_bollard_quadratic"),
    ("law", "jet_bilinear_speed"): _form("propulsor.law.jet_bilinear_speed"),
    ("law", "jet_momentum_flux_given"): _form("propulsor.law.jet_momentum_flux", flow_coefficient_method="given"),
    ("law", "jet_momentum_flux_bollard_point"): _form("propulsor.law.jet_momentum_flux",
                                                      flow_coefficient_method="from_bollard_point"),
    ("law", "jet_pump_system_match"): _form("propulsor.law.jet_pump_system_match"),
    ("law", "thrust_power_limited"): _form("propulsor.law.thrust_power_limited"),
    ("housing", "none"): _form("propulsor.housing.none"),
    ("housing", "leg_normal_force"): _form("propulsor.housing.leg_normal_force"),
}

# type -> (coupling inputs {name: rows, None = any}, outputs {name: rows, None = any}), in this order
COUPLINGS = {
    "drive": ({"command": None, "drive_state": None, "shaft_torque": 1},
              {"drive_output": 1, "drive_state_dot": None, "thrust_scale": 1, "thrust_limit": 1}),
    "steering": ({"command": None, "servo_state": None}, {"deflection": 1, "servo_state_dot": None}),
    "mounting": ({"steering_angle": 1}, {"thrust_axis": 3, "thrust_position": 3}),
    "inflow": ({"nu_r": 6, "thrust_axis": 3, "thrust_position": 3},
               {"advance_speed": 1, "thrust_factor": 1, "torque_factor": 1}),
    "law": ({"drive_output": 1, "advance_speed": 1, "water_density": 1}, {"thrust": 1, "torque": 1}),
    "housing": ({"nu_r": 6, "thrust_axis": 3, "thrust_position": 3, "steering_angle": 1, "water_density": 1},
                {"tau": 6}),
}
DENSITY_SLOTS = ("law", "housing")        # the types whose coupling inputs include water_density
PRODUCER_INPUTS = ["command", "state", "nu_r", "water_density"]
PRODUCER_OUTPUTS = ["tau", "state_dot", "thrust", "steering_angle"]


def _decl(shape, unit, minimum=None, exclusive=False, maximum=None):
    return (shape, unit, minimum, exclusive, maximum)


S = (1, 1)
V3 = (3, 1)
POSITIVE = lambda unit: _decl(S, unit, 0.0, True)                    # noqa: E731
FRACTION = _decl(S, "1", 0.0, False, 1.0)
ANY = lambda unit: _decl(S, unit)                                    # noqa: E731

# Declarations as the tests expect them: name -> (shape, unit, minimum, minimum_exclusive, maximum), in order.
# The Wageningen geometry names, shapes and ranges are the existing block's (``wageningen_parameters``) and are
# compared against it, not retyped here.
WAGENINGEN = "wageningen_geometry"   # placeholder expanded by ``expected_declared``
_OPEN_WATER = {"diameter": POSITIVE("m"), WAGENINGEN: None, "max_advance_number": POSITIVE("1")}
DECLARED = {
    ("drive", "speed_command_rpm"): {"min_command": ANY("rpm"), "max_command": ANY("rpm")},
    ("drive", "speed_command_rad_s"): {"min_command": ANY("rad/s"), "max_command": ANY("rad/s")},
    ("drive", "power_limited_speed_command_rpm"): {
        "min_command": ANY("rpm"), "max_command": ANY("rpm"), "max_power": POSITIVE("W"),
        "efficiency": _decl(S, "1", 0.0, True, 1.0), "thrust_saturation": POSITIVE("N")},
    ("drive", "normalized_lag"): {"thrust_rise_time": POSITIVE("s"), "thrust_fall_time": POSITIVE("s")},
    ("drive", "speed_command_normalized"): {"max_shaft_speed": POSITIVE("rev/s"),
                                            "min_shaft_speed": _decl(S, "rev/s", None, False, 0.0)},
    ("drive", "thrust_command"): {},
    ("drive", "thrust_command_lag"): {"thrust_time_constant": POSITIVE("s")},
    ("steering", "fixed"): {},
    ("steering", "servo_ideal"): {"max_deflection": POSITIVE("rad")},
    ("steering", "servo_no_limits"): {},
    ("steering", "servo_lag_rate_angle"): {"max_deflection": POSITIVE("rad"), "max_rate": POSITIVE("rad/s"),
                                           "time_constant": POSITIVE("s")},
    ("steering", "servo_output_saturated"): {"max_deflection": POSITIVE("rad"), "max_rate": POSITIVE("rad/s"),
                                             "time_constant": POSITIVE("s")},
    ("steering", "jet_nozzle_clamp"): {"min_angle": ANY("rad"), "max_angle": ANY("rad"),
                                       "nozzle_velocity": POSITIVE("rad/s")},
    ("mounting", "fixed_euler"): {"position": _decl(V3, "m"), "orientation": _decl(V3, "rad")},
    ("mounting", "fixed_direction"): {"position": _decl(V3, "m"), "direction": _decl(V3, "1")},
    ("mounting", "azimuth"): {"position": _decl(V3, "m")},
    ("mounting", "jet_nozzle_direction"): {"position": _decl(V3, "m"), "trim_angle": ANY("rad"),
                                           "yaw_bias": ANY("rad")},
    ("mounting", "steered_rotation"): {"position": _decl(V3, "m"), "trim_angle": ANY("rad"),
                                       "yaw_bias": ANY("rad")},
    ("inflow", "none"): {},
    ("inflow", "speed_magnitude"): {"thrust_deduction": FRACTION, "wake_fraction": FRACTION},
    ("inflow", "surge_component"): {"thrust_deduction": FRACTION, "advance_speed_factor": ANY("1"),
                                    "wake_fraction": FRACTION},
    ("inflow", "horizontal_speed_floored"): {"speed_floor": POSITIVE("m/s")},
    ("law", "open_water_bollard"): _OPEN_WATER,
    ("law", "open_water_linearized"): _OPEN_WATER,
    ("law", "open_water_linearized_given"): {**_OPEN_WATER, "thrust_torque_coefficients": _decl((4, 1), "1")},
    ("law", "open_water_polynomial"): _OPEN_WATER,
    ("law", "open_water_signed_advance"): {"diameter": POSITIVE("m"), WAGENINGEN: None,
                                           "reverse_thrust_factor": POSITIVE("1"),
                                           "reverse_torque_factor": POSITIVE("1")},
    ("law", "quadratic_asymmetric"): {"positive_thrust_coefficient": POSITIVE("N*s^2/rad^2"),
                                      "negative_thrust_coefficient": POSITIVE("N*s^2/rad^2")},
    ("law", "normalized_linear"): {"max_thrust": POSITIVE("N"), "thrust_coefficient": ANY("1")},
    ("law", "jet_bollard_quadratic"): {"bollard_thrust_coefficient": POSITIVE("N*s^2/rev^2")},
    ("law", "jet_bilinear_speed"): {"speed_squared_coefficient": POSITIVE("N*s^2/rev^2"),
                                    "speed_advance_coefficient": ANY("N*s^2/(rev*m)")},
    ("law", "jet_momentum_flux_given"): {"nozzle_area": POSITIVE("m^2"), "flow_coefficient": POSITIVE("m^3/rev")},
    ("law", "jet_momentum_flux_bollard_point"): {"nozzle_area": POSITIVE("m^2"), "bollard_thrust": POSITIVE("N"),
                                                 "bollard_shaft_speed": POSITIVE("rev/s")},
    ("law", "jet_pump_system_match"): {
        "impeller_diameter": POSITIVE("m"), "nozzle_area": POSITIVE("m^2"),
        "pump_head_coefficient_at_zero_flow": POSITIVE("1"), "pump_head_slope": _decl(S, "1", 0.0),
        "nozzle_loss_coefficient": _decl(S, "1", 0.0), "inlet_loss_coefficient": FRACTION,
        "pump_efficiency": _decl(S, "1", 0.0, True, 1.0)},
    ("law", "thrust_power_limited"): {"max_thrust": POSITIVE("N"), "delivered_power": POSITIVE("W")},
    ("housing", "none"): {},
    ("housing", "leg_normal_force"): {"leg_area": POSITIVE("m^2"), "leg_aspect_ratio": POSITIVE("1")},
}
SKELETON_DECLARED = {"shaft_torque_scale": _decl(S, "1", 0.0)}

# Units of the pair drive / law (a proposal not yet accepted)
DRIVE_UNITS = {"speed_command_rpm": "rev/s", "speed_command_rad_s": "rev/s",
               "power_limited_speed_command_rpm": "rev/s", "normalized_lag": "1",
               "speed_command_normalized": "rev/s", "thrust_command": "N", "thrust_command_lag": "N"}
LAW_UNITS = {"normalized_linear": "1", "thrust_power_limited": "N"}   # every other law accepts "rev/s"

# The waterjet laws (options; none is yet the default jet law)
JET_LAWS = ("jet_bollard_quadratic", "jet_bilinear_speed", "jet_momentum_flux_given",
            "jet_momentum_flux_bollard_point", "jet_pump_system_match")

# Unit axes of BODY (FRD)
E_X = [1.0, 0.0, 0.0]
E_Y = [0.0, 1.0, 0.0]
E_Z = [0.0, 0.0, 1.0]


# --------------------------------------------------------------------------
# Modules and functions
# --------------------------------------------------------------------------
def module(relative):
    """Import ``PACKAGE.<relative>``; a missing module fails the test with the name of the part it waits for."""
    try:
        return importlib.import_module(f"{PACKAGE}.{relative}")
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed: {exc}")
        pytest.fail(f"not ported yet: {PACKAGE}.{relative} ({exc})")


def form_module(slot, variant):
    return module(FORMS[(slot, variant)].module)


def _stem_function(slot, variant, suffix):
    form = FORMS[(slot, variant)]
    mod = form_module(slot, variant)
    name = f"{form.stem}_{suffix}"
    if not hasattr(mod, name):
        pytest.fail(f"not ported yet: {mod.__name__}.{name}")
    return getattr(mod, name), form.kwargs


def part_declared(slot, variant):
    function, kwargs = _stem_function(slot, variant, "parameters")
    return function(**kwargs)


def part_block(slot, variant):
    """The part function with every own parameter open."""
    function, kwargs = _stem_function(slot, variant, "casadi")
    return function(**kwargs)


def part(slot, variant, values=None):
    """The part function with its numbers checked and frozen in (``values`` None: parameters left open)."""
    from more_transformations.more_casadi_transformations import freeze

    block = part_block(slot, variant)
    if values is None:
        return block
    return freeze(block, part_declared(slot, variant), values)


def skeleton():
    return module(SKELETON_MODULE)


def units_of(forms):
    """``{"drive": ..., "law": ...}`` from the unit constants of the drive and law modules."""
    return {"drive": form_module("drive", forms["drive"]).DRIVE_OUTPUT_UNIT,
            "law": form_module("law", forms["law"]).LAW_INPUT_UNIT}


def propulsor(parts, units, shaft_torque_scale=0.0):
    """A ``Propulsor`` from part functions, the skeleton's number checked and frozen in: a function of ``command,
    state, nu_r, water_density`` and any still-open part parameter."""
    from more_transformations.more_casadi_transformations import freeze

    s = skeleton()
    values = s.check_propulsor_values({"shaft_torque_scale": shaft_torque_scale})
    return freeze(s.propulsor_casadi(parts, units), s.propulsor_parameters(), values)


def propulsor_from_forms(forms, part_values=None, shaft_torque_scale=0.0):
    """``forms`` slot -> variant; ``part_values`` slot -> numbers (absent: the part has none, or its parameters stay
    open when the value is None)."""
    part_values = part_values or {}
    parts = {}
    for slot in SLOTS:
        variant = forms[slot]
        values = part_values.get(slot)
        parts[slot] = part(slot, variant, values if part_declared(slot, variant) else None)
    return propulsor(parts, units_of(forms), shaft_torque_scale)


def expected_declared(slot, variant):
    """The declaration of a variant as the tests expect it, the Wageningen geometry expanded from the existing block."""
    from more_dynamics.models.propeller.wageningen_kt_kq import wageningen_parameters

    expected = {}
    for name, spec in DECLARED[(slot, variant)].items():
        if name == WAGENINGEN:
            for d in wageningen_parameters():
                expected[d.name] = (d.shape, d.unit, d.minimum, d.minimum_exclusive, d.maximum)
        else:
            expected[name] = spec
    return expected


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------
def call(function, **inputs):
    """Numbers out of a CasADi function, every output as a flat array."""
    out = function(**{k: np.asarray(v, dtype=float) for k, v in inputs.items()})
    return {k: np.array(v, dtype=float).reshape(-1) for k, v in out.items()}


def run(function, command, nu_r, water_density, state=None):
    """All outputs of a propulsor function."""
    n = function.size1_in("state")
    state = np.zeros(n) if state is None else state
    return call(function, command=command, state=state, nu_r=nu_r, water_density=water_density)


def tau(function, command, nu_r, water_density, state=None):
    return run(function, command, nu_r, water_density, state)["tau"]


def max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def rng():
    return np.random.default_rng(SEED)


def random_nu_r(generator, n, speed=3.0, rate=1.0):
    return np.hstack([generator.uniform(-speed, speed, (n, 3)), generator.uniform(-rate, rate, (n, 3))])


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# --------------------------------------------------------------------------
# Reference data (all of it in this repository, found beside the tests)
# --------------------------------------------------------------------------
JET_CSV = "jet_nozzle_reference.csv"
JET_CSV_SHA256 = "1df0e56cc053837fe2b8bbc68f5494a25ebba1b454343e25616a5494bda505f7"
JET_PARTS_JSON = "jet_nozzle_test_parts.json"
PUBLISHED_JSON = "waterjet_published_values.json"


def jet_parts():
    """The three jet parameter sets (this library's part names)."""
    return json.loads((DATA_DIR / JET_PARTS_JSON).read_text())["sets"]


def jet_reference():
    """Rows recorded once from the shared repository's jet-nozzle plugin: ``{set: dict of arrays}``."""
    with (DATA_DIR / JET_CSV).open() as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    out = {}
    for name in dict.fromkeys(r["set"] for r in rows):
        chosen = [r for r in rows if r["set"] == name]
        out[name] = SimpleNamespace(
            state=np.array([[float(r["thrust_state"]), float(r["angle_state"])] for r in chosen]),
            command=np.array([[float(r["desired_thrust"]), float(r["desired_normalized_nozzle_angle"])]
                              for r in chosen]),
            output=np.array([[float(r[f"output_{i}"]) for i in range(6)] for r in chosen]),
            dot=np.array([[float(r["dot_0"]), float(r["dot_1"])] for r in chosen]))
    return out


def jet_forms_his():
    """The composition that is the shared repository's jet nozzle: his forms in every socket."""
    return {"drive": "normalized_lag", "steering": "jet_nozzle_clamp", "mounting": "jet_nozzle_direction",
            "inflow": "none", "law": "normalized_linear", "housing": "none"}


def jet_values_his(numbers):
    return {"drive": {"thrust_rise_time": numbers["thrust_rise_time"], "thrust_fall_time": numbers["thrust_fall_time"]},
            "steering": {k: numbers[k] for k in ("min_angle", "max_angle", "nozzle_velocity")},
            "mounting": {k: numbers[k] for k in ("position", "trim_angle", "yaw_bias")},
            "law": {"max_thrust": numbers["max_thrust"], "thrust_coefficient": numbers["thrust_coefficient"]}}


def jet_forms_fixed():
    """The default jet: his drive and law, our rotation for the direction and our servo for the nozzle angle."""
    return {"drive": "normalized_lag", "steering": "servo_lag_rate_angle", "mounting": "steered_rotation",
            "inflow": "none", "law": "normalized_linear", "housing": "none"}


def jet_values_fixed(numbers):
    """Our servo with the time constant and rate limit his angle-error clamp implies (1 s; his nozzle velocity), the
    angle limit his ``max_angle`` (valid for symmetric limits)."""
    return {"drive": jet_values_his(numbers)["drive"],
            "steering": {"max_deflection": numbers["max_angle"], "max_rate": numbers["nozzle_velocity"],
                         "time_constant": 1.0},
            "mounting": jet_values_his(numbers)["mounting"],
            "law": jet_values_his(numbers)["law"]}


def published():
    """Published numbers of the waterjet sources, each with its page (``waterjet_published_values.json``)."""
    return json.loads((DATA_DIR / PUBLISHED_JSON).read_text())


def published_value(name):
    return float(published()["values"][name]["value"])


def _csv_columns(path):
    header = path.read_text().splitlines()[0].split(",")
    values = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {name: values[:, i] for i, name in enumerate(header)}


def remus_propeller_values():
    """The propeller numbers of the frozen REMUS 100 parameter file (each with its ``remus100.m`` line there)."""
    file = json.loads((VEHICLE_DATA_DIR / "remus100" / "remus100_parameters_consistent.json").read_text())
    return {name: entry["value"] for name, entry in file["parameters"].items()}


def remus_composition(law, open_water_coefficients=None):
    """``(forms, part_values, shaft_torque_scale)`` of the fixed propeller: shaft-speed command in rpm, no steering,
    shaft along the first axis through the CO, speed-magnitude inflow, one of the three open-water laws. Numbers from
    the frozen REMUS 100 parameter file; ``roll_moment_scale`` is the shaft-torque scale (a comparison value there)."""
    v = remus_propeller_values()
    given = open_water_coefficients == "given"
    variant = {"linearized": "open_water_linearized_given" if given else "open_water_linearized",
               "bollard": "open_water_bollard", "full": "open_water_polynomial"}[law]
    forms = {"drive": "speed_command_rpm", "steering": "fixed", "mounting": "fixed_euler", "inflow": "speed_magnitude",
             "law": variant, "housing": "none"}
    law_values = {"diameter": v["propeller_diameter"], "pitch_diameter_ratio": v["pitch_diameter_ratio"],
                  "blade_area_ratio": v["blade_area_ratio"], "blade_count": v["blade_count"],
                  "max_advance_number": v["max_advance_number"]}
    if given:
        law_values["thrust_torque_coefficients"] = v["thrust_torque_coefficients"]
    values = {"drive": {"min_command": -v["max_shaft_speed"], "max_command": v["max_shaft_speed"]},
              "mounting": {"position": [0.0, 0.0, 0.0], "orientation": [0.0, 0.0, 0.0]},
              "inflow": {"thrust_deduction": v["thrust_deduction"], "wake_fraction": v["wake_fraction"]},
              "law": law_values}
    return forms, values, v["roll_moment_scale"], v["water_density"]


def remus_block_values(open_water_coefficients=None):
    """The same numbers under the propeller block's names, to build ``propeller_casadi`` from."""
    v = remus_propeller_values()
    values = {k: v[k] for k in ("propeller_diameter", "max_shaft_speed", "thrust_deduction", "wake_fraction",
                                "pitch_diameter_ratio", "blade_area_ratio", "blade_count", "max_advance_number",
                                "roll_moment_scale", "water_density")}
    values["position"] = [0.0, 0.0, 0.0]
    values["orientation"] = [0.0, 0.0, 0.0]
    if open_water_coefficients == "given":
        values["thrust_torque_coefficients"] = v["thrust_torque_coefficients"]
    return values


def matlab_remus_rows():
    """MATLAB running ``remus100.m`` (frozen): shaft command (rpm, column ``ui3``), relative velocity, ``X_prop``,
    ``K_prop``; the fins are in the wrench but not in these columns."""
    ref = _csv_columns(FORCE_DATA_DIR / "remus100_actuators_mss_current.csv")
    return SimpleNamespace(rpm=ref["ui3"], nu_r=np.column_stack([ref[f"nu_r_{i:02d}"] for i in range(1, 7)]),
                           x_prop=ref["X_prop"], k_prop=ref["K_prop"])


def _template_number(pattern):
    text = (FORCE_DATA_DIR / "generate_differential_thruster_template.m").read_text()
    match = re.search(pattern, text)
    assert match, pattern
    return float(match.group(1))


def two_thruster_values():
    """The two fixed thrusters' numbers, read from the in-repository template generator of the frozen MATLAB sweep
    (``generate_differential_thruster_template.m`` lines 29-37, written out from the MSS file of 2026-02) with the
    gravity of its latitude from ``more_transformations``: coefficients per thruster, lever arm, the thrust limits and
    the shaft-speed limits they imply. Returns ``(block_values, per_thruster)`` with the block's names and the
    propulsor parts' numbers for [left, right]."""
    from more_transformations.ecef_ned_transforms import ECEFNEDtransform

    k_pos = _template_number(r"k_pos\s*=\s*([0-9.]+)/2;") / 2.0   # the file writes k_pos = <number>/2 (one thruster)
    k_neg = _template_number(r"k_neg\s*=\s*([0-9.]+)/2;") / 2.0
    y = _template_number(r"y_pont\s*=\s*([0-9.]+);")
    f_fwd = _template_number(r"n_max\s*=\s*sqrt\(\(0\.5\*([0-9.]+)\*g\)/k_pos\)")
    f_rev = _template_number(r"n_min\s*=\s*-sqrt\(\(0\.5\*([0-9.]+)\*g\)/k_neg\)")
    latitude = _template_number(r"gravity\(deg2rad\(([0-9.]+)\)\)")
    g = ECEFNEDtransform.gravity(np.deg2rad(latitude), unit="rad")
    forward, reverse = 0.5 * f_fwd * g, 0.5 * f_rev * g
    block = {"positive_thrust_coefficients": [k_pos, k_pos], "negative_thrust_coefficients": [k_neg, k_neg],
             "thruster_positions": [[0.0, -y, 0.0], [0.0, y, 0.0]], "thruster_directions": [E_X, E_X],
             "max_forward_thrust": [forward, forward], "max_reverse_thrust": [reverse, reverse]}
    per_thruster = []
    for side in (-1.0, 1.0):
        per_thruster.append({
            "drive": {"min_command": -float(np.sqrt(reverse / k_neg)), "max_command": float(np.sqrt(forward / k_pos))},
            "mounting": {"position": [0.0, side * y, 0.0], "direction": E_X},
            "law": {"positive_thrust_coefficient": k_pos, "negative_thrust_coefficient": k_neg}})
    return block, per_thruster


TWO_THRUSTER_FORMS = {"drive": "speed_command_rad_s", "steering": "fixed", "mounting": "fixed_direction",
                      "inflow": "none", "law": "quadratic_asymmetric", "housing": "none"}


def matlab_two_thruster_rows():
    """The frozen MATLAB sweep of the template thrusters: commanded shaft speeds [left, right] (rad/s), the thrust of
    each, and the resulting ``tau_X``, ``tau_N``."""
    ref = _csv_columns(FORCE_DATA_DIR / "differential_thruster_matlab.csv")
    return SimpleNamespace(n=np.column_stack([ref["n_cmd_left"], ref["n_cmd_right"]]),
                           thrust=np.column_stack([ref["thrust_left"], ref["thrust_right"]]),
                           tau_x=ref["tau_X"], tau_n=ref["tau_N"])


OUTBOARD_FORMS = {"drive": "power_limited_speed_command_rpm", "steering": "servo_no_limits", "mounting": "azimuth",
                  "inflow": "surge_component", "law": "open_water_signed_advance", "housing": "none"}


def outboard_values():
    """``(block_values, part_values, water_density)`` of the steerable outboard of the frozen parameter set, the block
    names mapped to the part names: speed limit -> drive command limits, thrust limit -> drive, propwash factor ->
    wake fraction. The numbers are the frozen set's (test numbers of the Grethe class, not identified)."""
    block = json.loads((FORCE_DATA_DIR / "parameter_sets.json").read_text())["outboard_rpm"]
    parts = {"drive": {"min_command": -block["max_speed"], "max_command": block["max_speed"],
                       "max_power": block["max_power"], "efficiency": block["efficiency"],
                       "thrust_saturation": block["max_thrust"]},
             "mounting": {"position": block["position"]},
             "inflow": {"thrust_deduction": block["thrust_deduction"],
                        "advance_speed_factor": block["advance_speed_factor"],
                        "wake_fraction": block["propwash_factor"]},
             "law": {"diameter": block["propeller_diameter"], "pitch_diameter_ratio": block["pitch_diameter_ratio"],
                     "blade_area_ratio": block["blade_area_ratio"], "blade_count": block["blade_count"],
                     "reverse_thrust_factor": block["reverse_thrust_factor"],
                     "reverse_torque_factor": block["reverse_torque_factor"]}}
    return block, parts, block["water_density"]
