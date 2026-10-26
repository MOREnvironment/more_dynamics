"""The forms the propulsor gained after its first contract: the normalised-command drive, the steering convention of
the azimuth mounting and the normal force of an outboard leg.

* ``speed_command_normalized``: a signed command in [-1, 1] becomes a shaft speed through the unit's own limits,
  ``n = u n_max`` ahead and ``n = -u n_min`` astern (``n_min`` negative), the inverse of the signed normalisation
  ``u = n / n_max`` ahead and ``n / |n_min|`` astern; output rev/s, no state, scale 1, no thrust limit.
* ``steering_convention`` of the azimuth: ``port_positive`` is the default formula (a positive angle turns the thrust
  force to starboard, so the bow of a boat with the unit aft of the CO to port), ``starboard_positive`` the mirror;
  a thrust point at or ahead of the CO is refused for both (the boat's turn is defined for a unit aft of it).
* ``leg_normal_force`` (housing): the normal force of the outboard leg of the MATLAB model file ``grethe.m``
  (branch ``grethe-vehicle`` @ 01cd570, lines 265-268 and 394-404). G1: the housing's wrench equals the difference
  of the frozen ``tau_out_leg`` and ``tau_out`` columns of ``grethe_m_derivative.csv`` on its 300 states; G3: the
  model file's two ``leg`` runs (turning circle, steering square wave) stepped with RK4 through the vehicle with the
  leg wrench added to the outboard wrench equal the frozen trajectories within the tolerances of
  ``grethe_m_trajectory_tolerances_leg.csv``.
  The leg's area and aspect ratio are the model file's estimates (``A_leg = 0.04`` m^2, ``AR_leg = 1.5``, lines
  265-266), typed from those lines.

Author:    Enio Krizman
Date:      2026-10-23
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from propulsor_contract import N_RANDOM, call, max_diff, part, rng

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vehicle_models"))
import test_monohull_grethe as grethe  # noqa: E402
import vehicle_contract as vc  # noqa: E402

from more_dynamics.models.outboard.outboard_motor import outboard_leg_normal_force  # noqa: E402

A_LEG, AR_LEG = 0.04, 1.5         # grethe.m 265-266 [estimates]
DELTA_MAX = np.deg2rad(35.0)      # grethe.m 260
RHO = 1025.0                      # MSS constants, the density of the frozen runs
G1_TOLERANCE = 1e-9


# --------------------------------------------------------------------------
# speed_command_normalized
# --------------------------------------------------------------------------
NORMALIZED = {"max_shaft_speed": 9.5, "min_shaft_speed": -6.0}


def test_a_normalised_command_maps_to_the_shaft_speed_through_the_limits():
    f = part("drive", "speed_command_normalized", NORMALIZED)
    g = rng()
    for u in g.uniform(-1.5, 1.5, N_RANDOM):
        out = call(f, command=[u], drive_state=[], shaft_torque=g.uniform(-5, 5))
        clipped = float(np.clip(u, -1.0, 1.0))
        expected = clipped * NORMALIZED["max_shaft_speed"] if clipped >= 0 else -clipped * NORMALIZED["min_shaft_speed"]
        assert abs(out["drive_output"][0] - expected) <= 1e-14 * max(1.0, abs(expected)), u
        assert out["drive_state_dot"].size == 0 and out["thrust_scale"][0] == 1.0 and np.isinf(out["thrust_limit"][0])


def test_the_normalised_command_reaches_each_limit_at_plus_and_minus_one_and_is_zero_at_zero():
    f = part("drive", "speed_command_normalized", NORMALIZED)
    speed = lambda u: call(f, command=[u], drive_state=[], shaft_torque=0.0)["drive_output"][0]  # noqa: E731
    assert speed(1.0) == NORMALIZED["max_shaft_speed"] and speed(-1.0) == NORMALIZED["min_shaft_speed"] and speed(0.0) == 0.0
    assert speed(7.0) == speed(1.0) and speed(-7.0) == speed(-1.0)


def test_the_allocators_normalisation_is_inverted_by_the_drive():
    """u = n / n_max ahead and n / |n_min| astern, then the drive: n again."""
    f = part("drive", "speed_command_normalized", NORMALIZED)
    for n in (-6.0, -2.5, 0.0, 3.0, 9.5):
        u = n / NORMALIZED["max_shaft_speed"] if n >= 0 else n / abs(NORMALIZED["min_shaft_speed"])
        assert abs(call(f, command=[u], drive_state=[], shaft_torque=0.0)["drive_output"][0] - n) <= 1e-14


def test_a_positive_astern_limit_is_refused():
    from more_transformations.more_casadi_transformations import freeze
    from propulsor_contract import part_block, part_declared
    with pytest.raises(ValueError, match="min_shaft_speed"):
        freeze(part_block("drive", "speed_command_normalized"), part_declared("drive", "speed_command_normalized"),
               {"max_shaft_speed": 9.5, "min_shaft_speed": 1.0})


# --------------------------------------------------------------------------
# azimuth steering convention
# --------------------------------------------------------------------------
def _azimuth(convention):
    from more_dynamics.models.propulsor.mounting import azimuth
    from more_transformations.more_casadi_transformations import freeze
    return freeze(azimuth.azimuth_casadi(steering_convention=convention),
                  azimuth.azimuth_parameters(steering_convention=convention), {"position": [-2.3, 0.0, 0.55]})


def test_port_positive_is_the_default_formula_and_starboard_positive_its_mirror():
    for a in rng().uniform(-np.pi, np.pi, 30):
        port = call(_azimuth("port_positive"), steering_angle=[a])["thrust_axis"]
        starboard = call(_azimuth("starboard_positive"), steering_angle=[a])["thrust_axis"]
        assert max_diff(port, [np.cos(a), np.sin(a), 0.0]) <= 1e-15
        assert max_diff(starboard, [np.cos(a), -np.sin(a), 0.0]) <= 1e-15
    default = part("mounting", "azimuth", {"position": [-2.3, 0.0, 0.55]})
    assert max_diff(call(default, steering_angle=[0.4])["thrust_axis"], call(_azimuth("port_positive"), steering_angle=[0.4])["thrust_axis"]) == 0.0


@pytest.mark.parametrize("convention", ["port_positive", "starboard_positive"])
@pytest.mark.parametrize("position", [[0.0, 0.0, 0.5], [0.4, 0.0, 0.5]])
def test_a_thrust_point_at_or_ahead_of_the_co_is_refused_for_the_convention(convention, position):
    from more_dynamics.models.propulsor.mounting.azimuth import check_azimuth_values
    with pytest.raises(ValueError, match="position"):
        check_azimuth_values({"position": position}, steering_convention=convention)


def test_an_unknown_convention_is_refused_naming_the_choices():
    from more_dynamics.models.propulsor.mounting.azimuth import azimuth_casadi
    with pytest.raises(ValueError, match="steering_convention must be one of"):
        azimuth_casadi(steering_convention="up_positive")


@pytest.mark.parametrize("convention, sign", [("port_positive", -1.0), ("starboard_positive", 1.0)])
def test_the_boats_turn_for_a_unit_aft_of_the_co_is_what_the_name_says(convention, sign):
    """A positive angle with ``port_positive`` yaws the bow to port (negative N for the unit aft), with
    ``starboard_positive`` to starboard."""
    from more_dynamics.models.propulsor.propulsor import propulsor_casadi, propulsor_parameters
    from more_transformations.more_casadi_transformations import freeze
    from propulsor_contract import SLOTS, run, units_of
    forms = {"drive": "normalized_lag", "steering": "servo_no_limits", "mounting": "azimuth", "inflow": "none",
             "law": "normalized_linear", "housing": "none"}
    numbers = {"drive": {"thrust_rise_time": 1.0, "thrust_fall_time": 1.0},
               "law": {"max_thrust": 100.0, "thrust_coefficient": 1.0}}
    parts = {slot: part(slot, forms[slot], numbers.get(slot)) for slot in SLOTS if slot != "mounting"}
    parts["mounting"] = _azimuth(convention)
    f = freeze(propulsor_casadi({slot: parts[slot] for slot in SLOTS}, units_of(forms)), propulsor_parameters(),
               {"shaft_torque_scale": 0.0})
    out = run(f, [1.0, 0.3], np.zeros(6), RHO, state=[1.0])
    assert out["tau"][5] * sign > 0.0


# --------------------------------------------------------------------------
# leg_normal_force
# --------------------------------------------------------------------------
LEG = {"leg_area": A_LEG, "leg_aspect_ratio": AR_LEG}
POSITION = grethe.json.loads((grethe.DATA / "grethe_m_outboard.json").read_text())["position"]


def _leg_part():
    return part("housing", "leg_normal_force", LEG)


def _leg_wrench(nu, delta):
    """The housing part's wrench at the outboard point (axis not read by this form)."""
    return call(_leg_part(), nu_r=nu, thrust_axis=[1.0, 0.0, 0.0], thrust_position=POSITION,
                steering_angle=[delta], water_density=[RHO])["tau"]


def test_G1_the_leg_wrench_equals_the_model_files_frozen_leg_columns():
    table = grethe._table("grethe_m_derivative.csv")
    nu = grethe._columns(table, "x", 12)[:, :6]
    base = grethe._columns(table, "tau_out", 6)
    with_leg = grethe._columns(table, "tau_out_leg", 6)
    worst, bites = 0.0, 0
    for i in range(len(nu)):
        delta = float(np.clip(table["ui_02"][i], -DELTA_MAX, DELTA_MAX))
        leg = _leg_wrench(nu[i], delta)
        worst = max(worst, max_diff(leg, with_leg[i] - base[i]))
        bites += float(np.max(np.abs(leg))) > 1.0
    assert worst <= G1_TOLERANCE, worst
    assert bites > 50, bites


def test_the_leg_wrench_is_the_library_function_applied_to_the_same_numbers():
    g = rng()
    for _ in range(50):
        nu, delta = g.uniform(-3, 3, 6), g.uniform(-0.6, 0.6)
        direct = outboard_leg_normal_force(nu, delta, POSITION, RHO, A_LEG, AR_LEG)
        assert max_diff(_leg_wrench(nu, delta), np.array(direct, dtype=float).reshape(-1)) <= 1e-12


def test_the_leg_opposes_the_flow_across_the_motor_axis_and_vanishes_for_flow_along_it():
    """The normal force acts against the normal flow component (damping); with the flow along the motor axis the normal
    component is zero and only the guard under the root remains (no force)."""
    along = _leg_wrench([2.0, 0.0, 0, 0, 0, 0], 0.0)
    assert abs(along[1]) <= 1e-12 and abs(along[0]) <= 1e-12
    side = _leg_wrench([0.0, 1.0, 0, 0, 0, 0], 0.0)
    assert side[1] < 0.0                                             # flow to starboard pushes the leg to port
    assert _leg_wrench([0.0, -1.0, 0, 0, 0, 0], 0.0)[1] > 0.0


def test_the_leg_force_scales_with_density_and_area_and_the_lift_slope_of_the_aspect_ratio():
    nu = [1.5, 0.8, 0, 0, 0, 0.1]
    base = call(_leg_part(), nu_r=nu, thrust_axis=[1, 0, 0], thrust_position=POSITION, steering_angle=[0.2],
                water_density=[RHO])["tau"]
    dense = call(_leg_part(), nu_r=nu, thrust_axis=[1, 0, 0], thrust_position=POSITION, steering_angle=[0.2],
                 water_density=[1.1 * RHO])["tau"]
    assert max_diff(dense, 1.1 * base) <= 1e-12
    big = part("housing", "leg_normal_force", {"leg_area": 2 * A_LEG, "leg_aspect_ratio": AR_LEG})
    twice = call(big, nu_r=nu, thrust_axis=[1, 0, 0], thrust_position=POSITION, steering_angle=[0.2], water_density=[RHO])["tau"]
    assert max_diff(twice, 2 * base) <= 1e-12
    slender = part("housing", "leg_normal_force", {"leg_area": A_LEG, "leg_aspect_ratio": 3.0})
    ratio = (6.13 * 3.0 / (3.0 + 2.25)) / (6.13 * AR_LEG / (AR_LEG + 2.25))
    longer = call(slender, nu_r=nu, thrust_axis=[1, 0, 0], thrust_position=POSITION, steering_angle=[0.2], water_density=[RHO])["tau"]
    assert max_diff(longer, ratio * base) <= 1e-12


# (the model file's leg runs are gated on the part itself: tests/vehicle_models/test_monohull_grethe.py, set ``leg``)


# --------------------------------------------------------------------------
# the plugin's literal parameter list
# --------------------------------------------------------------------------
def test_the_plugins_parameter_list_holds_every_number_of_every_variant_once_and_the_selectors():
    """The list is written out in the plugin (a literal name and ``unset()`` each, the form the placeholder scan
    reads); it must equal what the variants declare."""
    import importlib
    vc.plugins()
    module = importlib.import_module("more_dynamics.plugins.force_producers.propulsor")
    names = [d.name for d in module.Propulsor.PARAMETERS]
    selectors = ["drive", "steering", "mounting", "inflow", "law", "housing", "steering_convention"]
    assert names[:len(selectors)] == selectors
    assert names[len(selectors):] == ["shaft_torque_scale", *module._number_names()]
    assert len(names) == len(set(names))
    for slot, variants in module.FORMS.items():
        assert set(variants) >= {v for (s, v) in __import__("propulsor_contract").FORMS if s == slot}
