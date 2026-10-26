"""The parts of a unit commanded by its thrust: the drives ``thrust_command`` and ``thrust_command_lag``, the inflow
``horizontal_speed_floored`` and the law ``thrust_power_limited`` (unit ``"N"`` between the drive and the law), and the
outboard they make with the steering servo, the azimuth mounting and the leg housing.

Oracles are written here from the stated equations with numpy, not taken from the parts:

* ``thrust_command``: the output is the command, no state, scale 1, no thrust limit;
* ``thrust_command_lag``: the output is the state, ``d(thrust)/dt = (command - thrust) / T``;
* ``horizontal_speed_floored``: ``U = max(sqrt(u_r^2 + v_r^2), floor)``, factors 1, finite derivative at rest;
* ``thrust_power_limited``: ``T = clip(command, -T_lim, T_lim)``, ``T_lim = min(T_max, P / U)``;
* the composition with an ideal servo equals the earlier steerable-outboard throttle block
  (``models/outboard/outboard_motor.outboard_motor_throttle_casadi``) on the wrench, and the lag form changes it only
  through the state.

Author:    Enio Krizman
Date:      2026-10-26
"""
import numpy as np
import pytest

from propulsor_contract import N_RANDOM, call, max_diff, part, propulsor_from_forms, rng, run

TIME_CONSTANT = {"thrust_time_constant": 0.3}
POWER_LAW = {"max_thrust": 511.5, "delivered_power": 1120.0}
FLOOR = {"speed_floor": 0.1}


def _drive(variant, values=None):
    return part("drive", variant, values)


# --------------------------------------------------------------------------
# drive
# --------------------------------------------------------------------------
def test_thrust_command_hands_the_command_on_with_no_state_scale_or_limit():
    f = _drive("thrust_command", {})
    g = rng()
    for c, q in zip(g.uniform(-900.0, 900.0, N_RANDOM), g.uniform(-5.0, 5.0, N_RANDOM)):
        out = call(f, command=[c], drive_state=[], shaft_torque=q)
        assert out["drive_output"][0] == c
        assert out["drive_state_dot"].size == 0 and out["thrust_scale"][0] == 1.0
        assert np.isinf(out["thrust_limit"][0]) and out["thrust_limit"][0] > 0


def test_thrust_command_lag_is_a_first_order_lag_of_the_thrust():
    f = _drive("thrust_command_lag", TIME_CONSTANT)
    g = rng()
    for c, state in zip(g.uniform(-900.0, 900.0, N_RANDOM), g.uniform(-600.0, 600.0, N_RANDOM)):
        out = call(f, command=[c], drive_state=[state], shaft_torque=0.0)
        assert out["drive_output"][0] == state  # the actual thrust is the state
        assert out["drive_state_dot"][0] == pytest.approx((c - state) / 0.3, rel=1e-14, abs=1e-12)
        assert out["thrust_scale"][0] == 1.0 and np.isinf(out["thrust_limit"][0])


@pytest.mark.parametrize("bad", [0.0, -0.1])
def test_thrust_command_lag_refuses_a_time_constant_that_is_not_positive(bad):
    from more_transformations.more_casadi_transformations import check_values
    from propulsor_contract import part_declared
    with pytest.raises(ValueError, match="thrust_time_constant"):
        check_values(part_declared("drive", "thrust_command_lag"), {"thrust_time_constant": bad})


# --------------------------------------------------------------------------
# inflow
# --------------------------------------------------------------------------
def _inflow(nu_r, values=FLOOR):
    f = part("inflow", "horizontal_speed_floored", values)
    return call(f, nu_r=nu_r, thrust_axis=[1, 0, 0], thrust_position=[-2.24, 0, 0.55])


def test_the_horizontal_speed_ignores_heave_and_the_rates_and_is_floored():
    g = rng()
    for nu_r in np.hstack([g.uniform(-3.0, 3.0, (N_RANDOM, 3)), g.uniform(-1.0, 1.0, (N_RANDOM, 3))]):
        out = _inflow(nu_r)
        assert out["advance_speed"][0] == pytest.approx(max(np.hypot(nu_r[0], nu_r[1]), 0.1), rel=1e-14)
        assert out["thrust_factor"][0] == 1.0 and out["torque_factor"][0] == 1.0
    assert _inflow([0.0, 0.0, 5.0, 0.0, 0.0, 2.0])["advance_speed"][0] == 0.1  # at rest the floor stands in


def test_the_horizontal_speed_has_a_finite_derivative_at_rest():
    import casadi as ca
    f = part("inflow", "horizontal_speed_floored", FLOOR)
    nu_r = ca.SX.sym("nu_r", 6)
    out = f(nu_r=nu_r, thrust_axis=ca.SX([1, 0, 0]), thrust_position=ca.SX([-2.24, 0, 0.55]))
    jacobian = ca.Function("j", [nu_r], [ca.jacobian(out["advance_speed"], nu_r)])
    assert np.isfinite(np.asarray(jacobian(np.zeros(6)))).all()


# --------------------------------------------------------------------------
# law
# --------------------------------------------------------------------------
def _law(command, advance, values=POWER_LAW):
    f = part("law", "thrust_power_limited", values)
    return call(f, drive_output=command, advance_speed=advance, water_density=1025.0)


def test_the_thrust_is_the_command_limited_by_the_maximum_thrust_and_the_power_over_the_speed():
    g = rng()
    for command, advance in zip(g.uniform(-900.0, 900.0, N_RANDOM), g.uniform(0.1, 3.0, N_RANDOM)):
        limit = min(511.5, 1120.0 / advance)
        out = _law(command, advance)
        assert out["thrust"][0] == pytest.approx(np.clip(command, -limit, limit), rel=1e-14, abs=1e-12)
        assert out["torque"][0] == 0.0


def test_the_power_cap_acts_above_the_speed_where_it_meets_the_maximum_thrust_and_not_below():
    crossing = 1120.0 / 511.5  # 2.19 m/s
    assert _law(900.0, crossing - 0.2)["thrust"][0] == pytest.approx(511.5)  # the thrust limit
    assert _law(900.0, crossing + 0.5)["thrust"][0] == pytest.approx(1120.0 / (crossing + 0.5))  # the power limit
    assert _law(-900.0, crossing + 0.5)["thrust"][0] == pytest.approx(-1120.0 / (crossing + 0.5))  # symmetric astern


def test_the_law_and_its_drive_pass_a_thrust_in_newtons_and_refuse_a_shaft_speed_drive():
    forms = {"drive": "speed_command_rpm", "steering": "fixed", "mounting": "fixed_direction", "inflow": "none",
             "law": "thrust_power_limited", "housing": "none"}
    with pytest.raises(ValueError, match="units differ"):
        propulsor_from_forms(forms, {"drive": {"min_command": -1.0, "max_command": 1.0},
                                     "mounting": {"position": [0, 0, 0], "direction": [1, 0, 0]},
                                     "law": POWER_LAW})


# --------------------------------------------------------------------------
# the outboard
# --------------------------------------------------------------------------
OUTBOARD = {"drive": "thrust_command", "steering": "servo_ideal", "mounting": "azimuth",
            "inflow": "horizontal_speed_floored", "law": "thrust_power_limited", "housing": "none"}
POSITION = [-2.24, 0.0, 0.55]
DELTA_MAX = 0.6108652381980153


def _outboard(forms=None, extra=None):
    values = {"steering": {"max_deflection": DELTA_MAX}, "mounting": {"position": POSITION},
              "inflow": FLOOR, "law": POWER_LAW, "drive": {}, **(extra or {})}
    return propulsor_from_forms(forms or OUTBOARD, values)


def test_the_outboard_with_an_ideal_servo_equals_the_earlier_throttle_block_on_the_wrench():
    """``outboard_motor_throttle_casadi`` (thrust, angle -> wrench with the power cap), the block the Grethe gates used
    before the outboard was a part: the same wrench for the same thrust command, angle command and flow."""
    from more_dynamics.models.outboard.outboard_motor import outboard_motor_throttle_casadi
    block = outboard_motor_throttle_casadi()
    f = _outboard()
    g = rng()
    worst = 0.0
    for nu_r, thrust, delta in zip(np.hstack([g.uniform(-1.0, 3.0, (N_RANDOM, 2)), g.uniform(-0.3, 0.3, (N_RANDOM, 1)),
                                              g.uniform(-0.3, 0.3, (N_RANDOM, 3))]),
                                   g.uniform(-600.0, 600.0, N_RANDOM), g.uniform(-0.9, 0.9, N_RANDOM)):
        mine = run(f, [thrust, delta], nu_r, 1025.0)
        reference = np.asarray(block(thrust=thrust, delta=np.clip(delta, -DELTA_MAX, DELTA_MAX), nu_r=nu_r,
                                     max_thrust=511.5, max_propulsive_power=1120.0, position=POSITION)["tau"]).ravel()
        worst = max(worst, max_diff(mine["tau"], reference))
    assert worst <= 1e-12, worst


def test_the_steering_angle_command_is_limited_by_the_servo_and_the_thrust_by_the_power():
    f = _outboard()
    out = run(f, [900.0, 2.0], [3.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1025.0)
    assert out["steering_angle"][0] == DELTA_MAX
    assert out["thrust"][0] == pytest.approx(1120.0 / 3.0)


def test_the_lag_form_has_the_thrust_and_the_steering_angle_as_states_and_the_wrench_follows_the_states():
    forms = {**OUTBOARD, "drive": "thrust_command_lag", "steering": "servo_lag_rate_angle"}
    f = _outboard(forms, {"drive": TIME_CONSTANT,
                          "steering": {"max_deflection": DELTA_MAX, "max_rate": 1.0471975511965976,
                                       "time_constant": 0.5}})
    assert f.size1_in("state") == 2
    nu_r = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    out = run(f, [400.0, 0.3], nu_r, 1025.0, state=[200.0, 0.1])
    assert out["state_dot"][0] == pytest.approx((400.0 - 200.0) / 0.3)
    assert out["state_dot"][1] == pytest.approx(min((0.3 - 0.1) / 0.5, 1.0471975511965976))
    same = run(_outboard(), [200.0, 0.1], nu_r, 1025.0)  # the ideal form commanded with the actual values
    assert max_diff(out["tau"], same["tau"]) <= 1e-12
