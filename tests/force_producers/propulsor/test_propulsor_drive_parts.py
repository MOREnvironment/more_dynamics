"""The drive parts of the propulsor: command -> shaft output, state, thrust scale, thrust limit.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

* ``speed_command`` (rpm and rad/s variants): the command saturated at its
  limits and converted to rev/s; no state, scale 1, no thrust limit
  (remus100.m 115 for rpm; otter.m 220 for rad/s; conversion 1 rpm = 1/60
  rev/s, 1 rad/s = 1/(2 pi) rev/s). The load torque is ignored.
* ``power_limited_speed_command``: the same speed, plus a cap on the thrust
  when the shaft power ``|2 pi n Q|`` exceeds ``eta P_max`` (the efficiency
  floored at 1e-3), scale ``eta P_max / P`` (otherwise 1), and a thrust limit
  that is the declared saturation. A construction of the steerable outboard
  block (no MSS file holds it), stated there line by line.
* ``normalized_lag``: a normalised-effort state clamped to [-1, 1], pulled to
  the clamped command with separate rise and fall time constants; the output
  is the clamped state (unit "1").

Oracles are written here from the stated equations with numpy, not taken from
the parts.

Author:    Enio Krizman
Date:      2026-10-10
"""

import numpy as np
import pytest

from propulsor_contract import N_RANDOM, call, max_diff, part, rng

LIMITS = {"min_command": -1500.0, "max_command": 1200.0}


def _drive(variant, values=None, **_):
    return part("drive", variant, values if values is not None else LIMITS)


# --------------------------------------------------------------------------
# speed_command
# --------------------------------------------------------------------------
@pytest.mark.parametrize("variant, per_unit", [("speed_command_rpm", 1.0 / 60.0),
                                               ("speed_command_rad_s", 1.0 / (2.0 * np.pi))])
def test_speed_command_saturates_and_converts_to_revolutions_per_second(variant, per_unit):
    f = _drive(variant)
    g = rng()
    for c, q in zip(g.uniform(-2500.0, 2500.0, N_RANDOM), g.uniform(-5.0, 5.0, N_RANDOM)):
        out = call(f, command=[c], drive_state=[], shaft_torque=q)
        expected = np.clip(c, LIMITS["min_command"], LIMITS["max_command"]) * per_unit
        assert abs(out["drive_output"][0] - expected) <= 1e-14 * max(1.0, abs(expected)), (c, q)
        assert out["drive_state_dot"].size == 0
        assert out["thrust_scale"][0] == 1.0
        assert np.isinf(out["thrust_limit"][0]) and out["thrust_limit"][0] > 0


@pytest.mark.parametrize("variant, command", [("speed_command_rpm", 60.0), ("speed_command_rad_s", 2.0 * np.pi)])
def test_one_revolution_per_second_in_either_unit(variant, command):
    out = call(_drive(variant, {"min_command": -1e3, "max_command": 1e3}), command=[command], drive_state=[],
               shaft_torque=0.0)
    assert abs(out["drive_output"][0] - 1.0) <= 1e-14


def test_speed_command_limits_are_asymmetric_when_given_so():
    f = _drive("speed_command_rad_s", {"min_command": -100.0, "max_command": 50.0})
    assert call(f, command=[1e6], drive_state=[], shaft_torque=0.0)["drive_output"][0] == pytest.approx(50.0 / (2 * np.pi))
    assert call(f, command=[-1e6], drive_state=[], shaft_torque=0.0)["drive_output"][0] == pytest.approx(-100.0 / (2 * np.pi))


# --------------------------------------------------------------------------
# power_limited_speed_command
# --------------------------------------------------------------------------
POWER = {**LIMITS, "max_power": 2000.0, "efficiency": 0.6, "thrust_saturation": 520.0}


def _power_oracle(command, torque, values):
    n = np.clip(command, values["min_command"], values["max_command"]) / 60.0
    power = abs(2.0 * np.pi * n * torque)
    limit = values["max_power"] * max(values["efficiency"], 1e-3)
    return n, (limit / power if power > limit else 1.0)


def test_power_limited_drive_scales_the_thrust_down_above_the_power_limit_and_not_below():
    f = part("drive", "power_limited_speed_command_rpm", POWER)
    g = rng()
    acted = untouched = 0
    for c, q in zip(g.uniform(-1500.0, 1500.0, N_RANDOM), g.uniform(-60.0, 60.0, N_RANDOM)):
        out = call(f, command=[c], drive_state=[], shaft_torque=q)
        n, scale = _power_oracle(c, q, POWER)
        assert abs(out["drive_output"][0] - n) <= 1e-14 * max(1.0, abs(n))
        assert abs(out["thrust_scale"][0] - scale) <= 1e-13, (c, q)
        assert out["thrust_limit"][0] == POWER["thrust_saturation"]
        assert out["drive_state_dot"].size == 0
        acted += scale < 1.0
        untouched += scale == 1.0
    assert acted > 20 and untouched > 20, (acted, untouched)   # both branches are exercised


def test_power_limited_drive_floors_the_efficiency_at_one_per_mille():
    low = {**POWER, "efficiency": 1e-4}
    f = part("drive", "power_limited_speed_command_rpm", low)
    out = call(f, command=[1000.0], drive_state=[], shaft_torque=20.0)
    n, scale = _power_oracle(1000.0, 20.0, low)
    assert scale < 1.0 and abs(out["thrust_scale"][0] - scale) <= 1e-13
    assert abs(out["thrust_scale"][0] - POWER["max_power"] * 1e-3 / abs(2.0 * np.pi * n * 20.0)) <= 1e-13


# --------------------------------------------------------------------------
# normalized_lag
# --------------------------------------------------------------------------
LAG = {"thrust_rise_time": 0.6, "thrust_fall_time": 0.3}


def test_normalized_lag_is_a_clamped_asymmetric_first_order_lag():
    f = part("drive", "normalized_lag", LAG)
    g = rng()
    rising = falling = 0
    for s, c in zip(g.uniform(-1.4, 1.4, N_RANDOM), g.uniform(-1.4, 1.4, N_RANDOM)):
        out = call(f, command=[c], drive_state=[s], shaft_torque=g.uniform(-3, 3))
        effort, ref = np.clip(s, -1.0, 1.0), np.clip(c, -1.0, 1.0)
        dot = (ref - effort) / (LAG["thrust_rise_time"] if ref > effort else LAG["thrust_fall_time"])
        assert abs(out["drive_output"][0] - effort) <= 1e-15
        assert abs(out["drive_state_dot"][0] - dot) <= 1e-13, (s, c)
        assert out["thrust_scale"][0] == 1.0 and np.isinf(out["thrust_limit"][0])
        rising += ref > effort
        falling += ref < effort
    assert rising > 50 and falling > 50


def test_normalized_lag_has_one_state_and_one_command():
    f = part("drive", "normalized_lag", LAG)
    assert f.size1_in("command") == 1 and f.size1_in("drive_state") == 1 and f.size1_out("drive_state_dot") == 1
