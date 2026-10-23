"""The steering parts of the propulsor: command -> deflection and state.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

* ``fixed``: a thruster that does not steer: no command, no state, deflection
  zero.
* The servo (module ``servo.servo``, one device with switches) in
  the steering socket: ``no_limits`` is the "direct" steering (the deflection
  is the command, no state); ``ideal`` saturates the command; ``lag_rate_
  angle`` is the first-order lag with a rate limit and the angle limit on the
  command. Its full behaviour is gated in ``tests/force_producers/fin``; here
  only that the propulsor's socket accepts it with the same equations.
* ``jet_nozzle_clamp``: the nozzle angle of the shared repository's jet
  nozzle: a normalised command in [-1, 1] mapped to ``max_angle`` (ahead) or
  ``min_angle`` (astern), a state clamped to [min, max], its derivative the
  angle error clamped to +-nozzle_velocity. Below the clamp this is a
  first-order lag of 1 s (the error is in rad and used as rad/s), above it a
  rate limit.
* The reduction: our servo with time constant 1 s, rate limit equal to the
  clamp's ``nozzle_velocity`` and angle limit ``max_angle`` equals the clamp
  inside the limits when the limits are symmetric (``min_angle =
  -max_angle``), the state is inside them, and the command is the normalised
  command times ``max_angle`` (rad instead of normalised). Each condition is
  shown to matter.

Author:    Enio Krizman
Date:      2026-10-10
"""

import numpy as np
import pytest

from propulsor_contract import N_RANDOM, SERVO_SETTINGS, call, max_diff, part, rng

CLAMP = {"min_angle": -0.5, "max_angle": 0.6, "nozzle_velocity": 0.3}   # small enough that the clamp bites


def _clamp_oracle(state, command, v):
    angle = np.clip(state, v["min_angle"], v["max_angle"])
    c = np.clip(command, -1.0, 1.0)
    ref = c * v["max_angle"] if c >= 0.0 else -c * v["min_angle"]
    return angle, float(np.clip(ref - angle, -v["nozzle_velocity"], v["nozzle_velocity"]))


def test_fixed_steering_has_no_command_no_state_and_no_deflection():
    f = part("steering", "fixed")
    assert f.size1_in("command") == 0 and f.size1_in("servo_state") == 0 and f.size1_out("servo_state_dot") == 0
    out = call(f, command=[], servo_state=[])
    assert out["deflection"][0] == 0.0


def test_the_servo_with_no_switches_is_a_direct_steering():
    f = part("steering", "servo_no_limits")
    g = rng()
    for c in g.uniform(-3.0, 3.0, 50):
        out = call(f, command=[c], servo_state=[])
        assert out["deflection"][0] == c and out["servo_state_dot"].size == 0


def test_the_ideal_servo_saturates_the_command():
    f = part("steering", "servo_ideal", {"max_deflection": 0.4})
    for c, expected in ((0.1, 0.1), (0.9, 0.4), (-0.9, -0.4)):
        assert call(f, command=[c], servo_state=[])["deflection"][0] == pytest.approx(expected, abs=1e-15)


def test_the_lag_rate_angle_servo_in_the_socket_follows_its_equations():
    numbers = {"max_deflection": 0.35, "max_rate": 0.5, "time_constant": 0.1}
    f = part("steering", "servo_lag_rate_angle", numbers)
    g = rng()
    for s, c in zip(g.uniform(-0.35, 0.35, N_RANDOM), g.uniform(-0.8, 0.8, N_RANDOM)):
        out = call(f, command=[c], servo_state=[s])
        target = np.clip(c, -numbers["max_deflection"], numbers["max_deflection"])
        dot = np.clip((target - s) / numbers["time_constant"], -numbers["max_rate"], numbers["max_rate"])
        assert abs(out["deflection"][0] - s) <= 1e-15 and abs(out["servo_state_dot"][0] - dot) <= 1e-13


def test_the_four_servo_settings_are_the_contracts_combinations():
    assert sorted(SERVO_SETTINGS) == ["ideal", "lag_rate_angle", "no_limits", "output_saturated"]


# --------------------------------------------------------------------------
# jet_nozzle_clamp
# --------------------------------------------------------------------------
def test_jet_nozzle_clamp_follows_the_clamped_angle_error():
    f = part("steering", "jet_nozzle_clamp", CLAMP)
    g = rng()
    pos = neg = beyond = 0
    for s, c in zip(g.uniform(-0.8, 0.9, N_RANDOM), g.uniform(-1.4, 1.4, N_RANDOM)):
        out = call(f, command=[c], servo_state=[s])
        angle, dot = _clamp_oracle(s, c, CLAMP)
        assert abs(out["deflection"][0] - angle) <= 1e-15
        assert abs(out["servo_state_dot"][0] - dot) <= 1e-13, (s, c)
        pos += c > 0
        neg += c < 0
        beyond += abs(dot) == CLAMP["nozzle_velocity"]
    assert pos > 50 and neg > 50 and beyond > 5


def test_jet_nozzle_clamp_maps_the_normalised_command_to_each_limit():
    f = part("steering", "jet_nozzle_clamp", CLAMP)
    # state at the target: no motion; a full command ahead / astern asks for max / min angle
    assert call(f, command=[1.0], servo_state=[CLAMP["max_angle"]])["servo_state_dot"][0] == pytest.approx(0.0, abs=1e-15)
    assert call(f, command=[-1.0], servo_state=[CLAMP["min_angle"]])["servo_state_dot"][0] == pytest.approx(0.0, abs=1e-15)
    assert call(f, command=[7.0], servo_state=[CLAMP["max_angle"]])["servo_state_dot"][0] == pytest.approx(0.0, abs=1e-15)


# --------------------------------------------------------------------------
# The reduction: our servo equals the clamp under stated conditions
# --------------------------------------------------------------------------
def _servo_like_clamp(clamp):
    return part("steering", "servo_lag_rate_angle",
                {"max_deflection": clamp["max_angle"], "max_rate": clamp["nozzle_velocity"], "time_constant": 1.0})


def test_our_servo_equals_the_clamp_inside_the_limits_under_the_stated_conditions():
    clamp = {"min_angle": -0.6, "max_angle": 0.6, "nozzle_velocity": 2.0}      # condition 1: symmetric limits
    ours, his = _servo_like_clamp(clamp), part("steering", "jet_nozzle_clamp", clamp)
    g = rng()
    worst = 0.0
    for s, c in zip(g.uniform(-0.6, 0.6, N_RANDOM), g.uniform(-1.0, 1.0, N_RANDOM)):  # condition 2: state inside
        a = call(ours, command=[c * clamp["max_angle"]], servo_state=[s])        # condition 3: command in rad
        b = call(his, command=[c], servo_state=[s])
        worst = max(worst, max_diff(a["deflection"], b["deflection"]), max_diff(a["servo_state_dot"], b["servo_state_dot"]))
    assert worst <= 1e-13, worst


def test_the_reduction_needs_the_one_second_time_constant():
    clamp = {"min_angle": -0.6, "max_angle": 0.6, "nozzle_velocity": 2.0}
    his = part("steering", "jet_nozzle_clamp", clamp)
    other = part("steering", "servo_lag_rate_angle",
                 {"max_deflection": 0.6, "max_rate": 2.0, "time_constant": 0.5})
    a = call(other, command=[0.3], servo_state=[0.0])["servo_state_dot"][0]
    b = call(his, command=[0.3 / 0.6], servo_state=[0.0])["servo_state_dot"][0]
    assert abs(a - b) > 0.1


def test_the_reduction_fails_for_asymmetric_limits():
    clamp = {"min_angle": -0.3, "max_angle": 0.6, "nozzle_velocity": 2.0}
    his = part("steering", "jet_nozzle_clamp", clamp)
    ours = _servo_like_clamp(clamp)
    # full astern: the clamp aims at min_angle = -0.3; our servo, given -max_angle in rad, aims at -0.6 (limit 0.6)
    a = call(ours, command=[-clamp["max_angle"]], servo_state=[0.0])["servo_state_dot"][0]
    b = call(his, command=[-1.0], servo_state=[0.0])["servo_state_dot"][0]
    assert abs(a - b) > 0.1


def test_the_reduction_needs_the_rate_limit_to_be_the_nozzle_velocity():
    clamp = {"min_angle": -0.6, "max_angle": 0.6, "nozzle_velocity": 0.05}     # the clamp bites at 0.05 rad/s
    his = part("steering", "jet_nozzle_clamp", clamp)
    wrong_rate = part("steering", "servo_lag_rate_angle", {"max_deflection": 0.6, "max_rate": 2.0, "time_constant": 1.0})
    a = call(wrong_rate, command=[0.6], servo_state=[0.0])["servo_state_dot"][0]
    b = call(his, command=[1.0], servo_state=[0.0])["servo_state_dot"][0]
    assert abs(a - b) > 0.1


def test_the_output_saturated_servo_limits_the_output_not_the_state_and_winds_up():
    numbers = {"max_deflection": 0.35, "max_rate": 0.5, "time_constant": 0.1}
    f = part("steering", "servo_output_saturated", numbers)
    g = rng()
    for s, c in zip(g.uniform(-0.6, 0.6, N_RANDOM), g.uniform(-0.9, 0.9, N_RANDOM)):
        out = call(f, command=[c], servo_state=[s])
        assert abs(out["deflection"][0] - np.clip(s, -0.35, 0.35)) <= 1e-15
        assert abs(out["servo_state_dot"][0] - np.clip((c - s) / 0.1, -0.5, 0.5)) <= 1e-13   # command not clamped first
    # wound up past the stop: the output is at the stop while the state keeps moving
    out = call(f, command=[-0.3], servo_state=[0.5])
    assert out["deflection"][0] == 0.35 and out["servo_state_dot"][0] < 0.0
