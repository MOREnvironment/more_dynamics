"""The sign of the steering angle of an azimuth-mounted propulsor (the outboard), pinned by a worked example.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

The outboard block's docstring calls the angle "port-positive" while its
formula gives ``f_y = T sin(delta)`` with ``y`` to starboard, so a positive
angle pushes the thrust point to starboard. Which of the two is wrong is
still to be decided; these tests pin the
*formula* through ``STEERING_SIGN`` in the contract, so a decision to flip
the sign is one constant and the one block change that follows it.

* The wrench of a point force: ``tau = [f; r x f]``, so the yaw moment is ``N =
  x f_y - y f_x`` (Fossen 2011, eq. 12.226, p. 400) and, with the thrust point
  aft of the origin (``x_P < 0``), a force to starboard turns the bow to port
  (negative yaw).
* A worked allocation example (it holds for a signed ``x_P``): ``x_P =
  -2.6 m``, requested ``X = 200 N`` and ``N = +100 N m`` need ``T = 203.7 N`` at
  ``delta = -10.9 deg`` with ``f_y = -38.5 N``; the composition applied with
  that thrust and angle returns those forces and that moment. This is the
  allocator-inverts-the-plant check: the inverse of the plant's law.
* Reversal of the angle reverses the sideways force and the yaw moment and
  keeps the surge force.

Author:    Enio Krizman
Date:      2026-10-10
"""

import math

import numpy as np

from propulsor_contract import STEERING_SIGN, max_diff, propulsor_from_forms, run

POSITION = [-2.6, 0.0, 0.55]
FORMS = {"drive": "normalized_lag", "steering": "servo_no_limits", "mounting": "azimuth", "inflow": "none",
         "law": "normalized_linear", "housing": "none"}


def _outboard_with_thrust(thrust):
    """A steered thrust of exactly ``thrust`` N at ``POSITION`` (effort 1 through a linear law)."""
    values = {"drive": {"thrust_rise_time": 1.0, "thrust_fall_time": 1.0}, "mounting": {"position": POSITION},
              "law": {"max_thrust": thrust, "thrust_coefficient": 1.0}}
    return propulsor_from_forms(FORMS, values)


def _apply(thrust, delta):
    out = run(_outboard_with_thrust(thrust), [1.0, delta], np.zeros(6), 1025.0, state=[1.0])
    return out["tau"]


def test_a_positive_angle_pushes_to_starboard_and_turns_the_bow_to_port_when_the_thrust_point_is_aft():
    tau = _apply(200.0, 0.3)
    assert tau[0] > 0.0 and tau[1] * STEERING_SIGN > 0.0
    # N = x f_y - y f_x with y = 0
    assert abs(tau[5] - POSITION[0] * tau[1]) <= 1e-10
    assert tau[5] * STEERING_SIGN < 0.0                              # x_P < 0: starboard force, negative yaw moment


def test_the_forces_are_t_cos_delta_and_signed_t_sin_delta():
    for thrust, delta in ((200.0, 0.3), (150.0, -0.45), (80.0, 1.2)):
        tau = _apply(thrust, delta)
        assert abs(tau[0] - thrust * math.cos(delta)) <= 1e-10
        assert abs(tau[1] - STEERING_SIGN * thrust * math.sin(delta)) <= 1e-10
        assert tau[2] == 0.0
        assert abs(tau[3] - (-POSITION[2]) * tau[1]) <= 1e-10         # roll: y f_z - z f_y = -z f_y


def test_reversing_the_angle_reverses_the_sideways_force_and_the_yaw_moment():
    a, b = _apply(200.0, 0.35), _apply(200.0, -0.35)
    assert abs(a[0] - b[0]) <= 1e-12
    assert max_diff(a[[1, 3, 5]], -b[[1, 3, 5]]) <= 1e-12


def test_a_worked_allocation_is_inverted_by_the_plant():
    """x_P = -2.6 m, X = 200 N, N = +100 N m  =>  f_y = N / x_P, T = hypot(X, f_y), delta = atan2(f_y, X) (signed
    point, N = x f_y - y f_x); the plant applied with that T and delta returns X and N."""
    x_p, x_req, n_req = POSITION[0], 200.0, 100.0
    f_y = n_req / x_p
    thrust, delta = math.hypot(x_req, f_y), math.atan2(f_y, x_req) / STEERING_SIGN
    assert abs(thrust - 203.7) < 0.05 and abs(math.degrees(delta) - (-10.9)) < 0.05 and abs(f_y - (-38.5)) < 0.05   # hand-worked figures
    tau = _apply(thrust, delta)
    assert abs(tau[0] - x_req) <= 1e-9 and abs(tau[1] - f_y) <= 1e-9 and abs(tau[5] - n_req) <= 1e-9
