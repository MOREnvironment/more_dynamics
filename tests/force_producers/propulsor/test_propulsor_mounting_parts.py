"""The mounting parts of the propulsor: steering angle -> thrust axis and thrust point (BODY, from the CO).

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

* ``fixed_euler``: the axis is the first column of ``R_zyx(orientation)``
  (Fossen 2011, eq. 2.18, p. 22); the point is the declared position; the
  steering angle is ignored.
* ``fixed_direction``: the axis is the declared direction normalised (any
  length); a zero direction is not a direction.
* ``azimuth``: the axis is ``[cos d, sin d, 0]`` about the body z axis
  (Fossen 2011, Table 12.3, p. 398); a positive angle points the thrust to
  starboard (+y): ``STEERING_SIGN`` in the contract, pinned in
  ``test_propulsor_outboard_sign.py``.
* ``jet_nozzle_direction``: the shared repository's direction ``[cos(a+b)
  cos t, sin(a+b), -sin t]`` for steering angle ``a``, yaw bias ``b`` and
  trim ``t`` (his file, lines 91-95): reproduced as it is, so it is **not a
  unit vector** when ``t`` is not zero and the nozzle is not straight; the
  excess is computed here (up to 0.54 % at t = 0.2 rad).
* ``steered_rotation``: the fixed direction: ``R_z(a + b) R_y(t) e_x =
  [cos(a+b) cos t, sin(a+b) cos t, -sin t]``, a unit vector at every angle;
  equal to the jet-nozzle direction at zero trim (Fossen 2011, eqs.
  2.15-2.18, p. 22).

Oracles use the numpy rotations of ``more_transformations``.

Author:    Enio Krizman
Date:      2026-10-10
"""

import numpy as np
import pytest

from more_transformations.matrix_transforms import MatrixTransforms as T

from propulsor_contract import E_X, N_RANDOM, STEERING_SIGN, call, max_diff, part, rng


def _mount(variant, numbers):
    return part("mounting", variant, numbers)


def _axis(f, angle):
    out = call(f, steering_angle=[angle])
    return out["thrust_axis"], out["thrust_position"]


def test_fixed_euler_is_the_first_column_of_the_zyx_rotation_and_ignores_the_steering_angle():
    g = rng()
    for _ in range(50):
        euler = g.uniform(-1.0, 1.0, 3)
        position = g.uniform(-2.0, 2.0, 3)
        f = _mount("fixed_euler", {"position": position.tolist(), "orientation": euler.tolist()})
        axis, point = _axis(f, g.uniform(-1.0, 1.0))
        assert max_diff(axis, T.Rzyx(euler)[:, 0]) <= 1e-14
        assert max_diff(point, position) == 0.0
        assert abs(np.linalg.norm(axis) - 1.0) <= 1e-14
    assert max_diff(_axis(_mount("fixed_euler", {"position": [0, 0, 0], "orientation": [0, 0, 0]}), 0.4)[0], E_X) <= 1e-15


@pytest.mark.parametrize("direction", [[1.0, 0.0, 0.0], [1.0, 0.05, 0.0], [0.0, 3.0, -4.0], [-2.0, 1.0, 0.5]])
def test_fixed_direction_is_the_normalised_direction(direction):
    f = _mount("fixed_direction", {"position": [0.1, -0.2, 0.3], "direction": direction})
    axis, point = _axis(f, 0.3)
    assert max_diff(axis, np.array(direction) / np.linalg.norm(direction)) <= 1e-15
    assert max_diff(point, [0.1, -0.2, 0.3]) == 0.0


def test_azimuth_points_the_thrust_to_starboard_for_a_positive_angle():
    f = _mount("azimuth", {"position": [-2.3, 0.0, 0.55]})
    g = rng()
    for a in g.uniform(-np.pi, np.pi, 50):
        axis, point = _axis(f, a)
        assert max_diff(axis, [np.cos(a), STEERING_SIGN * np.sin(a), 0.0]) <= 1e-15
        assert max_diff(point, [-2.3, 0.0, 0.55]) == 0.0
    assert _axis(f, 0.3)[0][1] * STEERING_SIGN > 0.0


JET_NUMBERS = {"position": [-1.2, 0.1, 0.3], "trim_angle": 0.2, "yaw_bias": 0.03}


def test_jet_nozzle_direction_is_his_formula_and_not_a_unit_vector_when_trimmed():
    f = _mount("jet_nozzle_direction", JET_NUMBERS)
    g = rng()
    longest = 0.0
    for a in g.uniform(-0.6, 0.6, N_RANDOM):
        axis, point = _axis(f, a)
        b, t = JET_NUMBERS["yaw_bias"], JET_NUMBERS["trim_angle"]
        expected = [np.cos(a + b) * np.cos(t), np.sin(a + b), -np.sin(t)]
        assert max_diff(axis, expected) <= 1e-15
        # |axis|^2 = cos^2(x) cos^2(t) + sin^2(x) + sin^2(t) = 1 + sin^2(t) sin^2(x), with x = a + b: longer than 1
        # except at t = 0 or x = 0 (the excess is the identity, not a fit)
        assert abs(np.linalg.norm(axis) ** 2 - (1.0 + np.sin(t) ** 2 * np.sin(a + b) ** 2)) <= 1e-14
        longest = max(longest, np.linalg.norm(axis) - 1.0)
        assert max_diff(point, JET_NUMBERS["position"]) == 0.0
    assert longest > 3e-3, longest            # about 0.7 % at the end of the nozzle range for 0.2 rad of trim


def test_jet_nozzle_direction_is_a_unit_vector_at_zero_trim():
    f = _mount("jet_nozzle_direction", {**JET_NUMBERS, "trim_angle": 0.0})
    for a in np.linspace(-0.6, 0.6, 25):
        assert abs(np.linalg.norm(_axis(f, a)[0]) - 1.0) <= 1e-15


def test_steered_rotation_is_rz_ry_ex_and_a_unit_vector_at_every_angle():
    f = _mount("steered_rotation", JET_NUMBERS)
    g = rng()
    for a in g.uniform(-1.5, 1.5, N_RANDOM):
        axis, point = _axis(f, a)
        expected = T.Rz(a + JET_NUMBERS["yaw_bias"]) @ T.Ry(JET_NUMBERS["trim_angle"]) @ np.array(E_X)
        assert max_diff(axis, expected) <= 1e-14
        assert abs(np.linalg.norm(axis) - 1.0) <= 1e-14
        assert max_diff(point, JET_NUMBERS["position"]) == 0.0


def test_steered_rotation_equals_the_jet_nozzle_direction_at_zero_trim_exactly():
    zero = {**JET_NUMBERS, "trim_angle": 0.0}
    a, b = _mount("steered_rotation", zero), _mount("jet_nozzle_direction", zero)
    for angle in np.linspace(-0.8, 0.8, 41):
        assert max_diff(_axis(a, angle)[0], _axis(b, angle)[0]) <= 1e-15


def test_steered_rotation_differs_from_the_jet_nozzle_direction_when_trimmed():
    a, b = _mount("steered_rotation", JET_NUMBERS), _mount("jet_nozzle_direction", JET_NUMBERS)
    differences = [max_diff(_axis(a, angle)[0], _axis(b, angle)[0]) for angle in np.linspace(-0.6, 0.6, 25)]
    assert max(differences) > 1e-3                                       # the fix is real
    assert min(np.linalg.norm(_axis(a, angle)[0]) for angle in np.linspace(-0.6, 0.6, 25)) > 1.0 - 1e-14
