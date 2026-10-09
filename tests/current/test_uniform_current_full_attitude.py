"""The current taken into the body axes with the full attitude and rotation rate
(``models/current``, ``uniform_current_full_attitude``).

For a current ``v_c^n = [V cos(beta), V sin(beta), w_c]`` constant in NED,
``v_c^b = R(Theta)^T v_c^n`` and ``d/dt v_c^b = -S(omega) v_c^b`` (Fossen 2011,
eqs. 8.138-8.141, 8.157, pp. 221-225; eq. 2.33, Sec. 2.2, pp. 20-27). The
expected values below are the closed forms of those equations for the stated
attitudes (``R_y``, ``R_x`` and ``R_z`` of eq. 2.15-2.18), written out in each
test; the derivative is also checked against a finite difference of the
velocity along the motion of the attitude (``Theta_dot = T(Theta) omega``,
eq. 2.28, p. 25), which needs no formula for the derivative at all.

The two MSS forms are the heading-only forms they reduce to: the velocity is
the same at ``phi = theta = 0``; the acceleration of ``yaw_rate_terms`` is the
same when also ``p = q = 0``, that of ``full_rotation_rate`` when ``phi = theta
= 0``.

Author:    Enio Krizman
Date:      2026-10-09
"""

import casadi as ca
import numpy as np
import pytest

from more_transformations.more_casadi_transformations import EulerNEDBodyTransforms

from more_dynamics.models.current.current import (
    horizontal_current_full_rotation_rate, horizontal_current_yaw_rate_terms, uniform_current_full_attitude)

DEG = np.pi / 180.0


def _full(eta, nu, speed, direction, vertical=0.0):
    out = uniform_current_full_attitude()(eta=eta, nu=nu, current_speed=speed, current_direction=direction,
                                          current_vertical_speed=vertical)
    return np.asarray(out["current_velocity"]).ravel(), np.asarray(out["current_acceleration"]).ravel()


def _pose(phi=0.0, theta=0.0, psi=0.0):
    return [0, 0, 0, phi, theta, psi]


def _rates(p=0.0, q=0.0, r=0.0):
    return [0, 0, 0, p, q, r]


@pytest.mark.parametrize("theta", [-20 * DEG, 20 * DEG])
def test_pitched_craft_in_a_current_from_astern(theta):
    """psi = 0, beta = 0 (``v^n = [V, 0, 0]``): ``R_y(theta)^T v^n = V [cos theta, 0, sin theta]`` (eq. 8.157 with
    eq. 2.17). The heading-only form would give ``[V, 0, 0]``."""
    v, _ = _full(_pose(theta=theta), _rates(), 0.5, 0.0)
    assert v == pytest.approx([0.5 * np.cos(theta), 0.0, 0.5 * np.sin(theta), 0, 0, 0], abs=1e-15)
    assert abs(v[2]) == pytest.approx(0.5 * np.sin(20 * DEG), abs=1e-15)  # 0.171 m/s of heave flow


def test_pitching_craft_current_acceleration():
    """theta = 20 deg, q = 0.1 rad/s: ``-S(omega) v^b`` with ``omega = [0, q, 0]`` is ``[-q v_z, 0, q v_x]``
    (eq. 8.141)."""
    theta, q = 20 * DEG, 0.1
    v, a = _full(_pose(theta=theta), _rates(q=q), 0.5, 0.0)
    assert a == pytest.approx([-q * v[2], 0.0, q * v[0], 0, 0, 0], abs=1e-15)
    assert a[:3] == pytest.approx([-0.5 * q * np.sin(theta), 0.0, 0.5 * q * np.cos(theta)], abs=1e-15)


@pytest.mark.parametrize("phi, p", [(10 * DEG, 0.3), (-10 * DEG, -0.3)])
def test_rolling_craft_in_a_beam_current(phi, p):
    """psi = 0, beta = 90 deg (``v^n = [0, V, 0]``): ``R_x(phi)^T v^n = V [0, cos phi, -sin phi]`` (eq. 8.157 with
    eq. 2.16); the acceleration ``-S(omega) v^b`` with ``omega = [p, 0, 0]`` is ``[0, p v_z, -p v_y]``
    (eq. 8.141)."""
    v, a = _full(_pose(phi=phi), _rates(p=p), 0.5, 0.5 * np.pi)
    assert v == pytest.approx([0.0, 0.5 * np.cos(phi), -0.5 * np.sin(phi), 0, 0, 0], abs=1e-15)
    assert a == pytest.approx([0.0, p * v[2], -p * v[1], 0, 0, 0], abs=1e-15)


def test_vertical_current_enters_with_the_attitude():
    """``w_c`` (down positive) is the third component of ``v^n``; level, it is the heave flow."""
    v, _ = _full(_pose(), _rates(), 0.0, 0.0, vertical=0.2)
    assert v == pytest.approx([0, 0, 0.2, 0, 0, 0], abs=1e-15)
    v, _ = _full(_pose(theta=0.5), _rates(), 0.0, 0.0, vertical=0.2)
    assert v[0] == pytest.approx(-0.2 * np.sin(0.5), abs=1e-15)  # R_y^T [0, 0, w] = w [-sin theta, 0, cos theta]
    assert v[2] == pytest.approx(0.2 * np.cos(0.5), abs=1e-15)


def test_acceleration_is_the_time_derivative_of_the_velocity_along_the_motion():
    """Finite difference of ``nu_c`` along ``Theta_dot = T(Theta) omega`` equals ``nu_c_dot``, over seeded attitudes
    (|theta| < 60 deg), rates and currents. Independent of any formula for the derivative."""
    rng = np.random.default_rng(20261009)
    h = 1e-5
    for _ in range(50):
        theta = np.array([rng.uniform(-1.2, 1.2), rng.uniform(-1.0, 1.0), rng.uniform(-np.pi, np.pi)])
        omega = rng.uniform(-0.5, 0.5, 3)
        speed, direction, vertical = rng.uniform(0.0, 2.0), rng.uniform(-np.pi, np.pi), rng.uniform(-0.5, 0.5)
        euler_rate = np.asarray(EulerNEDBodyTransforms.T_euler_fossen(theta[0], theta[1])).reshape(3, 3) @ omega
        up = _full(_pose(*(theta + h * euler_rate)), _rates(), speed, direction, vertical)[0]
        down = _full(_pose(*(theta - h * euler_rate)), _rates(), speed, direction, vertical)[0]
        _, analytic = _full(_pose(*theta), _rates(*omega), speed, direction, vertical)
        assert (up - down) / (2 * h) == pytest.approx(analytic, abs=1e-8)


def test_reduces_to_the_heading_only_forms_when_the_craft_is_level():
    """phi = theta = 0: the velocity equals the MSS heading-only velocity; the acceleration equals the yaw-rate form
    when p = q = 0 and the rotation-rate form for any p, q, r (the derivative of the heading-only vector is then
    exact: remus100.m 118-122, otter.m 113-118)."""
    rng = np.random.default_rng(7)
    for _ in range(20):
        psi = rng.uniform(-np.pi, np.pi)
        speed, direction, vertical = rng.uniform(0.0, 2.0), rng.uniform(-np.pi, np.pi), rng.uniform(-0.3, 0.3)
        r = rng.uniform(-0.5, 0.5)
        for name, form, rates in (("yaw_rate_terms", horizontal_current_yaw_rate_terms(), _rates(r=r)),
                                  ("full_rotation_rate", horizontal_current_full_rotation_rate(),
                                   _rates(*rng.uniform(-0.5, 0.5, 3)))):
            kwargs = dict(eta=_pose(psi=psi), nu=rates, current_speed=speed, current_direction=direction,
                          current_vertical_speed=vertical)
            mss = form(**kwargs)
            v, a = _full(_pose(psi=psi), rates, speed, direction, vertical)
            assert v == pytest.approx(np.asarray(mss["current_velocity"]).ravel(), abs=1e-14), name
            assert a == pytest.approx(np.asarray(mss["current_acceleration"]).ravel(), abs=1e-14), name


def test_the_two_mss_forms_differ_from_the_full_form_when_the_craft_is_pitched_and_pitching():
    """The shortcuts are not the physics at attitude: at theta = 20 deg and q = 0.1 the yaw-rate form gives zero
    acceleration where the full form gives about 5 cm/s^2 of heave."""
    eta, nu = _pose(theta=20 * DEG), _rates(q=0.1)
    _, a = _full(eta, nu, 0.5, 0.0)
    mss = horizontal_current_yaw_rate_terms()(eta=eta, nu=nu, current_speed=0.5, current_direction=0.0,
                                              current_vertical_speed=0.0)
    assert np.abs(np.asarray(mss["current_acceleration"])).max() == 0.0
    assert abs(a[2]) > 0.04
