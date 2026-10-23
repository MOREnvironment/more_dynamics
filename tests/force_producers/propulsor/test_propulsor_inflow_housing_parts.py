"""The inflow and housing parts of the propulsor.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

* inflow ``none``: advance speed 0, both factors 1.
* inflow ``speed_magnitude``: ``Va = (1 - w) |nu_r[0:3]|``, thrust factor
  ``1 - t``, torque factor 1 (MSS remus100.m 127, 150, 249). The speed has a
  zero (not NaN) derivative at standstill.
* inflow ``surge_component``: ``Va = c_Va (1 - w) u_r`` (signed: astern flow
  gives a negative advance speed), thrust factor ``1 - t`` (the steerable
  outboard block's construction, stated there).
* housing ``none``: zero wrench.

The fractions ``t`` and ``w`` are refused outside [0, 1] by the declaration
(``Parameter`` ranges); the checks are in ``test_propulsor_type_interface``.

Author:    Enio Krizman
Date:      2026-10-10
"""

import casadi as ca
import numpy as np

from propulsor_contract import N_RANDOM, call, max_diff, part, random_nu_r, rng


def _inflow(variant, numbers=None):
    return part("inflow", variant, numbers)


def _flow(f, nu_r):
    return call(f, nu_r=nu_r, thrust_axis=[1.0, 0.0, 0.0], thrust_position=[-1.0, 0.0, 0.0])


def test_no_inflow_has_no_advance_speed_and_unit_factors():
    out = _flow(_inflow("none"), [2.0, 1.0, 0.5, 0.1, 0.2, 0.3])
    assert out["advance_speed"][0] == 0.0 and out["thrust_factor"][0] == 1.0 and out["torque_factor"][0] == 1.0


def test_speed_magnitude_inflow():
    numbers = {"thrust_deduction": 0.1, "wake_fraction": 0.056}
    f = _inflow("speed_magnitude", numbers)
    g = rng()
    for nu_r in random_nu_r(g, N_RANDOM):
        out = _flow(f, nu_r)
        assert abs(out["advance_speed"][0] - (1.0 - 0.056) * np.linalg.norm(nu_r[:3])) <= 1e-14
        assert abs(out["thrust_factor"][0] - 0.9) <= 1e-15 and out["torque_factor"][0] == 1.0


def test_speed_magnitude_has_a_zero_derivative_at_standstill():
    f = _inflow("speed_magnitude", {"thrust_deduction": 0.1, "wake_fraction": 0.056})
    nu_r = ca.SX.sym("nu_r", 6)
    out = f(nu_r=nu_r, thrust_axis=ca.DM([1, 0, 0]), thrust_position=ca.DM([-1, 0, 0]))
    jac = ca.Function("j", [nu_r], [ca.jacobian(out["advance_speed"], nu_r)])(np.zeros(6))
    assert np.all(np.isfinite(np.array(jac))) and float(np.max(np.abs(np.array(jac)))) == 0.0


def test_surge_component_inflow_is_signed():
    numbers = {"thrust_deduction": 0.1, "advance_speed_factor": 0.9, "wake_fraction": 0.25}
    f = _inflow("surge_component", numbers)
    g = rng()
    for nu_r in random_nu_r(g, N_RANDOM):
        out = _flow(f, nu_r)
        assert abs(out["advance_speed"][0] - 0.9 * 0.75 * nu_r[0]) <= 1e-14
        assert abs(out["thrust_factor"][0] - 0.9) <= 1e-15 and out["torque_factor"][0] == 1.0
    assert _flow(f, [-2.0, 0, 0, 0, 0, 0])["advance_speed"][0] < 0.0


def test_surge_component_ignores_the_other_velocity_entries():
    f = _inflow("surge_component", {"thrust_deduction": 0.0, "advance_speed_factor": 1.0, "wake_fraction": 0.0})
    a = _flow(f, [1.5, 0, 0, 0, 0, 0])["advance_speed"][0]
    b = _flow(f, [1.5, 3.0, -2.0, 1.0, 1.0, 1.0])["advance_speed"][0]
    assert a == b == 1.5


def test_no_housing_adds_no_wrench():
    f = part("housing", "none")
    g = rng()
    for nu_r in random_nu_r(g, 20):
        out = call(f, nu_r=nu_r, thrust_axis=[1, 0, 0], thrust_position=[-1, 0, 0], steering_angle=0.3,
                   water_density=1025.0)
        assert max_diff(out["tau"], np.zeros(6)) == 0.0
