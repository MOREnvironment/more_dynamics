"""``saturate``, ``safe_norm`` and ``wrench`` (``models/shared/force_producer_common.py``).

Equation verified: a point force ``f`` at ``r`` (BODY, about the CO) gives the
generalized force ``tau = [f; r x f]`` (Fossen 2011, eq. 12.226, p. 400), units
N and N m. Independent check: the cross product of numpy; a sign case with
``r = (1, 0, 0)``, ``f = (0, 1, 0)`` giving a positive yaw moment
(right-handed axes, x forward, y starboard, z down).
"""
import casadi as ca
import numpy as np
import pytest

from more_dynamics.models.shared.force_producer_common import safe_norm, saturate, wrench


def _eval(expression, *symbols_and_values):
    symbols = [s for s, _ in symbols_and_values]
    values = [v for _, v in symbols_and_values]
    if not symbols:
        return np.asarray(ca.evalf(expression)).ravel()
    return np.asarray(ca.Function("f", symbols, [expression])(*values)).ravel()


def test_wrench_sign_case_yaw_moment_of_a_starboard_force_ahead_of_the_origin():
    tau = _eval(wrench(ca.SX([0.0, 1.0, 0.0]), ca.SX([1.0, 0.0, 0.0])))
    assert np.allclose(tau, [0.0, 1.0, 0.0, 0.0, 0.0, 1.0])


def test_wrench_matches_the_cross_product_on_random_values():
    rng = np.random.default_rng(3)
    f, r = ca.SX.sym("f", 3), ca.SX.sym("r", 3)
    expression = wrench(f, r)
    for _ in range(20):
        fv, rv = rng.normal(size=3), rng.normal(size=3)
        tau = _eval(expression, (f, fv), (r, rv))
        assert np.allclose(tau[:3], fv)
        assert np.allclose(tau[3:], np.cross(rv, fv), atol=1e-12)


def test_wrench_of_a_force_through_the_origin_has_no_moment():
    tau = _eval(wrench(ca.SX([3.0, -1.0, 2.0]), ca.SX([0.0, 0.0, 0.0])))
    assert np.allclose(tau[3:], 0.0)


@pytest.mark.parametrize("x, expected", [(-2.0, -1.0), (0.3, 0.3), (5.0, 1.0)])
def test_saturate_clips_to_the_interval(x, expected):
    s = ca.SX.sym("x")
    assert _eval(saturate(s, -1.0, 1.0), (s, x))[0] == expected


def test_saturate_has_the_order_of_numpy_clip_when_the_bounds_cross():
    s = ca.SX.sym("x")
    for x in (-3.0, 0.0, 3.0):
        assert _eval(saturate(s, 2.0, 1.0), (s, x))[0] == np.clip(x, 2.0, 1.0)


def test_saturate_is_elementwise():
    v = ca.SX.sym("v", 3)
    out = _eval(saturate(v, -1.0, 1.0), (v, [-2.0, 0.5, 2.0]))
    assert np.allclose(out, [-1.0, 0.5, 1.0])


NAN = float("nan")
REFUSAL = "saturate refuses a NaN command"


def test_saturate_mx_refuses_a_nan_scalar_at_evaluation():
    x = ca.MX.sym("x")
    f = ca.Function("f", [x], [saturate(x, -1.0, 1.0)])
    with pytest.raises(RuntimeError, match=REFUSAL):
        f(NAN)


def test_saturate_mx_refuses_a_nan_in_a_vector_at_evaluation():
    v = ca.MX.sym("v", 3)
    f = ca.Function("f", [v], [saturate(v, -1.0, 1.0)])
    for bad in ([NAN, 0.0, 0.0], [0.0, 2.0, NAN], [NAN, NAN, NAN]):
        with pytest.raises(RuntimeError, match=REFUSAL):
            f(bad)


def test_saturate_mx_finite_values_and_infinities_are_unchanged():
    v = ca.MX.sym("v", 4)
    f = ca.Function("f", [v], [saturate(v, [-1.0, -2.0, 0.0, -3.0], [1.0, 2.0, 0.5, 3.0])])
    x = [-2.0, 1.5, 0.25, np.inf]
    out = np.asarray(f(x)).ravel()
    assert np.allclose(out, np.clip(x, [-1.0, -2.0, 0.0, -3.0], [1.0, 2.0, 0.5, 3.0]))
    out = np.asarray(f([-np.inf, 0.0, 0.1, -np.inf])).ravel()
    assert np.allclose(out, [-1.0, 0.0, 0.1, -3.0])


def test_saturate_mx_derivative_is_unchanged_where_defined():
    x = ca.MX.sym("x", 3)
    jac = ca.Function("j", [x], [ca.jacobian(saturate(x, -1.0, 1.0), x)])
    assert np.allclose(np.asarray(jac([0.3, -2.0, 5.0])), np.diag([1.0, 0.0, 0.0]))
    assert np.allclose(np.asarray(jac([-0.5, 0.0, 0.99])), np.eye(3))


def test_saturate_sx_gives_nan_for_a_nan_command_and_is_unchanged_when_finite():
    v = ca.SX.sym("v", 3)
    f = ca.Function("f", [v], [saturate(v, -1.0, 1.0)])
    assert np.allclose(np.asarray(f([-2.0, 0.5, 2.0])).ravel(), [-1.0, 0.5, 1.0])
    out = np.asarray(f([NAN, 0.5, 2.0])).ravel()
    assert np.isnan(out[0]) and np.allclose(out[1:], [0.5, 1.0])


def test_saturate_of_a_number_refuses_nan_with_the_message():
    with pytest.raises(ValueError, match=REFUSAL):
        saturate(NAN, -1.0, 1.0)
    with pytest.raises(ValueError, match=REFUSAL):
        saturate([0.0, NAN], -1.0, 1.0)
    with pytest.raises(ValueError, match=REFUSAL):
        saturate(ca.DM([NAN]), -1.0, 1.0)


def test_saturate_of_a_finite_number_returns_the_clip():
    assert float(saturate(3.0, -1.0, 1.0)) == 1.0
    assert np.allclose(np.asarray(saturate([-2.0, 0.5], -1.0, 1.0)).ravel(), [-1.0, 0.5])


def test_safe_norm_is_the_euclidean_norm():
    v = ca.SX.sym("v", 3)
    assert _eval(safe_norm(v), (v, [3.0, 4.0, 12.0]))[0] == pytest.approx(13.0)


def test_safe_norm_has_a_finite_zero_derivative_at_the_origin():
    v = ca.SX.sym("v", 3)
    jac = ca.Function("j", [v], [ca.jacobian(safe_norm(v), v)])
    at_zero = np.asarray(jac([0.0, 0.0, 0.0]))
    assert np.all(np.isfinite(at_zero)) and np.allclose(at_zero, 0.0)
    # away from the origin the derivative is v / |v|
    at_x = np.asarray(jac([3.0, 4.0, 0.0])).ravel()
    assert np.allclose(at_x, [0.6, 0.8, 0.0])
