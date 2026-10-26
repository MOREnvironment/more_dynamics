"""Two forms of the hull loads of a monohull that the model file ``grethe.m`` has and the library had not:

* the residual coefficient with the range edge of the regression softened (``residual_coefficient_ramped``, surge form
  ``ittc_residual_ramped``): equal to the regression's own coefficient from ``Fn_V = 0.6`` upward, zero at rest, a smoothstep
  of ``Fn_V / 0.6`` below the edge, continuous at it; its gate against the model file is the ``ramp`` set of
  ``tests/vehicle_models/test_monohull_grethe.py`` (derivative and runs);
* the coupled sway-yaw law with the two diagonal terms computed from the time constants
  (``manoeuvring_damping = "time_constants_coupled"``): ``Y_v = -M22 / T_sway``, ``N_r = -M66 / T_yaw``, ``Y_r`` and ``N_v``
  given; equal to the typed coupled law (``linear_coupled``) with those four numbers, so the diagonal terms are not typed
  again (rule 16).

Oracles are written here from the stated equations with numpy, not taken from the parts. G1 against the model file is in
``test_monohull_grethe.py``.

Author:    Enio Krizman
Date:      2026-10-26
"""
import numpy as np
import pytest

from more_dynamics.models.damping.hydrodynamic_load_parts import (coupled_manoeuvring_damping_surface,
                                                                  coupled_manoeuvring_damping_surface_parameters)
from more_dynamics.models.surge_resistance import ra14_residual_resistance as ra14

DISPLACED_VOLUME = 0.7863  # m^3, a 806 kg hull in sea water
LENGTH = 5.2
UNIT_SPEED = np.sqrt(9.81 * DISPLACED_VOLUME ** (1.0 / 3.0))  # the speed at Fn_V = 1


def _value(function, speed):
    return float(function(DISPLACED_VOLUME, LENGTH, speed))


def test_the_ramped_residual_equals_the_regression_residual_from_the_range_edge_upward():
    for froude in (0.6, 0.65, 0.8, 1.0):
        speed = froude * UNIT_SPEED
        assert _value(ra14.residual_coefficient_ramped, speed) == pytest.approx(
            _value(ra14.residual_coefficient, speed), rel=1e-13)


def test_the_ramped_residual_is_zero_at_rest_and_a_smoothstep_of_the_froude_number_below_the_edge():
    edge = _value(ra14.residual_coefficient_ramped, ra14.VALIDITY_FROUDE_NUMBER * UNIT_SPEED)
    assert abs(_value(ra14.residual_coefficient_ramped, 0.0)) <= 1e-9 * abs(edge)
    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        speed = q * ra14.VALIDITY_FROUDE_NUMBER * UNIT_SPEED
        assert _value(ra14.residual_coefficient_ramped, speed) == pytest.approx((3 * q**2 - 2 * q**3) * edge, rel=1e-12)


def test_the_ramped_residual_has_no_step_at_the_range_edge_where_the_plain_one_has():
    edge_speed = ra14.VALIDITY_FROUDE_NUMBER * UNIT_SPEED
    below, above = edge_speed * (1.0 - 1e-9), edge_speed * (1.0 + 1e-9)
    assert abs(_value(ra14.residual_coefficient_ramped, above) - _value(ra14.residual_coefficient_ramped, below)) <= 1e-6 * abs(
        _value(ra14.residual_coefficient_ramped, above))
    step = _value(ra14.residual_coefficient, above) - _value(ra14.residual_coefficient, below)
    assert step > 0.0 and abs(step) > 0.1 * abs(_value(ra14.residual_coefficient_ramped, above))  # the step it removes


def _mass_and_restoring(seed):
    g = np.random.default_rng(seed)
    a = g.normal(size=(6, 6))
    mass = a @ a.T + 6.0 * np.eye(6)  # symmetric positive definite
    restoring = np.diag(np.r_[0.0, 0.0, g.uniform(5e3, 2e4), g.uniform(5e3, 3e4), g.uniform(5e3, 3e4), 0.0])
    return mass, restoring


def test_the_time_constants_coupled_law_equals_the_typed_coupled_law_with_the_same_four_numbers():
    typed = coupled_manoeuvring_damping_surface()
    computed = coupled_manoeuvring_damping_surface(from_time_constants=True)
    assert [d.name for d in coupled_manoeuvring_damping_surface_parameters(from_time_constants=True)] == [
        "max_forward_thrust", "max_speed", "damping_ratios", "time_constants", "sway_yaw_coupling"]
    g = np.random.default_rng(3)
    worst = 0.0
    for seed in range(40):
        mass, restoring = _mass_and_restoring(seed)
        nu_r = g.uniform(-2.0, 2.0, 6)
        time_constants, coupling = g.uniform(1.0, 20.0, 2), g.uniform(-3000.0, 3000.0, 2)
        common = dict(nu_r=nu_r, mass_matrix=mass, restoring_matrix=restoring, max_forward_thrust=0.0, max_speed=2.4,
                      damping_ratios=[0.3, 0.05, 0.5])
        derived = np.array([-mass[1, 1] / time_constants[0], coupling[0], coupling[1], -mass[5, 5] / time_constants[1]])
        out_typed = typed(**common, sway_yaw_derivatives=derived)
        out_computed = computed(**common, time_constants=time_constants, sway_yaw_coupling=coupling)
        for name in ("tau", "D", "damping_derivatives"):
            worst = max(worst, float(np.max(np.abs(np.asarray(out_typed[name]) - np.asarray(out_computed[name])))))
    assert worst <= 1e-12 * 1e5, worst


def test_the_time_constants_coupled_law_follows_the_time_constants_it_is_given():
    """Perturbation: halving T_yaw doubles ``-N_r`` (the diagonal term is computed, not read from elsewhere)."""
    computed = coupled_manoeuvring_damping_surface(from_time_constants=True)
    mass, restoring = _mass_and_restoring(1)
    common = dict(nu_r=np.zeros(6), mass_matrix=mass, restoring_matrix=restoring, max_forward_thrust=0.0, max_speed=2.4,
                  damping_ratios=[0.3, 0.05, 0.5], sway_yaw_coupling=[0.0, 2500.0])
    slow = np.asarray(computed(**common, time_constants=[8.0, 16.0])["damping_derivatives"]).ravel()
    fast = np.asarray(computed(**common, time_constants=[8.0, 8.0])["damping_derivatives"]).ravel()
    assert fast[5] == pytest.approx(2.0 * slow[5], rel=1e-14) and fast[1] == slow[1]  # N_r changes, Y_v does not
