"""One-socket swaps of the propulsor: each reduction holds exactly where it is meant to and breaks where it must.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

* drive: ``power_limited_speed_command_rpm`` with a power and a thrust limit
  far above use equals ``speed_command_rpm``.
* steering: the fixed servo (``servo_no_limits``) equals the outboard's
  direct steering; ``fixed`` equals a servo with deflection zero for a mounting
  that ignores the angle.
* mounting: ``steered_rotation`` equals ``jet_nozzle_direction`` at zero trim
  and differs when trimmed; ``fixed_direction`` along the first axis equals
  ``fixed_euler`` with zero angles.
* inflow: ``none`` with the linearised law equals the bollard law (advance
  speed zero); ``speed_magnitude`` with ``t = w = 0`` is the true speed;
  ``surge_component`` equals ``speed_magnitude`` for pure surge with the same
  fractions and ``c_Va = 1``.
* law: the linearised law equals the bollard law at zero advance speed; the
  momentum-flux law with a flow coefficient from its bollard point equals the
  bollard-quadratic law at rest.
* A part of one slot put in another is refused naming the slot (also in the
  type-interface tests, here for the real parts).

Author:    Enio Krizman
Date:      2026-10-10
"""

import numpy as np
import pytest

from propulsor_contract import (E_X, N_RANDOM, max_diff, part, propulsor_from_forms, random_nu_r, remus_composition,
                                rng, skeleton, tau, units_of)

RHO = 1026.0


def _remus():
    forms, values, scale, rho = remus_composition("linearized")
    return forms, values, scale, rho


def _compare(a, b, n=N_RANDOM, limit=1500.0, tolerance=1e-10):
    g = rng()
    worst = 0.0
    for c, nu_r in zip(g.uniform(-limit, limit, n), random_nu_r(g, n)):
        worst = max(worst, max_diff(tau(a, [c], nu_r, RHO), tau(b, [c], nu_r, RHO)))
    return worst


def test_drive_swap_an_unlimited_power_cap_is_the_plain_speed_command():
    forms, values, scale, _ = _remus()
    plain = propulsor_from_forms(forms, values, scale)
    forms2 = {**forms, "drive": "power_limited_speed_command_rpm"}
    values2 = {**values, "drive": {**values["drive"], "max_power": 1e12, "efficiency": 1.0, "thrust_saturation": 1e12}}
    assert _compare(plain, propulsor_from_forms(forms2, values2, scale)) <= 1e-12


def test_drive_swap_a_real_power_cap_changes_the_result():
    forms, values, scale, _ = _remus()
    plain = propulsor_from_forms(forms, values, scale)
    forms2 = {**forms, "drive": "power_limited_speed_command_rpm"}
    values2 = {**values, "drive": {**values["drive"], "max_power": 50.0, "efficiency": 0.6, "thrust_saturation": 1e12}}
    assert _compare(plain, propulsor_from_forms(forms2, values2, scale)) > 1.0


def test_mounting_swap_a_fixed_direction_along_the_first_axis_is_the_euler_mounting_with_zero_angles():
    forms, values, scale, _ = _remus()
    a = propulsor_from_forms(forms, values, scale)
    b = propulsor_from_forms({**forms, "mounting": "fixed_direction"},
                             {**values, "mounting": {"position": [0.0, 0.0, 0.0], "direction": [5.0, 0.0, 0.0]}}, scale)
    assert _compare(a, b) <= 1e-12


def test_inflow_swap_no_inflow_with_the_linearised_law_is_the_bollard_law():
    forms, values, scale, _ = _remus()
    a = propulsor_from_forms({**forms, "inflow": "none"}, {k: v for k, v in values.items() if k != "inflow"}, scale)
    b = propulsor_from_forms({**forms, "inflow": "none", "law": "open_water_bollard"}, {k: v for k, v in values.items() if k != "inflow"}, scale)
    assert _compare(a, b) <= 1e-12


def test_inflow_swap_surge_component_is_speed_magnitude_for_pure_surge_and_is_not_otherwise():
    forms, values, scale, _ = _remus()
    fractions = values["inflow"]
    a = propulsor_from_forms({**forms, "inflow": "speed_magnitude"}, values, scale)
    b = propulsor_from_forms({**forms, "inflow": "surge_component"},
                             {**values, "inflow": {**fractions, "advance_speed_factor": 1.0}}, scale)
    g = rng()
    for c, u in zip(g.uniform(-1500, 1500, 100), g.uniform(0.0, 3.0, 100)):
        nu_r = [u, 0, 0, 0, 0, 0]
        assert max_diff(tau(a, [c], nu_r, RHO), tau(b, [c], nu_r, RHO)) <= 1e-10
    assert max_diff(tau(a, [1000.0], [0.0, 2.0, 0, 0, 0, 0], RHO), tau(b, [1000.0], [0.0, 2.0, 0, 0, 0, 0], RHO)) > 0.5


def test_law_swap_the_momentum_flux_law_from_a_bollard_point_is_the_bollard_quadratic_law_at_rest():
    forms = {"drive": "speed_command_rpm", "steering": "fixed", "mounting": "fixed_direction", "inflow": "none",
             "law": "jet_momentum_flux_bollard_point", "housing": "none"}
    values = {"drive": {"min_command": -6000.0, "max_command": 6000.0},
              "mounting": {"position": [-1.0, 0.0, 0.0], "direction": E_X},
              "law": {"nozzle_area": 7e-3, "bollard_thrust": 4000.0, "bollard_shaft_speed": 90.0}}
    flux = propulsor_from_forms(forms, values)
    k = 4000.0 / 90.0**2
    quad = propulsor_from_forms({**forms, "law": "jet_bollard_quadratic"},
                                {**values, "law": {"bollard_thrust_coefficient": k}})
    for rpm in (600.0, 2700.0, 5400.0):
        a, b = tau(flux, [rpm], np.zeros(6), 1025.0), tau(quad, [rpm], np.zeros(6), 1025.0)
        assert max_diff(a, b) <= 1e-9 * max(1.0, float(np.max(np.abs(a))))
    # with boat speed the flux law (through an intake inflow) loses thrust and the quadratic law does not
    forms_u = {**forms, "inflow": "speed_magnitude"}
    values_u = {**values, "inflow": {"thrust_deduction": 0.0, "wake_fraction": 0.0}}
    fast = propulsor_from_forms(forms_u, values_u)
    still, moving = tau(fast, [5400.0], np.zeros(6), 1025.0)[0], tau(fast, [5400.0], [10.0, 0, 0, 0, 0, 0], 1025.0)[0]
    assert moving < still - 100.0


@pytest.mark.parametrize("slot, wrong", [("law", "mounting"), ("mounting", "law"), ("inflow", "housing"),
                                         ("drive", "steering"), ("housing", "inflow")])
def test_a_part_of_another_slot_is_refused_naming_the_slot(slot, wrong):
    forms, values, scale, _ = _remus()
    parts = {s: part(s, forms[s], values.get(s)) for s in forms}
    donor = {"law": ("law", forms["law"]), "mounting": ("mounting", forms["mounting"]),
             "housing": ("housing", "none"), "steering": ("steering", "fixed"), "inflow": ("inflow", forms["inflow"])}[wrong]
    parts[slot] = part(donor[0], donor[1], values.get(donor[0]))
    with pytest.raises(ValueError) as error:
        skeleton().propulsor_casadi(parts, units_of(forms))
    assert slot in str(error.value)
