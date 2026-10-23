"""The propeller, the two fixed thrusters and the steerable outboard as compositions of the propulsor.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

Each composition equals the existing block of its device (the second
implementation, G2) and, where a frozen MATLAB row set exists, MATLAB running
MSS (G1):

* Fixed propeller (shaft-speed command in rpm, no steering, shaft along the
  first body axis through the CO, speed-magnitude inflow, shaft torque scale):
  the three open-water laws equal ``propeller_casadi`` with the same numbers on
  random states beyond the saturation limits (``G2_TOLERANCE``); the
  linearised law with the printed coefficients equals the frozen MATLAB
  ``X_prop`` and ``K_prop`` rows of ``remus100.m`` (``tau`` = [(1-t) X_prop, 0,
  0, K_prop/10, 0, 0], ``remus100.m`` 249-252, ``G1_TOLERANCE``).
* Two fixed thrusters (shaft-speed command in rad/s, fixed direction, quadratic
  asymmetric law, no inflow): two propulsor instances, their wrenches summed
  by the vehicle, equal ``differential_thruster_casadi`` and the frozen
  MATLAB sweep of ``tau_X``, ``tau_N`` (the ``otter.m`` thrust law written
  out in the template generator, whose numbers the contract reads).
* Steerable outboard (rpm command, direct steering angle, azimuth mounting at
  a point, surge-component inflow, signed-advance Wageningen law with reverse
  factors, shaft-power cap and thrust limit in the drive): equal to
  ``outboard_motor_rpm_casadi(power_limit=True)`` on random states in which
  the power cap acts in many and the thrust limit in some.

Generality: each device is also built with a second parameter set
that shares none of the first set's numbers; the equality must hold there too.
Second sets are test constructions (arbitrary numbers, no craft's).

Author:    Enio Krizman
Date:      2026-10-10
"""

import numpy as np
import pytest

from propulsor_contract import (G1_TOLERANCE, G2_TOLERANCE, N_RANDOM, OUTBOARD_FORMS, TWO_THRUSTER_FORMS, call,
                                matlab_remus_rows, matlab_two_thruster_rows, max_diff, outboard_values,
                                propulsor, propulsor_from_forms, remus_block_values, remus_composition,
                                random_nu_r, rng, run, tau, two_thruster_values, units_of, part)

from more_dynamics.models.outboard.outboard_motor import outboard_motor_rpm_casadi, outboard_motor_rpm_parameters
from more_dynamics.models.propeller.propeller import propeller_casadi, propeller_parameters
from more_dynamics.models.thruster.differential_thruster import (check_differential_thruster_values,
                                                                  differential_thruster_casadi,
                                                                  differential_thruster_parameters)
from more_transformations.more_casadi_transformations import check_values, freeze

# a second fixed propeller sharing nothing with the first: other geometry, other mounting, tilted shaft
SECOND_PROPELLER = {"propeller_diameter": 0.3, "max_shaft_speed": 800.0, "thrust_deduction": 0.15, "wake_fraction": 0.1,
                    "pitch_diameter_ratio": 0.9, "blade_area_ratio": 0.55, "blade_count": 4, "max_advance_number": 0.8,
                    "roll_moment_scale": 0.2, "water_density": 1025.0, "position": [-0.8, 0.0, 0.05],
                    "orientation": [0.0, 0.1, 0.05]}


def _propeller_block(values, open_water_model, coefficients=None):
    declared = propeller_parameters(open_water_coefficients="given" if coefficients else "polynomial")
    return freeze(propeller_casadi(open_water_model=open_water_model,
                                   open_water_coefficients="given" if coefficients else "polynomial"),
                  declared, check_values(declared, values))


def _propeller_composition(values, law, coefficients=None):
    variant = {"linearized": "open_water_linearized_given" if coefficients else "open_water_linearized",
               "bollard": "open_water_bollard", "full": "open_water_polynomial"}[law]
    forms = {"drive": "speed_command_rpm", "steering": "fixed", "mounting": "fixed_euler",
             "inflow": "speed_magnitude", "law": variant, "housing": "none"}
    law_values = {k: values[k] for k in ("pitch_diameter_ratio", "blade_area_ratio", "blade_count",
                                         "max_advance_number")}
    law_values["diameter"] = values["propeller_diameter"]
    if coefficients:
        law_values["thrust_torque_coefficients"] = coefficients
    part_values = {"drive": {"min_command": -values["max_shaft_speed"], "max_command": values["max_shaft_speed"]},
                   "mounting": {"position": values["position"], "orientation": values["orientation"]},
                   "inflow": {"thrust_deduction": values["thrust_deduction"], "wake_fraction": values["wake_fraction"]},
                   "law": law_values}
    return propulsor_from_forms(forms, part_values, values["roll_moment_scale"])


def _remus_values():
    return remus_block_values()


# --------------------------------------------------------------------------
# The fixed propeller
# --------------------------------------------------------------------------
@pytest.mark.parametrize("which", ["first", "second"])
@pytest.mark.parametrize("law", ["linearized", "bollard", "full"])
def test_G2_the_fixed_propeller_equals_the_propeller_block_in_each_law(law, which):
    values = _remus_values() if which == "first" else SECOND_PROPELLER
    block, composition = _propeller_block(values, law), _propeller_composition(values, law)
    g = rng()
    limit = values["max_shaft_speed"]
    for n, nu_r in zip(g.uniform(-1.2 * limit, 1.2 * limit, N_RANDOM), random_nu_r(g, N_RANDOM)):
        expected = call(block, n=[n], nu_r=nu_r)["tau"]
        got = tau(composition, [n], nu_r, values["water_density"])
        assert max_diff(got, expected) <= G2_TOLERANCE * max(1.0, float(np.max(np.abs(expected)))), (law, n)


def test_G2_the_fixed_propeller_with_given_coefficients_equals_the_block_with_them():
    values = _remus_values()
    coefficients = [0.4566, 0.07, 0.1798, 0.0312]                       # the numbers the frozen REMUS file prints
    from propulsor_contract import remus_propeller_values
    assert list(remus_propeller_values()["thrust_torque_coefficients"]) == coefficients
    block = _propeller_block({**values, "thrust_torque_coefficients": coefficients}, "linearized", coefficients)
    composition = _propeller_composition(values, "linearized", coefficients)
    g = rng()
    for n, nu_r in zip(g.uniform(-1800, 1800, N_RANDOM), random_nu_r(g, N_RANDOM)):
        expected = call(block, n=[n], nu_r=nu_r)["tau"]
        assert max_diff(tau(composition, [n], nu_r, values["water_density"]), expected) <= G2_TOLERANCE * max(
            1.0, float(np.max(np.abs(expected))))


def test_G1_the_fixed_propeller_equals_the_matlab_rows_of_the_linearised_law():
    forms, part_values, scale, rho = remus_composition("linearized", open_water_coefficients="given")
    f = propulsor_from_forms(forms, part_values, scale)
    ref = matlab_remus_rows()
    t_prop = part_values["inflow"]["thrust_deduction"]
    worst = 0.0
    for rpm, nu_r, x_prop, k_prop in zip(ref.rpm, ref.nu_r, ref.x_prop, ref.k_prop):
        got = tau(f, [rpm], nu_r, rho)
        expected = [(1.0 - t_prop) * x_prop, 0.0, 0.0, scale * k_prop, 0.0, 0.0]   # remus100.m 249-252
        worst = max(worst, max_diff(got, expected))
    assert worst <= G1_TOLERANCE, worst
    assert len(ref.rpm) == 1036 and np.any(np.abs(ref.rpm) > 1525.0)             # beyond the saturation, both signs


def test_the_shaft_torque_scale_is_the_only_thing_between_a_free_and_a_reacting_shaft():
    values = _remus_values()
    free = _propeller_composition({**values, "roll_moment_scale": 0.0}, "linearized")
    full = _propeller_composition({**values, "roll_moment_scale": 1.0}, "linearized")
    a, b = tau(free, [900.0], [1.5, 0, 0, 0, 0, 0], values["water_density"]), tau(full, [900.0], [1.5, 0, 0, 0, 0, 0],
                                                                                values["water_density"])
    assert max_diff(a[[0, 1, 2, 4, 5]], b[[0, 1, 2, 4, 5]]) <= 1e-12 and a[3] == 0.0 and b[3] > 0.0


# --------------------------------------------------------------------------
# Two fixed thrusters
# --------------------------------------------------------------------------
SECOND_THRUSTERS = {"positive_thrust_coefficients": [0.09, 0.08], "negative_thrust_coefficients": [0.05, 0.04],
                    "thruster_positions": [[-1.2, -0.9, 0.3], [-1.2, 0.9, 0.3]],
                    "thruster_directions": [[1.0, 0.05, 0.0], [1.0, -0.05, 0.0]],
                    "max_forward_thrust": [400.0, 380.0], "max_reverse_thrust": [200.0, 180.0]}


def _thrusters_from_block_values(values):
    per = []
    for i in range(2):
        k_pos, k_neg = values["positive_thrust_coefficients"][i], values["negative_thrust_coefficients"][i]
        per.append({"drive": {"min_command": -float(np.sqrt(values["max_reverse_thrust"][i] / k_neg)),
                              "max_command": float(np.sqrt(values["max_forward_thrust"][i] / k_pos))},
                    "mounting": {"position": values["thruster_positions"][i],
                                 "direction": values["thruster_directions"][i]},
                    "law": {"positive_thrust_coefficient": k_pos, "negative_thrust_coefficient": k_neg}})
    return per


def _pair(per_thruster):
    return [propulsor_from_forms(TWO_THRUSTER_FORMS, values) for values in per_thruster]


def _summed(pair, n):
    return sum(tau(f, [n[i]], np.zeros(6), 1025.0) for i, f in enumerate(pair))


def test_G1_two_fixed_thrusters_equal_the_matlab_sweep():
    block, per = two_thruster_values()
    pair = _pair(per)
    ref = matlab_two_thruster_rows()
    worst = 0.0
    for n, tau_x, tau_n in zip(ref.n, ref.tau_x, ref.tau_n):
        got = _summed(pair, n)
        worst = max(worst, max_diff(got[[0, 5]], [tau_x, tau_n]), max_diff(got[[1, 2, 3, 4]], 0.0))
    assert worst <= G1_TOLERANCE, worst


@pytest.mark.parametrize("which", ["first", "second"])
def test_G2_two_fixed_thrusters_equal_the_differential_thruster_block(which):
    values = two_thruster_values()[0] if which == "first" else SECOND_THRUSTERS
    declared = differential_thruster_parameters()
    block = freeze(differential_thruster_casadi(), declared, check_differential_thruster_values(values))
    pair = _pair(_thrusters_from_block_values(values))
    g = rng()
    n_max = np.sqrt(np.array(values["max_forward_thrust"]) / np.array(values["positive_thrust_coefficients"]))
    for _ in range(N_RANDOM):
        n = g.uniform(-1.3, 1.3, 2) * n_max
        expected = call(block, n=n, nu_r=np.zeros(6))["tau"]
        assert max_diff(_summed(pair, n), expected) <= G2_TOLERANCE * max(1.0, float(np.max(np.abs(expected))))


def test_two_thrusters_take_nothing_from_the_water_or_the_boat_speed():
    pair = _pair(two_thruster_values()[1])
    a = [tau(f, [20.0], np.zeros(6), 1000.0) for f in pair]
    b = [tau(f, [20.0], [3.0, 1.0, 0, 0.2, 0, 0.1], 1030.0) for f in pair]
    assert max_diff(np.array(a), np.array(b)) == 0.0


# --------------------------------------------------------------------------
# The steerable outboard
# --------------------------------------------------------------------------
SECOND_OUTBOARD = {"max_thrust": 1500.0, "max_power": 6000.0, "position": [-3.5, 0.4, 0.7], "max_speed": 2000.0,
                   "propeller_diameter": 0.4, "pitch_diameter_ratio": 1.1, "blade_count": 4,
                   "reverse_thrust_factor": 0.4, "reverse_torque_factor": 0.7}


def _outboard_pair(block_values):
    declared = outboard_motor_rpm_parameters()
    block = freeze(outboard_motor_rpm_casadi(power_limit=True), declared, check_values(declared, block_values))
    parts = {"drive": {"min_command": -block_values["max_speed"], "max_command": block_values["max_speed"],
                       "max_power": block_values["max_power"], "efficiency": block_values["efficiency"],
                       "thrust_saturation": block_values["max_thrust"]},
             "mounting": {"position": block_values["position"]},
             "inflow": {"thrust_deduction": block_values["thrust_deduction"],
                        "advance_speed_factor": block_values["advance_speed_factor"],
                        "wake_fraction": block_values["propwash_factor"]},
             "law": {"diameter": block_values["propeller_diameter"],
                     "pitch_diameter_ratio": block_values["pitch_diameter_ratio"],
                     "blade_area_ratio": block_values["blade_area_ratio"], "blade_count": block_values["blade_count"],
                     "reverse_thrust_factor": block_values["reverse_thrust_factor"],
                     "reverse_torque_factor": block_values["reverse_torque_factor"]}}
    return block, propulsor_from_forms(OUTBOARD_FORMS, parts)


@pytest.mark.parametrize("which", ["first", "second"])
def test_G2_the_steerable_outboard_equals_the_outboard_block_the_cap_and_limit_acting(which):
    first = outboard_values()[0]
    values = first if which == "first" else {**first, **SECOND_OUTBOARD}
    block, composition = _outboard_pair(values)
    g = rng()
    capped = limited = astern = 0
    for _ in range(N_RANDOM):
        n = g.uniform(-1.2, 1.2) * values["max_speed"]
        delta = g.uniform(-0.6, 0.6)
        nu_r = random_nu_r(g, 1)[0]
        expected = call(block, n=[n], delta=[delta], nu_r=nu_r)["tau"]
        out = run(composition, [n, delta], nu_r, values["water_density"])
        assert max_diff(out["tau"], expected) <= G2_TOLERANCE * max(1.0, float(np.max(np.abs(expected)))), (n, delta)
        limited += abs(out["thrust"][0]) >= values["max_thrust"] - 1e-9
        astern += n < 0
    assert astern > 50
    # the cap acts on many of the samples; asserted through the block's own thrust below its uncapped law
    sample = run(composition, [values["max_speed"], 0.0], np.zeros(6), values["water_density"])
    assert sample["thrust"][0] <= values["max_thrust"] + 1e-9


def test_the_outboards_power_cap_and_thrust_limit_each_act_somewhere():
    """Two operating points chosen to separate them: a low-power outboard is power-limited at speed, a high-power one is
    thrust-limited."""
    first = outboard_values()[0]
    for override, label in (({"max_power": 300.0, "max_thrust": 1e6}, "power"), ({"max_power": 1e9, "max_thrust": 40.0}, "thrust")):
        block, composition = _outboard_pair({**first, **override})
        free = _outboard_pair({**first, "max_power": 1e9, "max_thrust": 1e6})[1]
        n = 0.9 * first["max_speed"]
        limited = run(composition, [n, 0.0], np.zeros(6), first["water_density"])["thrust"][0]
        unlimited = run(free, [n, 0.0], np.zeros(6), first["water_density"])["thrust"][0]
        assert limited < unlimited - 1.0, (label, limited, unlimited)
        if label == "thrust":
            assert abs(limited - 40.0) <= 1e-9
