"""Part-swap reduction tests of the fin skeleton: two compositions that
differ in one socket are equal in the limit that socket's physics names.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).
The tests never reach into a part; they swap parts, the way a user upgrades
a fin.

Arrows:

* inflow ``translational`` -> ``rigid_point`` together with flow angle
  ``none`` -> ``small_angle`` (deflection-only -> effective angle, section
  ``quadratic_drag`` in both): equal in straight motion, ``nu_r = (u, 0, 0,
  0, 0, 0)`` of either sign (Prestero 2001, eq. 4.42, p. 32: the flow angle
  vanishes with ``v``, ``w``, ``p``, ``q``, ``r``); different otherwise.
* section ``quadratic_drag`` <-> ``linear_section`` (no zero-lift drag):
  the same lift; the difference is the drag ``1/2 rho A U^2 lift_slope
  alpha^2`` along ``-chord_axis`` and its moment (MSS ``remus100.m``
  238-239 vs Prestero eq. 4.37, p. 31).
* servo ``lag_rate_angle`` at rest (deflection = saturated command) ->
  ``ideal``: the same force, zero rate.
* section ``lifting_line`` -> ``linear_section`` and interference
  ``slender_body`` -> ``none`` (second fidelity -> effective angle with
  Prestero's section): written, waiting on their sources.

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 238-245.
[Prestero 2001] Prestero, T. (2001). *Verification of a six-degree of
    freedom simulation model for the REMUS autonomous underwater vehicle*.
    MSc thesis, MIT/WHOI. Ch. 4, eqs. 4.37, 4.40-4.43, pp. 31-33.
[Murray-Smith 2016] Murray-Smith, D. J. (2016). Inverse simulation methods
    applied to investigations of actuator nonlinearities in ship steering.
    *Simulation Notes Europe* 26(4), 245-256. Fig. 1, p. 246.

Author:    Enio Krizman
Date:      2026-10-08
"""

import numpy as np
import pytest

from fin_parts_contract import (
    E_X,
    G2_TOLERANCE,
    N_RANDOM,
    NEG_Y,
    NEG_Z,
    WAITING_REASON,
    call,
    fin_from_forms,
    geometry,
    max_diff,
    producer_set,
    random_nu_r,
    rng,
    tau,
)

RHO = 1026.0           # test numbers, not a vehicle's
MAX_DEFLECTION = 0.35
GEOMETRY = [([-0.8, 0.0, 0.0], NEG_Y, 0.0133, 0.5), ([-0.78, 0.0, 0.0], NEG_Z, 0.0127, 0.7)]


def _set(inflow, flow_angle, section="quadratic_drag", zero_lift_drag=0.0):
    fins = []
    for position, lift, area, slope in GEOMETRY:
        values = {"lift_slope": slope}
        if section == "linear_section":
            values["zero_lift_drag"] = zero_lift_drag
        fins.append(fin_from_forms({"servo": "ideal", "inflow": inflow, "flow_angle": flow_angle,
                                    "interference": "none", "section": section},
                                   {"servo": {"max_deflection": MAX_DEFLECTION}, "section": values},
                                   geometry(position, lift, area, E_X)))
    return producer_set(fins, np.eye(2))


def test_effective_angle_equals_deflection_only_in_straight_motion():
    deflection_only = _set("translational", "none")
    effective = _set("rigid_point", "small_angle")
    g = rng()
    for command, u in zip(g.uniform(-0.5, 0.5, (N_RANDOM, 2)), g.uniform(-3, 3, N_RANDOM)):
        nu_r = [u, 0, 0, 0, 0, 0]
        assert max_diff(tau(effective, command, nu_r, RHO), tau(deflection_only, command, nu_r, RHO)) <= G2_TOLERANCE


def test_effective_angle_differs_from_deflection_only_in_general_motion():
    deflection_only = _set("translational", "none")
    effective = _set("rigid_point", "small_angle")
    g = rng()
    worst = max(max_diff(tau(effective, c, v, RHO), tau(deflection_only, c, v, RHO))
                for c, v in zip(g.uniform(-0.5, 0.5, (100, 2)), random_nu_r(g, 100)))
    assert worst > 1e-3


@pytest.mark.parametrize("inflow, flow_angle", [("translational", "none"), ("rigid_point", "small_angle")])
def test_section_swap_changes_the_drag_term_only(inflow, flow_angle):
    """tau(quadratic_drag) - tau(linear_section, no drag) = the wrench of
    -1/2 rho A U^2 lift_slope alpha^2 chord_axis at the fin."""
    quadratic, linear = _set(inflow, flow_angle), _set(inflow, flow_angle, "linear_section")
    singles = [fin_from_forms({"servo": "ideal", "inflow": inflow, "flow_angle": flow_angle,
                               "interference": "none", "section": "quadratic_drag"},
                              {"servo": {"max_deflection": MAX_DEFLECTION}, "section": {"lift_slope": slope}},
                              geometry(position, lift, area, E_X))
               for position, lift, area, slope in GEOMETRY]
    g = rng()
    nus = random_nu_r(g, N_RANDOM)
    nus[:, 0] = np.where(np.abs(nus[:, 0]) < 2.0, np.copysign(2.0, nus[:, 0]), nus[:, 0])  # |u| >= 2 m/s
    for command, nu_r in zip(g.uniform(-0.5, 0.5, (N_RANDOM, 2)), nus):
        expected = np.zeros(6)
        for single, c, (position, lift, area, slope) in zip(singles, command, GEOMETRY):
            out = call(single, command=c, state=np.zeros(0), nu_r=nu_r, water_density=RHO)
            r = np.array(position)
            v_fin = nu_r[0:3] + (np.cross(nu_r[3:6], r) if inflow == "rigid_point" else 0.0)
            u_c, u_n = v_fin @ np.array(E_X), v_fin @ np.array(lift)
            speed_squared = u_c ** 2 if flow_angle == "small_angle" else u_c ** 2 + u_n ** 2
            f = -0.5 * RHO * area * speed_squared * slope * out["angle_of_attack"][0] ** 2 * np.array(E_X)
            expected += np.concatenate([f, np.cross(r, f)])
        diff = tau(quadratic, command, nu_r, RHO) - tau(linear, command, nu_r, RHO)
        assert max_diff(diff, expected) <= G2_TOLERANCE * max(1.0, float(np.max(np.abs(expected))))


def test_lag_rate_angle_at_rest_equals_ideal():
    def one(servo, values):
        return fin_from_forms({"servo": servo, "inflow": "rigid_point", "flow_angle": "small_angle",
                               "interference": "none", "section": "quadratic_drag"},
                              {"servo": values, "section": {"lift_slope": 0.5}},
                              geometry([-0.8, 0, 0], NEG_Y, 0.0133, E_X))
    lagged = one("lag_rate_angle", {"max_deflection": MAX_DEFLECTION, "max_rate": 0.5, "time_constant": 0.1})
    ideal = one("ideal", {"max_deflection": MAX_DEFLECTION})
    g = rng()
    for c, nu_r in zip(g.uniform(-0.5, 0.5, N_RANDOM), random_nu_r(g, N_RANDOM)):
        rest = [float(np.clip(c, -MAX_DEFLECTION, MAX_DEFLECTION))]
        out = call(lagged, command=c, state=rest, nu_r=nu_r, water_density=RHO)
        assert max_diff(out["tau"], tau(ideal, c, nu_r, RHO)) <= G2_TOLERANCE
        assert abs(out["state_dot"][0]) <= G2_TOLERANCE


@pytest.mark.skip(reason=WAITING_REASON)
def test_second_fidelity_fin_reduces_to_the_effective_angle_fin():
    """lifting_line + slender_body -> linear_section + none as the aspect
    ratio grows and the body radius over the semispan vanishes: two swaps,
    each testable alone (the section and interference modules hold the
    one-socket tests)."""
    raise AssertionError("contract not fixed")
