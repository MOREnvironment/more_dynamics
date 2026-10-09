"""Gate tests of the deflection-only fin composition: ``LiftingFin`` with
servo ``ideal``, inflow ``translational``, flow angle ``none``,
interference ``none``, section ``quadratic_drag``; fins grouped by a
``ForceProducerSet``.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).
The skeleton's own equations (the only physics it holds):
``alpha = k_delta delta - K_beta beta``; ``f = 1/2 rho A U^2 (C_L lift_axis -
C_D chord_axis)``; ``tau = [f; r x f]`` (Fossen 2011, eq. 12.226, p. 400).

Gates:

* G1: the set of a rudder (lift along -y) and a stern plane (lift along -z)
  on ``remus100.m``'s own numbers equals MATLAB running MSS (``tau`` of
  ``remus100.m`` with the propeller at rest, lines 238-254) on every row.
* G1: four fins, each with half the MSS area, off the centre line in pairs
  (a test construction: +-0.08 m), under one command map ``[[1, 0], [1, 0],
  [0, 1], [0, 1]]`` (the pairs move together) equal the same MATLAB rows:
  the roll moments of a pair cancel.
* G2: equal to today's block ``fins_casadi`` (default sign convention), the
  reference this composition replaces.
* G4: a stern-plane area 1 % off is detected.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 12,
    eq. 12.226, p. 400.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 98, 109,
    113-114, 179-189, 234-254.

Author:    Enio Krizman
Date:      2026-10-08
"""

import numpy as np

from fin_parts_contract import (
    E_X,
    G1_TOLERANCE,
    G2_TOLERANCE,
    G4_FACTOR,
    N_RANDOM,
    NEG_Y,
    NEG_Z,
    fin_from_forms,
    geometry,
    max_diff,
    mss_columns,
    mss_constant,
    mss_fins,
    producer_set,
    random_nu_r,
    rng,
    tau,
)

DEFLECTION_ONLY = {"servo": "ideal", "inflow": "translational", "flow_angle": "none",
                   "interference": "none", "section": "quadratic_drag"}


def mss_numbers():
    ref = mss_fins()
    return {name: mss_constant(ref, name)
            for name in ("rho", "delta_max", "A_r", "A_s", "CL_delta_r", "CL_delta_s", "x_r", "x_s")}


def deflection_only_fin(position, lift_axis, area, lift_slope, max_deflection):
    return fin_from_forms(DEFLECTION_ONLY,
                          {"servo": {"max_deflection": max_deflection}, "section": {"lift_slope": lift_slope}},
                          geometry(position, lift_axis, area, E_X))


def two_fin_set(n, stern_area_factor=1.0):
    rudder = deflection_only_fin([n["x_r"], 0, 0], NEG_Y, n["A_r"], n["CL_delta_r"], n["delta_max"])
    stern = deflection_only_fin([n["x_s"], 0, 0], NEG_Z, stern_area_factor * n["A_s"], n["CL_delta_s"], n["delta_max"])
    return producer_set([rudder, stern], np.eye(2))


def four_fin_set(n, offset=0.08):
    fins = [deflection_only_fin([n["x_r"], 0, -offset], NEG_Y, n["A_r"] / 2, n["CL_delta_r"], n["delta_max"]),
            deflection_only_fin([n["x_r"], 0, offset], NEG_Y, n["A_r"] / 2, n["CL_delta_r"], n["delta_max"]),
            deflection_only_fin([n["x_s"], -offset, 0], NEG_Z, n["A_s"] / 2, n["CL_delta_s"], n["delta_max"]),
            deflection_only_fin([n["x_s"], offset, 0], NEG_Z, n["A_s"] / 2, n["CL_delta_s"], n["delta_max"])]
    return producer_set(fins, [[1, 0], [1, 0], [0, 1], [0, 1]])


def _worst_against_matlab(fin_set, rho):
    ref = mss_fins()
    commands = np.column_stack([ref["ui1"], ref["ui2"]])
    nu_r, expected = mss_columns(ref, "nu_r", 6), mss_columns(ref, "tau", 6)
    return max(max_diff(tau(fin_set, c, v, rho), e) for c, v, e in zip(commands, nu_r, expected))


def test_G1_two_fin_set_equals_matlab_remus100():
    n = mss_numbers()
    assert _worst_against_matlab(two_fin_set(n), n["rho"]) <= G1_TOLERANCE


def test_G1_four_fins_under_one_command_map_equal_matlab_remus100():
    n = mss_numbers()
    assert _worst_against_matlab(four_fin_set(n), n["rho"]) <= G1_TOLERANCE


def test_G4_stern_plane_area_plus_1_percent_is_detected():
    n = mss_numbers()
    assert _worst_against_matlab(two_fin_set(n, stern_area_factor=1.01), n["rho"]) > G4_FACTOR * G1_TOLERANCE


def test_G2_two_fin_set_equals_todays_fins_block():
    from more_dynamics.models.fin.fins import fins_casadi, fins_parameters
    from more_transformations.more_casadi_transformations import freeze

    n = mss_numbers()
    values = {"rudder_area": n["A_r"], "stern_plane_area": n["A_s"],
              "rudder_lift_coefficient": n["CL_delta_r"], "stern_plane_lift_coefficient": n["CL_delta_s"],
              "rudder_position": n["x_r"], "stern_plane_position": n["x_s"],
              "max_deflection": n["delta_max"], "water_density": n["rho"]}
    block = freeze(fins_casadi(), fins_parameters(), values)
    fin_set = two_fin_set(n)
    g = rng()
    for command, nu_r in zip(g.uniform(-0.6, 0.6, (N_RANDOM, 2)), random_nu_r(g, N_RANDOM)):
        expected = np.array(block(delta=command, nu_r=nu_r)["tau"], dtype=float).reshape(-1)
        assert max_diff(tau(fin_set, command, nu_r, n["rho"]), expected) <= G2_TOLERANCE


def test_deflection_only_set_has_no_state():
    fin_set = two_fin_set(mss_numbers())
    assert fin_set.size1_in("state") == 0
    assert fin_set.name_in()[:4] == ["command", "state", "nu_r", "water_density"]
