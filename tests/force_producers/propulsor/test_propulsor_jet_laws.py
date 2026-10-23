"""The waterjet thrust laws J1-J4 of the law socket: higher-fidelity options of the jet.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).
Which law is the default jet law, and whether a reversing bucket belongs in
the model, are open design choices: every law here is an option beside the
constant law and none is asserted to be the default, and no bucket is
contracted. Every law gives zero thrust for a non-positive shaft speed (a
waterjet needs a bucket for reverse thrust: Varyani and Krishnankutty 2009,
J. Mar. Eng. Technol. 13, p. 14); the thrust of the momentum-flux form may
turn negative with boat speed when the intake flow exceeds the jet speed (ram
drag, the form's own sign).

* ``jet_bollard_quadratic`` (J1): ``T = k n^2`` for ``n > 0``. Anchor: the
  published bollard law ``F = 1.89e-4 N^2`` with ``N`` in rpm (Tajima et al.
  2015, eq. 7 and Fig. 3, p. 378), i.e. ``k = 1.89e-4 * 3600 N s^2/rev^2``:
  1000, 3000, 5000 rpm give ``k N^2``.
* ``jet_bilinear_speed`` (J2): ``T = b1 n^2 + b2 n Va``. Anchor: the thrust
  ratio of Han et al. (2017) eq. 7 with Table III: ``1 + (b11 k2 / b11 k1)
  u`` at full throttle.
* ``jet_momentum_flux`` (J3): ``T = rho Q (Q/A_n - Va)`` with ``Q = c_Q n``
  (Bulten 2006, eq. 2.52, p. 33, with ``v_out = Q/A_n``, ``v_in = Va``),
  ``c_Q`` given or from a bollard point (``T_b = rho c_Q^2 n_b^2 / A_n``, the
  inverse of Kandasamy et al. 2010 eq. 3 as the memo quotes it).
  - With ``Q`` proportional to ``n`` the law is *exactly* bilinear: J3 equals J2 with
    ``b1 = rho c_Q^2 / A_n``, ``b2 = -rho c_Q`` and equals J1 at ``Va = 0``
    with ``k = rho c_Q^2 / A_n``.
  - Anchor: the ideal bollard identity ``T = (4 rho A_n P^2)^(1/3)`` for a jet
    taking shaft power ``P`` (Loghis & Xiros 2022, App. A, eqs. A2-A8, pp.
    23-24), at the published pump point of Li & Kao (2023): nozzle diameter,
    5950 rpm, rotor torque (Tables 1, 3, 6).
* ``jet_pump_system_match`` (J4): ``Q`` from the intersection of a linear pump
  head line ``psi = psi_0 - s phi`` (``phi = Q / (Omega D^3)``, ``psi = g H /
  (Omega D)^2``; Bulten 2006, eqs. 2.15-2.16, p. 22) and the system head
  ``H_R = v_out^2 (1 + phi_n) / 2g - v_in^2 (1 - eps) / 2g`` (eq. 2.56, p. 35,
  ``h_j`` = 0); ``T = rho Q (v_out - Va)``; shaft torque ``rho g H Q / (eta_p
  Omega)`` (eq. 2.20, p. 22). The oracle solves the intersection by bisection
  (not the closed form of the part).

Author:    Enio Krizman
Date:      2026-10-10
"""

import math

import numpy as np
import pytest

from propulsor_contract import JET_LAWS, N_RANDOM, call, max_diff, part, published_value, rng

G = 9.81   # only inside the oracle's pump head; the part's own head ratio is dimensionless


def _law(f, n, va, rho=1025.0):
    out = call(f, drive_output=[n], advance_speed=[va], water_density=[rho])
    return out["thrust"][0], out["torque"][0]


def test_every_waterjet_law_is_listed_in_the_contract():
    assert set(JET_LAWS) == {"jet_bollard_quadratic", "jet_bilinear_speed", "jet_momentum_flux_given",
                             "jet_momentum_flux_bollard_point", "jet_pump_system_match"}


@pytest.mark.parametrize("variant, numbers", [
    ("jet_bollard_quadratic", {"bollard_thrust_coefficient": 0.68}),
    ("jet_bilinear_speed", {"speed_squared_coefficient": 0.68, "speed_advance_coefficient": -0.01}),
    ("jet_momentum_flux_given", {"nozzle_area": 7e-3, "flow_coefficient": 2e-3}),
    ("jet_momentum_flux_bollard_point", {"nozzle_area": 7e-3, "bollard_thrust": 4000.0, "bollard_shaft_speed": 90.0}),
    ("jet_pump_system_match", {"impeller_diameter": 0.15, "nozzle_area": 7e-3, "pump_head_coefficient_at_zero_flow": 2.0,
                               "pump_head_slope": 1.0, "nozzle_loss_coefficient": 0.02, "inlet_loss_coefficient": 0.2,
                               "pump_efficiency": 0.8}),
])
def test_no_waterjet_law_pushes_astern_and_a_stopped_pump_gives_nothing(variant, numbers):
    f = part("law", variant, numbers)
    g = rng()
    for n in g.uniform(-80.0, 0.0, 30):
        t, q = _law(f, n, g.uniform(0.0, 20.0))
        assert t == 0.0 and q == 0.0, (variant, n)
    assert _law(f, 0.0, 5.0) == (0.0, 0.0)


# --------------------------------------------------------------------------
# J1
# --------------------------------------------------------------------------
def test_J1_the_published_bollard_coefficient_gives_k_n_squared():
    k = published_value("bollard_coefficient_rpm") * 3600.0       # N/rpm^2 -> N s^2/rev^2
    f = part("law", "jet_bollard_quadratic", {"bollard_thrust_coefficient": k})
    for rpm in (1000.0, 3000.0, 5000.0):
        t, _ = _law(f, rpm / 60.0, 0.0)
        assert abs(t - published_value("bollard_coefficient_rpm") * rpm**2) <= 1e-9 * max(1.0, t)
    t_fast, _ = _law(f, 50.0, 7.5)                                  # J1 has no speed effect
    assert abs(t_fast - _law(f, 50.0, 0.0)[0]) == 0.0


# --------------------------------------------------------------------------
# J2
# --------------------------------------------------------------------------
@pytest.mark.parametrize("u", [1.0, 2.0, 3.5])
def test_J2_reproduces_the_published_thrust_speed_ratio_at_full_throttle(u):
    b1, b2 = published_value("han_thrust_square_coefficient"), published_value("han_thrust_speed_coefficient")
    f = part("law", "jet_bilinear_speed", {"speed_squared_coefficient": b1, "speed_advance_coefficient": b2})
    ratio = _law(f, 1.0, u)[0] / _law(f, 1.0, 0.0)[0]
    assert abs(ratio - (1.0 + b2 / b1 * u)) <= 1e-14
    assert 0.9 < ratio < 1.0                                        # the identified boat loses a few per cent


# --------------------------------------------------------------------------
# J3
# --------------------------------------------------------------------------
A_N, C_Q = 7e-3, 2e-3


def test_J3_is_rho_q_times_jet_speed_minus_intake_speed():
    f = part("law", "jet_momentum_flux_given", {"nozzle_area": A_N, "flow_coefficient": C_Q})
    g = rng()
    for n, va, rho in zip(g.uniform(1.0, 120.0, N_RANDOM), g.uniform(0.0, 30.0, N_RANDOM), g.uniform(990, 1030, N_RANDOM)):
        q = C_Q * n
        t, torque = _law(f, n, va, rho)
        assert abs(t - rho * q * (q / A_N - va)) <= 1e-12 * max(1.0, abs(t)) and torque == 0.0


def test_J3_thrust_falls_linearly_with_intake_speed_and_turns_into_drag_above_the_jet_speed():
    f = part("law", "jet_momentum_flux_given", {"nozzle_area": A_N, "flow_coefficient": C_Q})
    n = 60.0
    v_out = C_Q * n / A_N
    ts = [_law(f, n, va)[0] for va in (0.0, 0.25 * v_out, 0.5 * v_out, v_out, 1.5 * v_out)]
    assert all(a > b for a, b in zip(ts, ts[1:]))
    assert abs(ts[3]) <= 1e-9 and ts[4] < 0.0                       # zero net thrust when the intake runs at jet speed
    slopes = np.diff(ts[:3]) / np.diff([0.0, 0.25 * v_out, 0.5 * v_out])
    assert abs(slopes[0] - slopes[1]) <= 1e-9 * abs(slopes[0]) and abs(slopes[0] + 1025.0 * C_Q * n) <= 1e-6


def test_J3_equals_J2_with_the_derived_coefficients_exactly_and_J1_at_zero_speed():
    rho = 1025.0
    j3 = part("law", "jet_momentum_flux_given", {"nozzle_area": A_N, "flow_coefficient": C_Q})
    j2 = part("law", "jet_bilinear_speed", {"speed_squared_coefficient": rho * C_Q**2 / A_N,
                                            "speed_advance_coefficient": -rho * C_Q})
    j1 = part("law", "jet_bollard_quadratic", {"bollard_thrust_coefficient": rho * C_Q**2 / A_N})
    g = rng()
    for n, va in zip(g.uniform(0.1, 120.0, N_RANDOM), g.uniform(0.0, 30.0, N_RANDOM)):
        assert abs(_law(j3, n, va)[0] - _law(j2, n, va)[0]) <= 1e-11 * max(1.0, abs(_law(j3, n, va)[0]))
        assert abs(_law(j3, n, 0.0)[0] - _law(j1, n, 0.0)[0]) <= 1e-11 * max(1.0, _law(j3, n, 0.0)[0])


def test_J3_from_a_bollard_point_returns_that_bollard_thrust_and_scales_with_the_shaft_speed_squared():
    numbers = {"nozzle_area": A_N, "bollard_thrust": 4000.0, "bollard_shaft_speed": 90.0}
    f = part("law", "jet_momentum_flux_bollard_point", numbers)
    for rho in (1000.0, 1025.0, 1030.0):
        assert abs(_law(f, 90.0, 0.0, rho)[0] - 4000.0) <= 1e-9
    assert abs(_law(f, 45.0, 0.0)[0] - 1000.0) <= 1e-9
    assert _law(f, 90.0, 10.0)[0] < 4000.0
    # the same flow coefficient as the 'given' form: c_Q = sqrt(T_b A_n / rho) / n_b
    rho = 1025.0
    given = part("law", "jet_momentum_flux_given", {"nozzle_area": A_N, "flow_coefficient": math.sqrt(4000.0 * A_N / rho) / 90.0})
    for n, va in ((30.0, 3.0), (90.0, 10.0), (100.0, 0.0)):
        assert abs(_law(f, n, va, rho)[0] - _law(given, n, va, rho)[0]) <= 1e-9


def test_J3_ideal_bollard_identity_at_the_published_jet_ski_pump_point():
    """Loghis & Xiros 2022 eqs. A2-A8: a lossless jet absorbing P gives T = (4 rho A P^2)^(1/3). The law, with the flow
    that makes its jet power rho Q^3 / (2 A^2) equal the pump's published shaft power, gives the same thrust."""
    rho = 1025.0
    area = math.pi / 4.0 * published_value("jet_ski_pump_nozzle_outlet_diameter") ** 2
    power = 2.0 * math.pi * published_value("jet_ski_pump_shaft_speed_rpm") / 60.0 * published_value("jet_ski_pump_rotor_torque")
    q = (2.0 * area**2 * power / rho) ** (1.0 / 3.0)               # jet power rho Q^3 / (2 A^2) = P
    f = part("law", "jet_momentum_flux_given", {"nozzle_area": area, "flow_coefficient": q})   # n = 1 rev/s
    expected = (4.0 * rho * area * power**2) ** (1.0 / 3.0)
    assert abs(_law(f, 1.0, 0.0, rho)[0] - expected) <= 1e-9 * expected
    assert 1.0e4 < expected < 1.4e4                                 # the order the pump's kW and nozzle give


# --------------------------------------------------------------------------
# J4
# --------------------------------------------------------------------------
J4 = {"impeller_diameter": 0.15, "nozzle_area": 7e-3, "pump_head_coefficient_at_zero_flow": 2.0, "pump_head_slope": 1.0,
      "nozzle_loss_coefficient": 0.02, "inlet_loss_coefficient": 0.2, "pump_efficiency": 0.8}


def _j4_oracle(numbers, n, va, rho):
    """Intersection of the pump line and the system head by bisection (Bulten 2006 eqs. 2.15, 2.16, 2.20, 2.56)."""
    d, area, psi0, s = numbers["impeller_diameter"], numbers["nozzle_area"], numbers["pump_head_coefficient_at_zero_flow"], numbers["pump_head_slope"]
    omega = 2.0 * math.pi * n

    def head_pump(q):                                               # g H of the pump line
        phi = q / (omega * d**3)
        return (psi0 - s * phi) * (omega * d) ** 2

    def head_system(q):                                             # g H_R, h_j = 0
        return (q / area) ** 2 * (1.0 + numbers["nozzle_loss_coefficient"]) / 2.0 - va**2 * (1.0 - numbers["inlet_loss_coefficient"]) / 2.0

    low, high = 0.0, 5.0                                            # bisection: pump head falls, system head rises with Q
    for _ in range(200):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if head_pump(middle) - head_system(middle) > 0.0 else (low, middle)
    q = 0.5 * (low + high)
    thrust = rho * q * (q / area - va)
    torque = rho * head_pump(q) * q / (numbers["pump_efficiency"] * omega)
    return thrust, torque, q


def test_J4_equals_the_root_found_intersection_of_pump_line_and_system_head():
    f = part("law", "jet_pump_system_match", J4)
    g = rng()
    for n, va, rho in zip(g.uniform(10.0, 100.0, 200), g.uniform(0.0, 25.0, 200), g.uniform(990, 1030, 200)):
        t, torque = _law(f, n, va, rho)
        et, eq, _ = _j4_oracle(J4, n, va, rho)
        assert abs(t - et) <= 1e-8 * max(1.0, abs(et)) and abs(torque - eq) <= 1e-8 * max(1.0, abs(eq)), (n, va)


def test_J4_flat_pump_line_at_rest_is_a_bollard_square_law_with_the_derived_coefficient():
    numbers = {**J4, "pump_head_slope": 0.0, "nozzle_loss_coefficient": 0.0}
    f = part("law", "jet_pump_system_match", numbers)
    rho = 1025.0
    # Q^2 / (2 A^2) = psi0 (Omega D)^2  =>  T = rho Q^2 / A = 2 rho A psi0 (Omega D)^2
    k = 2.0 * rho * numbers["nozzle_area"] * numbers["pump_head_coefficient_at_zero_flow"] * (2.0 * math.pi * numbers["impeller_diameter"]) ** 2
    for n in (20.0, 60.0, 100.0):
        assert abs(_law(f, n, 0.0, rho)[0] - k * n**2) <= 1e-9 * k * n**2


def test_J4_thrust_falls_with_the_boat_speed_and_its_flow_rises_with_it_through_the_ram_term():
    f = part("law", "jet_pump_system_match", {**J4, "pump_head_slope": 0.0})
    thrusts = [_law(f, 80.0, va)[0] for va in (0.0, 5.0, 10.0, 15.0, 20.0, 25.0)]
    assert all(a > b for a, b in zip(thrusts, thrusts[1:])), thrusts
    flows = [_j4_oracle({**J4, "pump_head_slope": 0.0}, 80.0, va, 1025.0)[2] for va in (0.0, 10.0, 20.0)]
    assert flows[0] < flows[1] < flows[2]                            # the ram pressure helps the pump (Bulten p. 35)


def test_J4_loads_the_drive_the_torque_rises_with_flow_and_is_zero_for_a_stopped_pump():
    f = part("law", "jet_pump_system_match", J4)
    torques = [_law(f, n, 5.0)[1] for n in (20.0, 40.0, 80.0)]
    assert all(t > 0 for t in torques) and torques[0] < torques[1] < torques[2]
