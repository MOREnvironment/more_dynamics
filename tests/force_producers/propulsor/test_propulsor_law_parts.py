"""The open-water, quadratic and normalised law parts of the propulsor: drive output and advance speed -> thrust, torque.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

K_T(J) and K_Q(J) come from the frozen MATLAB rows of MSS ``wageningen.m``
(``tests/force_producers/data/wageningen_kt_kq_mss_current.csv``): with
density, diameter and shaft speed all 1 the thrust of a law is K_T(J) and the
torque K_Q(J), J = advance speed.

* ``open_water_signed_advance``: signed ``J = Va / (n D)``, the polynomial at
  that J (no clip by default: MSS), reverse factors on K_T, K_Q for n < 0,
  zero force when stopped (the steerable outboard block's construction).
* ``open_water_polynomial``: ``J = Va / (|n| D)`` saturated to ``[0, 2
  Ja_max]``; K_T(0), K_Q(0) for ``n <= 0``.
* ``open_water_linearized``: MSS ``remus100.m`` 167-176 with K_T0, K_Q0, K_Tmax,
  K_Qmax from the polynomial at J = 0 and Ja_max (or the printed four given).
* ``open_water_bollard``: ``remus100.m`` 175-176, K_T0, K_Q0.
* ``quadratic_asymmetric``: ``T = k w |w|``, ``w = 2 pi n``, ``k_pos`` ahead,
  ``k_neg`` astern (otter.m 223-229); torque 0.
* ``normalized_linear``: ``T = T_max c effort``, torque 0, indifferent to the
  advance speed and the density (unit "1" input).

Equations: T = rho D^4 K_T n|n|, Q = rho D^5 K_Q n|n| (Fossen 2011, eq.
12.265, p. 411).

Author:    Enio Krizman
Date:      2026-10-10
"""

import numpy as np
import pytest

from propulsor_contract import (FORCE_DATA_DIR, G1_TOLERANCE, G2_TOLERANCE, N_RANDOM, _csv_columns, call, max_diff, part,
                                remus_propeller_values, rng)

GEOMETRIES = [(1.0, 0.718, 3.0), (0.83, 0.55, 4.0), (0.83, 0.718, 3.0)]


def _rows():
    return _csv_columns(FORCE_DATA_DIR / "wageningen_kt_kq_mss_current.csv")


def _geometry_rows(geometry):
    ref = _rows()
    mask = (ref["PD"] == geometry[0]) & (ref["AEAO"] == geometry[1]) & (ref["z"] == geometry[2])
    assert mask.sum() > 1000, geometry
    return {k: v[mask] for k, v in ref.items()}


def _geometry_values(geometry, **extra):
    return {"diameter": 1.0, "pitch_diameter_ratio": geometry[0], "blade_area_ratio": geometry[1],
            "blade_count": geometry[2], **extra}


def _law(f, n, va, rho=1.0):
    out = call(f, drive_output=[n], advance_speed=[va], water_density=[rho])
    return out["thrust"][0], out["torque"][0]


# --------------------------------------------------------------------------
# G1: K_T, K_Q of MSS wageningen.m through the laws
# --------------------------------------------------------------------------
@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_G1_signed_advance_law_gives_the_matlab_polynomial_at_every_row(geometry):
    ref = _geometry_rows(geometry)
    f = part("law", "open_water_signed_advance",
             _geometry_values(geometry, reverse_thrust_factor=1.0, reverse_torque_factor=1.0))
    worst = 0.0
    for j, kt, kq in zip(ref["J"], ref["KT"], ref["KQ"]):
        t, q = _law(f, 1.0, j)
        worst = max(worst, abs(t - kt), abs(q - kq))
    assert worst <= G1_TOLERANCE, worst


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_G1_signed_advance_law_astern_applies_the_reverse_factors(geometry):
    ref = _geometry_rows(geometry)
    f = part("law", "open_water_signed_advance",
             _geometry_values(geometry, reverse_thrust_factor=0.5, reverse_torque_factor=0.6))
    worst = 0.0
    for j, kt, kq in list(zip(ref["J"], ref["KT"], ref["KQ"]))[::7]:
        t, q = _law(f, -1.0, -j)                              # n = -1: J = Va / (n D) = j
        worst = max(worst, abs(t + 0.5 * kt), abs(q + 0.6 * kq))
    assert worst <= G1_TOLERANCE, worst


def test_the_signed_advance_law_gives_no_force_when_stopped():
    f = part("law", "open_water_signed_advance", _geometry_values(GEOMETRIES[2], reverse_thrust_factor=0.5,
                                                                  reverse_torque_factor=0.6))
    assert _law(f, 0.0, 3.0, 1025.0) == (0.0, 0.0)


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_G1_polynomial_law_gives_the_matlab_polynomial_on_the_forward_branch(geometry):
    ref = _geometry_rows(geometry)
    f = part("law", "open_water_polynomial", _geometry_values(geometry, max_advance_number=1.25))   # clip at J = 2.5
    j0 = float(ref["J"][np.argmin(np.abs(ref["J"]))])
    kt0, kq0 = float(ref["KT"][ref["J"] == j0][0]), float(ref["KQ"][ref["J"] == j0][0])
    assert j0 == 0.0
    worst = 0.0
    for j, kt, kq in zip(ref["J"], ref["KT"], ref["KQ"]):
        t, q = _law(f, 1.0, j)
        if j < 0.0:                                             # J is saturated to [0, 2 Ja_max]
            worst = max(worst, abs(t - kt0), abs(q - kq0))
        else:
            worst = max(worst, abs(t - kt), abs(q - kq))
    assert worst <= G1_TOLERANCE, worst


def test_the_polynomial_law_uses_the_zero_advance_coefficients_astern_and_saturates_the_advance_number():
    geometry = GEOMETRIES[0]
    ref = _geometry_rows(geometry)
    kt0 = float(ref["KT"][ref["J"] == 0.0][0])
    kq0 = float(ref["KQ"][ref["J"] == 0.0][0])
    f = part("law", "open_water_polynomial", _geometry_values(geometry, max_advance_number=0.6632))
    t, q = _law(f, -2.0, 1.5, 1.0)
    assert abs(t - kt0 * (-4.0)) <= G1_TOLERANCE and abs(q - kq0 * (-4.0)) <= G1_TOLERANCE      # K(0) n|n|, n = -2
    k_at_2jmax = float(ref["KT"][np.argmin(np.abs(ref["J"] - 2 * 0.6632))])                     # nearest row, near 1.33
    t_big, _ = _law(f, 1.0, 50.0)                                                                # J saturated at 2 Ja_max
    assert abs(t_big - k_at_2jmax) < 0.02


# --------------------------------------------------------------------------
# G2: the three REMUS-style laws against the stated equations
# --------------------------------------------------------------------------
def _remus_k():
    ref = _geometry_rows(GEOMETRIES[0])
    pick = lambda j, name: float(ref[name][np.argmin(np.abs(ref["J"] - j))])  # noqa: E731
    jmax = remus_propeller_values()["max_advance_number"]
    assert abs(ref["J"][np.argmin(np.abs(ref["J"] - jmax))] - jmax) < 1e-12      # the rows hold J = Ja_max itself
    return pick(0.0, "KT"), pick(0.0, "KQ"), pick(jmax, "KT"), pick(jmax, "KQ")


def _linearized_oracle(k, rho, d, jmax, n, va):
    kt0, kq0, ktm, kqm = k
    quad = abs(n) * n
    if n > 0:
        return (rho * d**4 * (kt0 * quad + (ktm - kt0) / jmax * (va / d) * abs(n)),
                rho * d**5 * (kq0 * quad + (kqm - kq0) / jmax * (va / d) * abs(n)))
    return rho * d**4 * kt0 * quad, rho * d**5 * kq0 * quad


def _remus_law_values(**extra):
    v = remus_propeller_values()
    return {"diameter": v["propeller_diameter"], "pitch_diameter_ratio": v["pitch_diameter_ratio"],
            "blade_area_ratio": v["blade_area_ratio"], "blade_count": v["blade_count"],
            "max_advance_number": v["max_advance_number"], **extra}


@pytest.mark.parametrize("variant, coefficients", [("open_water_linearized", None),
                                                   ("open_water_linearized_given", "given")])
def test_G2_linearized_law_is_remus100_m_167_176(variant, coefficients):
    v = remus_propeller_values()
    k = _remus_k() if coefficients is None else tuple(v["thrust_torque_coefficients"])
    extra = {"thrust_torque_coefficients": v["thrust_torque_coefficients"]} if coefficients else {}
    f = part("law", variant, _remus_law_values(**extra))
    g = rng()
    forward = astern = 0
    for n, va, rho in zip(g.uniform(-25.0, 25.0, N_RANDOM), g.uniform(0.0, 3.0, N_RANDOM), g.uniform(1000, 1030, N_RANDOM)):
        t, q = _law(f, n, va, rho)
        et, eq = _linearized_oracle(k, rho, v["propeller_diameter"], v["max_advance_number"], n, va)
        assert abs(t - et) <= G2_TOLERANCE * max(1.0, abs(et)) and abs(q - eq) <= G2_TOLERANCE * max(1.0, abs(eq)), (n, va)
        forward += n > 0
        astern += n < 0
    assert forward > 50 and astern > 50


def test_G2_bollard_law_is_remus100_m_175_176():
    v = remus_propeller_values()
    kt0, kq0, _, _ = _remus_k()
    f = part("law", "open_water_bollard", _remus_law_values())
    g = rng()
    for n, va in zip(g.uniform(-25.0, 25.0, N_RANDOM), g.uniform(0.0, 3.0, N_RANDOM)):
        t, q = _law(f, n, va, 1026.0)
        assert abs(t - 1026.0 * v["propeller_diameter"] ** 4 * kt0 * abs(n) * n) <= G2_TOLERANCE * max(1.0, abs(t))
        assert abs(q - 1026.0 * v["propeller_diameter"] ** 5 * kq0 * abs(n) * n) <= G2_TOLERANCE * max(1.0, abs(q))


def test_the_linearized_law_is_the_bollard_law_at_zero_advance_speed_and_astern():
    values = _remus_law_values()
    lin, boll = part("law", "open_water_linearized", values), part("law", "open_water_bollard", values)
    g = rng()
    for n in g.uniform(-25.0, 25.0, 50):
        assert max_diff(_law(lin, n, 0.0, 1026.0), _law(boll, n, 0.0, 1026.0)) <= 1e-12
    for n, va in zip(g.uniform(-25.0, -1.0, 50), g.uniform(0.0, 3.0, 50)):
        assert max_diff(_law(lin, n, va, 1026.0), _law(boll, n, va, 1026.0)) <= 1e-12


def test_the_linearized_law_gains_thrust_with_advance_speed_less_than_bollard_says():
    """Physical direction: K_T falls with J, so the linearised thrust at a given n falls as the advance speed rises."""
    lin = part("law", "open_water_linearized", _remus_law_values())
    thrusts = [_law(lin, 20.0, va, 1026.0)[0] for va in (0.0, 0.5, 1.0, 2.0)]
    assert all(a > b for a, b in zip(thrusts, thrusts[1:])), thrusts


# --------------------------------------------------------------------------
# quadratic_asymmetric, normalized_linear
# --------------------------------------------------------------------------
def test_G2_quadratic_asymmetric_is_the_otter_thrust_law():
    k_pos, k_neg = 0.01108, 0.006445        # test numbers; only the structure is gated here (the frozen rows are in the composition tests)
    f = part("law", "quadratic_asymmetric", {"positive_thrust_coefficient": k_pos, "negative_thrust_coefficient": k_neg})
    g = rng()
    for n, va, rho in zip(g.uniform(-20.0, 20.0, N_RANDOM), g.uniform(-3, 3, N_RANDOM), g.uniform(900, 1100, N_RANDOM)):
        w = 2.0 * np.pi * n
        expected = (k_pos if w > 0 else k_neg) * w * abs(w)
        t, q = _law(f, n, va, rho)
        assert abs(t - expected) <= 1e-12 * max(1.0, abs(expected)) and q == 0.0
    assert _law(f, 0.0, 1.0)[0] == 0.0


def test_G2_normalized_linear_is_the_jet_nozzle_thrust_law_and_indifferent_to_flow():
    f = part("law", "normalized_linear", {"max_thrust": 800.0, "thrust_coefficient": 0.9})
    g = rng()
    for effort, va, rho in zip(g.uniform(-1.0, 1.0, 50), g.uniform(-5, 20, 50), g.uniform(900, 1100, 50)):
        t, q = _law(f, effort, va, rho)
        assert abs(t - 800.0 * 0.9 * effort) <= 1e-12 and q == 0.0
