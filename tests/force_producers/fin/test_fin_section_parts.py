"""Gate tests of the fin section parts (angle of attack -> lift and drag
coefficients): ``quadratic_drag``, ``linear_section`` now, ``lifting_line``
written and waiting.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).
Coupling input ``angle_of_attack`` (1, rad); outputs ``lift_coefficient``
and ``drag_coefficient`` (1 each). The skeleton turns them into the force
``f = 1/2 rho A U^2 (C_L lift_axis - C_D chord_axis)``.

* ``quadratic_drag``: ``C_L = lift_slope * alpha``, ``C_D = lift_slope *
  alpha^2`` (MSS ``remus100.m`` 238-245: ``X_r = -1/2 rho U_rh^2 A_r
  CL_delta_r delta_r^2``, ``Y_r = -1/2 rho U_rh^2 A_r CL_delta_r delta_r``;
  stern planes alike). G1 against MATLAB: the coefficients MATLAB's forces
  imply, ``C_L = -Y_r / q``, ``C_D = -X_r / q``, ``q = 1/2 rho U_rh^2 A_r``
  (rudder lifting along -y), and ``C_L = -Z_s / q``, ``C_D = -X_s / q`` with
  ``U_rv``, ``A_s`` (stern planes lifting along -z).
* ``linear_section``: ``C_L = lift_slope * alpha``, ``C_D =
  zero_lift_drag`` (Prestero 2001 eq. 4.37, p. 31, is lift only: with
  ``zero_lift_drag = 0`` the section is his).
* ``lifting_line`` (waiting): ``lift_slope`` from the effective aspect ratio,
  ``c_L_alpha = [1 / (2 a_tilde pi) + 1 / (pi AR_e)]^-1`` (Prestero 2001,
  eq. 4.38, p. 31, after Hoerner, **Hoerner unread**), ``AR_e = 2 b^2 / S``
  (eq. 4.39); its induced drag (lifting-line references, unread) fixes the
  rest of its contract. The two equations of Prestero's page are checked
  against his own Table A.5 now, without the part.

The lift slope is positive in every form: the sign of a fin's lift is its
``lift_axis`` (geometry of the skeleton), not the sign of a coefficient.

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 179-189,
    238-245.
[Prestero 2001] Prestero, T. (2001). *Verification of a six-degree of
    freedom simulation model for the REMUS autonomous underwater vehicle*.
    MSc thesis, MIT/WHOI. Ch. 4, eqs. 4.37-4.39, p. 31; App. A, Table A.5,
    p. 103.

Author:    Enio Krizman
Date:      2026-10-08
"""

import numpy as np
import pytest

from fin_parts_contract import (
    DECLARED,
    G1_TOLERANCE,
    G2_TOLERANCE,
    G4_FACTOR,
    N_RANDOM,
    WAITING_REASON,
    call,
    mss_constant,
    mss_fins,
    part,
    part_declared,
    prestero,
    prestero_value,
    rng,
)

FORMS = ("quadratic_drag", "linear_section")


@pytest.mark.parametrize("form", FORMS)
def test_section_declares_its_primitives(form):
    declared = part_declared("section", form)
    expected = DECLARED[("section", form)]
    assert [d.name for d in declared] == list(expected)
    for d in declared:
        shape, unit, minimum, exclusive = expected[d.name]
        assert (d.shape, d.unit, d.minimum, d.minimum_exclusive) == (shape, unit, minimum, exclusive), d
        assert d.meaning
        assert "water_density" not in d.name


@pytest.mark.parametrize("form", FORMS)
def test_section_refuses_values_outside_their_ranges(form):
    from more_transformations.more_casadi_transformations import check_values

    declared = part_declared("section", form)
    good = {d.name: 0.5 for d in declared}
    check_values(declared, good)
    with pytest.raises(ValueError, match="lift_slope"):
        check_values(declared, {**good, "lift_slope": 0.0})
    with pytest.raises(ValueError, match="lift_slope"):
        check_values(declared, {**good, "lift_slope": -3.0})
    if "zero_lift_drag" in good:
        check_values(declared, {**good, "zero_lift_drag": 0.0})
        with pytest.raises(ValueError, match="zero_lift_drag"):
            check_values(declared, {**good, "zero_lift_drag": -0.01})


# --------------------------------------------------------------------------
# quadratic_drag
# --------------------------------------------------------------------------
def _mss_coefficients(ref, fin):
    """(alpha, C_L, C_D, lift_slope) implied by MATLAB's forces, rows with a
    section-plane speed squared of at least 0.01 m^2/s^2."""
    rho = mss_constant(ref, "rho")
    if fin == "rudder":
        speed, area, slope = ref["U_rh"], mss_constant(ref, "A_r"), mss_constant(ref, "CL_delta_r")
        alpha, drag, lift = ref["delta_r"], ref["X_r"], ref["Y_r"]
    else:
        speed, area, slope = ref["U_rv"], mss_constant(ref, "A_s"), mss_constant(ref, "CL_delta_s")
        alpha, drag, lift = ref["delta_s"], ref["X_s"], ref["Z_s"]
    keep = speed ** 2 >= 0.01
    q = 0.5 * rho * speed[keep] ** 2 * area
    return alpha[keep], -lift[keep] / q, -drag[keep] / q, slope


def _worst_against_mss(fin, factor=1.0):
    ref = mss_fins()
    alpha, c_l, c_d, slope = _mss_coefficients(ref, fin)
    section = part("section", "quadratic_drag", {"lift_slope": factor * slope})
    worst = 0.0
    for a, l, d in zip(alpha, c_l, c_d):
        out = call(section, angle_of_attack=a)
        worst = max(worst, abs(out["lift_coefficient"][0] - l), abs(out["drag_coefficient"][0] - d))
    return worst, len(alpha)


@pytest.mark.parametrize("fin", ["rudder", "stern_plane"])
def test_G1_quadratic_drag_equals_the_coefficients_of_mss(fin):
    worst, rows = _worst_against_mss(fin)
    assert rows > 900
    assert worst <= G1_TOLERANCE


def test_G4_quadratic_drag_lift_slope_plus_1_percent_is_detected():
    assert _worst_against_mss("rudder", 1.01)[0] > G4_FACTOR * G1_TOLERANCE


def test_G2_quadratic_drag_equals_its_equation():
    section = part("section", "quadratic_drag", {"lift_slope": 0.7})
    for a in rng().uniform(-0.8, 0.8, N_RANDOM):
        out = call(section, angle_of_attack=a)
        assert abs(out["lift_coefficient"][0] - 0.7 * a) <= G2_TOLERANCE
        assert abs(out["drag_coefficient"][0] - 0.7 * a ** 2) <= G2_TOLERANCE


# --------------------------------------------------------------------------
# linear_section
# --------------------------------------------------------------------------
def test_G2_linear_section_equals_prestero_eq_4_37():
    """C_L = c_L_alpha alpha with Prestero's c_L_alpha (Table A.5, p. 103)
    and no fin drag (eq. 4.37, p. 31)."""
    slope = prestero_value("lift_slope")
    section = part("section", "linear_section", {"lift_slope": slope, "zero_lift_drag": 0.0})
    for a in rng().uniform(-0.8, 0.8, N_RANDOM):
        out = call(section, angle_of_attack=a)
        assert abs(out["lift_coefficient"][0] - slope * a) <= G2_TOLERANCE
        assert out["drag_coefficient"][0] == 0.0


def test_G2_linear_section_zero_lift_drag_is_constant():
    section = part("section", "linear_section", {"lift_slope": 2.0, "zero_lift_drag": 0.012})
    for a in rng().uniform(-0.8, 0.8, 50):
        assert abs(call(section, angle_of_attack=a)["drag_coefficient"][0] - 0.012) <= G2_TOLERANCE


# --------------------------------------------------------------------------
# Prestero's lift slope from the aspect ratio (source check, no part needed)
# --------------------------------------------------------------------------
def test_G5_prestero_effective_aspect_ratio_eq_4_39_reproduces_table_A5():
    span, area = prestero_value("fin_span"), prestero_value("fin_area")
    half = prestero()["printed_half_unit"]
    effective = 2.0 * span ** 2 / area                          # (Prestero 2001, eq. 4.39, p. 31)
    assert abs(effective - prestero_value("effective_aspect_ratio")) <= half["Table A.5 AR_e"]


def test_G5_prestero_lift_slope_eq_4_38_reproduces_table_A5():
    a_tilde, aspect = prestero_value("lift_slope_parameter"), prestero_value("effective_aspect_ratio")
    half = prestero()["printed_half_unit"]
    slope = 1.0 / (1.0 / (2.0 * a_tilde * np.pi) + 1.0 / (np.pi * aspect))   # (eq. 4.38, p. 31)
    assert abs(slope - prestero_value("lift_slope")) <= half["Table A.5 c_L_alpha"]


# --------------------------------------------------------------------------
# lifting_line (waiting)
# --------------------------------------------------------------------------
@pytest.mark.skip(reason=WAITING_REASON)
def test_lifting_line_lift_slope_equals_prestero_eq_4_38():
    """lifting_line at Prestero's AR_e, a_tilde gives c_L_alpha 3.12 (Table
    A.5): its C_L equals linear_section's with that slope."""
    raise AssertionError("contract not fixed")


@pytest.mark.skip(reason=WAITING_REASON)
def test_lifting_line_in_the_slender_limit_is_a_linear_section():
    """AR_e -> infinity: induced drag -> 0 and the slope -> 2 pi a_tilde; the
    section then equals linear_section (the reduction arrow F2 -> F1)."""
    raise AssertionError("contract not fixed")
