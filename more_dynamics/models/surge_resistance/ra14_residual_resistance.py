"""The RA14 residual-resistance prior for a hard-chine hull (Radojcic,
Zgradic, Kalajdzic, Simic 2014, Appendix 1, p. 24): ``R/Delta`` and
``S/V^(2/3)`` as cubics in the volumetric Froude number ``Fn_V``, each
coefficient a polynomial in the slenderness ``s = L / V^(1/3)``,
Froude-scaled from the paper's own 100 000 lb basis (p. 11: ITTC-1957
friction, ``C_A = 0``). Shared by the surge resistance's residual
coefficient (``models/surge_resistance``) and the hydrostatics' static
wetted-surface option (``models/restoring``), so the table is typed once.

``residual_coefficient(displaced_volume, length, speed)`` gives ``C_R(Fn_V)``:
scale the hull (displaced volume, length, speed) to RA14's standard craft
(same ``s``, same ``Fn_V``), read ``R/Delta`` and ``S/V^(2/3)`` there, form
``C_T,std = R_std / (1/2 rho_std S_std U_std^2)``, and subtract the standard
craft's own ITTC friction ``C_F,std`` (RA14 p. 11 basis) -- the residual that
transfers. The scaling procedure is a construction of A-64
(`agents-more/20_sources/monohull_structure_for_grethe.md` hidden assumption
3), not an equation printed in the paper. Below ``Fn_V = 0.6`` the regression
does not cover the hull and ``C_R`` is clamped to 0 (A-64, same place).

``wetted_surface_coefficient(displaced_volume, length)`` reads the same
``S/V^(2/3)`` regression at ``Fn_V = 0`` (the cubic's constant term, ``D(s)``
only): a hydrostatic quantity has no speed of its own.

Coefficients and the scaling procedure are carried over with citation from
``agents-more/60_working/A64_scripts/a64_hull_regime_resistance.py`` (A-64,
transcribed and verified there from the cited page; rule 4: not retyped from
the paper here).

References
----------
[RA14] Radojcic, D., Zgradic, A., Kalajdzic, M., Simic, A. (2014). Resistance
    prediction for hard chine hulls in the pre-planing regime. Polish
    Maritime Research 21(2):9-26. Appendix 1, p. 24 (the simple model: R/Delta,
    S/V^(2/3) cubic in Fn_V, coefficients quartic/cubic in L/V^(1/3)); basis,
    p. 11 (Delta = 100000 lb, rho = 1026 kg/m^3, nu = 1.1907e-6 m^2/s,
    ITTC-1957 friction, C_A = 0); validity, p. 17.
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and
    Motion Control. Wiley. Eq. 6.83, p. 125 (the ITTC-1957 friction line, used
    for the standard craft's own friction).

Author:    Enio Krizman
Date:      2026-10-09
"""

import casadi as ca

# RA14 Appendix 1, p. 24: R/Delta, S/V^(2/3), each cubic in Fn_V with
# coefficients A, B, C, D quartic/cubic in s = L/V^(1/3). Copied with
# citation from agents-more/60_working/A64_scripts/a64_hull_regime_resistance.py
# (A-64, transcribed and verified there); not retyped from the paper here
# (rule 4).
R_COEFFICIENTS = {
    "A": (-0.0048694, 0.1057838, -0.8432151, 2.8994541, -3.5683179),
    "B": (0.0301221, -0.6562651, 5.2443776, -18.0560600, 22.1656778),
    "C": (-0.0508810, 1.1119496, -8.9006694, 30.6066779, -37.2112557),
    "D": (0.0209140, -0.4573261, 3.6580101, -12.5431467, 15.14353),
}
S_COEFFICIENTS = {
    "A": (0.0197989, -0.2876721, 1.3944044, -2.1175915),
    "B": (-0.1612880, 2.3295698, -11.3511782, 18.0264798),
    "C": (0.4226973, -6.1030996, 29.8350864, -49.7163001),
    "D": (-0.4432887, 6.7794188, -33.1423777, 58.8177152),
}

# ITTC-1957 friction line (Fossen 2011, eq. 6.83, p. 125), the form the
# standard craft's own friction uses (RA14 p. 11 basis: ITTC-1957, C_A = 0).
ITTC_FRICTION_FACTOR = 0.075
ITTC_LOG10_OFFSET = 2.0

GRAVITY = 9.81  # RA14 p. 11 basis
STANDARD_WATER_DENSITY = 1026.0  # RA14 p. 11 basis: rho_s
STANDARD_KINEMATIC_VISCOSITY = 1.1907e-6  # RA14 p. 11 basis: nu_s
STANDARD_MASS = 100000.0 * 0.45359237  # RA14 p. 11 basis: Delta = 100000 lb

# A numerical floor of this one scaling step only, not RA14's own (A-64's
# construction, avoids log10(0) at the standard craft's Reynolds number).
_STANDARD_REYNOLDS_FLOOR = 1.0
# Numerical device: avoids division by zero in the Froude scaling at rest.
# The surge-damping blend multiplies this term by zero there (its own
# sigma = 1 - tanh(|u_r| / u_cross) = 1 at u_r = 0), so the result is
# unaffected; the same role as surge_damping.py's Reynolds-number log offset
# (``REYNOLDS_LOG_OFFSET``).
_SPEED_FLOOR = 1e-9

VALIDITY_FROUDE_NUMBER = 0.6  # below it the regression does not cover the hull (A-64 hidden assumption 3)


def _poly(coefficients, x):
    value = ca.SX(0.0)
    n = len(coefficients)
    for i, c in enumerate(coefficients):
        value = value + c * x ** (n - 1 - i)
    return value


def _cubic_in_froude(coefficients, s, froude):
    a, b, c, d = (_poly(coefficients[k], s) for k in "ABCD")
    return a * froude**3 + b * froude**2 + c * froude + d


def residual_coefficient(displaced_volume, length, speed):
    """``C_R(Fn_V)``: the RA14 residual coefficient (module docstring),
    0 below ``VALIDITY_FROUDE_NUMBER``."""
    s = length / displaced_volume ** (1.0 / 3.0)
    speed = ca.fmax(ca.fabs(speed), _SPEED_FLOOR)
    froude = speed / ca.sqrt(GRAVITY * displaced_volume ** (1.0 / 3.0))
    r_over_delta = _cubic_in_froude(R_COEFFICIENTS, s, froude)
    s_coefficient = _cubic_in_froude(S_COEFFICIENTS, s, froude)
    volume_standard = STANDARD_MASS / STANDARD_WATER_DENSITY
    scale = (volume_standard / displaced_volume) ** (1.0 / 3.0)
    speed_standard = speed * ca.sqrt(scale)
    length_standard = length * scale
    area_standard = s_coefficient * volume_standard ** (2.0 / 3.0)
    resistance_standard = r_over_delta * STANDARD_MASS * GRAVITY
    total_coefficient_standard = resistance_standard / (
        0.5 * STANDARD_WATER_DENSITY * area_standard * speed_standard**2)
    reynolds_standard = ca.fmax(speed_standard * length_standard / STANDARD_KINEMATIC_VISCOSITY,
                               _STANDARD_REYNOLDS_FLOOR)
    friction_standard = ITTC_FRICTION_FACTOR / (ca.log10(reynolds_standard) - ITTC_LOG10_OFFSET) ** 2
    residual = total_coefficient_standard - friction_standard
    return ca.if_else(froude < VALIDITY_FROUDE_NUMBER, 0.0, residual)


def wetted_surface_coefficient(displaced_volume, length):
    """``S / V^(2/3)`` at ``Fn_V = 0`` (the cubic's constant term, ``D(s)``
    only): the hydrostatic value, no speed of its own (module docstring)."""
    s = length / displaced_volume ** (1.0 / 3.0)
    return _poly(S_COEFFICIENTS["D"], s)
