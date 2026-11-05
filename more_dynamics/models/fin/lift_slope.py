"""Lift slope of a fin from its geometry (CasADi), the ``low_aspect_ratio`` form of ``lift_slope_form``.

``a = dC_L/d(alpha)`` per radian of a fin of effective aspect ratio ``AR_e``:

    AR_e = k_AR b^2 / S                               (Prestero 2001, eq. 4.39, p. 31; k_AR = 2: the body of the vehicle
                                                       acts as a reflection plane, so the wing is its own mirror image)
    a    = [ 1 / (2 pi a~) + 1 / (pi AR_e) ]^-1       (Hoerner's form, Prestero 2001, eq. 4.38, p. 31; a~ is the empirical
                                                       section factor, of the order 0.9)

``b`` is the span (m), ``S`` the planform area of the one fin (m^2); ``a`` and ``AR_e`` are dimensionless. The limits
check the form: for ``AR_e -> infinity`` the slope tends to ``2 pi a~`` (the section slope); for ``a~ -> infinity`` it tends
to ``pi AR_e`` (a slender fin).

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom simulation model for the REMUS autonomous
    underwater vehicle. MIT/WHOI MSc thesis. Eqs. 4.38-4.39, p. 31 (read as text).

Author:    Enio Krizman
"""
import math

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.shared.wiring import function_from

LOW_ASPECT_RATIO_LIFT_SLOPE_PARAMETERS = (
    Parameter("fin_span", (1, 1), "m", "span b of one fin", 0.0, minimum_exclusive=True),
    Parameter("fin_area", (1, 1), "m^2", "planform area S of one fin", 0.0, minimum_exclusive=True),
    Parameter("aspect_ratio_factor", (1, 1), "1", "k_AR of AR_e = k_AR b^2 / S; 2 for a fin on a body that acts as a "
              "reflection plane", 0.0, minimum_exclusive=True),
    Parameter("section_efficiency", (1, 1), "1", "the empirical factor a~ of Hoerner's low aspect ratio form, about 0.9",
              0.0, minimum_exclusive=True),
)


def low_aspect_ratio_lift_slope():
    """``(fin_span, fin_area, aspect_ratio_factor, section_efficiency) -> lift_slope``, 1/rad (module docstring)."""
    span, area = ca.SX.sym("fin_span"), ca.SX.sym("fin_area")
    factor, efficiency = ca.SX.sym("aspect_ratio_factor"), ca.SX.sym("section_efficiency")
    aspect_ratio = factor * span**2 / area  # AR_e (Prestero 2001, eq. 4.39, p. 31)
    slope = 1.0 / (1.0 / (2.0 * math.pi * efficiency) + 1.0 / (math.pi * aspect_ratio))  # (Prestero 2001, eq. 4.38, p. 31)
    return function_from("low_aspect_ratio_lift_slope", {"fin_span": span, "fin_area": area, "aspect_ratio_factor": factor,
                                                         "section_efficiency": efficiency}, {"lift_slope": slope})
