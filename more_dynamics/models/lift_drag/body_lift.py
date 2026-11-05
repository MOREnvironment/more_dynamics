"""Body lift of a slender body of revolution as one point force (CasADi), its
slope from the slenderness and its pressure centre from the length.

Frame: BODY FRD about the origin CO (x forward, y starboard, z down); the hull is
centred on the origin with the nose at ``x = +l/2``. ``nu_r = [u, v, w, p, q, r]``
is the velocity relative to the water, SI units; the force and moment are N and
N m about the CO.

Equations (keys in References):

* The lift force of a body of revolution in cross-flow is
  ``L = -1/2 rho d^2 c_ydb u w`` in the vertical plane (Prestero 2001, eq. 4.33,
  p. 31, from eq. 4.27, p. 30), with ``c_ydb`` the lift slope per radian and ``d``
  the maximum diameter. Its coefficient is ``Z_uwl = -1/2 rho d^2 c_ydb`` and,
  by the symmetry of revolution, ``Y_uvl = Z_uwl`` in the horizontal plane
  (eq. 4.34, p. 31). Both planes: ``F = [0, -k u v, -k u w]`` with
  ``k = 1/2 rho d^2 c_ydb``.
* The force acts at one point on the axis, the pressure centre
  ``r_cp = (x_cp, 0, 0)``: ``tau = [F; r_cp x F]`` (Fossen 2011, eq. 12.226,
  p. 400). The moments are ``M = -x_cp F_z`` and ``N = x_cp F_y``.
  **Deviation from the printed eq. 4.36 (p. 31):** it prints
  ``M_uwl = -N_uvl = -1/2 rho d^2 c_ydb x_cp``, which is ``M = +x_cp F_z``. A
  force along z at ``x_cp`` has the moment ``-x_cp F_z`` (the cross product of
  the arm and the force), so the printed sign is a typesetting error of the
  thesis; the opposite sign is used here, as the ratio ``M_wl / Z_wl = -x_cp``
  of the linearised coefficients of the same thesis (Table D.1, p. 112).
* Slender-body lift slope: ``c_ydb = (l/d) c_yb 180/pi``, with ``c_yb`` the slope
  per degree (Prestero 2001, eqs. 4.30 and 4.32, p. 30; the value 0.003 of
  eq. 4.31 is Hoerner's, stated for ``6.7 <= l/d <= 10``). Outside that range of
  ``l/d`` the slope is not stated by the source: ``check_slenderness`` refuses
  it, naming the ratio (a hull outside the range takes the ``given`` slope).
* Pressure centre: at the fraction ``f`` of the length from the nose,
  ``x_cp = l/2 - f l`` (Prestero 2001, eq. 4.35, p. 31: ``f = 0.65``, Hoerner's
  range 0.6 to 0.7).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and
    Motion Control. Wiley. Eq. 12.226, p. 400.
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc
    thesis. Eqs. 4.27-4.36, pp. 30-31 (read as text, pp. 30-31); Table D.1, p. 112
    (not re-read).

Author:    Enio Krizman
"""

import math

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.shared.force_producer_common import wrench
from more_dynamics.models.shared.wiring import function_from

SLENDERNESS_RANGE = (6.7, 10.0)  # l/d for which the slope per degree is stated (Prestero 2001, eq. 4.31, p. 30)

SLENDER_BODY_LIFT_SLOPE_PARAMETERS = (
    Parameter("body_lift_slope_per_degree", (1, 1), "1/deg",
              "slope c_yb per degree of the body lift coefficient in units of the slenderness l/d", 0.0,
              minimum_exclusive=True),
)
FRACTION_PRESSURE_CENTRE_PARAMETERS = (
    Parameter("pressure_centre_fraction", (1, 1), "1",
              "distance of the pressure centre of the body lift from the nose as a fraction of the length", 0.0, 1.0,
              minimum_exclusive=True, maximum_exclusive=True),
)
GIVEN_LIFT_SLOPE_PARAMETERS = (
    Parameter("body_lift_slope", (1, 1), "1/rad", "lift slope c_ydb of the body, on the squared diameter", 0.0,
              minimum_exclusive=True),
)
GIVEN_PRESSURE_CENTRE_PARAMETERS = (
    Parameter("pressure_centre", (1, 1), "m", "x of the pressure centre of the body lift from the CO, BODY"),
)


def check_slenderness(length, diameter):
    """Refuse a hull outside ``6.7 <= l/d <= 10`` (Prestero 2001, eq. 4.31, p. 30), naming the ratio; returns it."""
    ratio = float(length) / float(diameter)
    low, high = SLENDERNESS_RANGE
    if not low <= ratio <= high:
        raise ValueError(f"slender_body_hoerner is stated for {low:g} <= l/d <= {high:g}; this hull has l/d = {ratio:.4g}: "
                         "give the slope with body_lift_slope_method 'given'")
    return ratio


def slender_body_lift_slope():
    """``(length, span, body_lift_slope_per_degree) -> body_lift_slope``: ``c_ydb = (l/d) c_yb 180/pi``, 1/rad
    (Prestero 2001, eqs. 4.30, 4.32, p. 30); ``span`` is the diameter of the body of revolution."""
    length, span = ca.SX.sym("length"), ca.SX.sym("span")
    per_degree = ca.SX.sym("body_lift_slope_per_degree")
    slope = (length / span) * per_degree * 180.0 / math.pi  # (Prestero 2001, eqs. 4.30, 4.32, p. 30)
    return function_from("slender_body_lift_slope", {"length": length, "span": span,
                                                    "body_lift_slope_per_degree": per_degree},
                         {"body_lift_slope": slope})


def fraction_of_length_pressure_centre():
    """``(length, pressure_centre_fraction) -> pressure_centre``: ``x_cp = l/2 - f l``, m, BODY, the nose at ``+l/2``
    (Prestero 2001, eq. 4.35, p. 31)."""
    length, fraction = ca.SX.sym("length"), ca.SX.sym("pressure_centre_fraction")
    return function_from("fraction_of_length_pressure_centre", {"length": length, "pressure_centre_fraction": fraction},
                         {"pressure_centre": length / 2.0 - fraction * length})  # (Prestero 2001, eq. 4.35, p. 31)


def point_force_body_lift():
    """``(nu_r, span, water_density, body_lift_slope, pressure_centre) -> tau``: the body lift of both planes as one point
    force at ``(x_cp, 0, 0)``, ``tau = [F; r x F]`` with ``F = [0, -k u v, -k u w]``, ``k = 1/2 rho d^2 c_ydb`` (module
    docstring; the moment has the sign of the cross product, not of the printed eq. 4.36)."""
    nu_r = ca.SX.sym("nu_r", 6)
    span, density = ca.SX.sym("span"), ca.SX.sym("water_density")
    slope, centre = ca.SX.sym("body_lift_slope"), ca.SX.sym("pressure_centre")
    gain = 0.5 * density * span**2 * slope  # k = 1/2 rho d^2 c_ydb (Prestero 2001, eq. 4.34, p. 31)
    force = ca.vertcat(0.0, -gain * nu_r[0] * nu_r[1], -gain * nu_r[0] * nu_r[2])  # (eqs. 4.33-4.34; Y_uvl = Z_uwl)
    tau = wrench(force, ca.vertcat(centre, 0.0, 0.0))  # (Fossen 2011, eq. 12.226, p. 400)
    return function_from("point_force_body_lift", {"nu_r": nu_r, "span": span, "water_density": density,
                                                   "body_lift_slope": slope, "pressure_centre": centre},
                         {"tau": tau})
