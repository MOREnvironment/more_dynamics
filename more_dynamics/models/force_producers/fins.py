"""Rudder and stern planes: quadratic lift from deflection (CasADi).

``fins_parameters()`` declares the numbers (name, shape, SI unit, meaning,
admissible range); ``fins_casadi(convention=...)`` builds ``(delta, nu_r,
<parameters by name>) -> tau``.

Equations (keys in References), ``delta`` saturated at ``+-delta_max``
(MSS ``remus100.m`` 113-114):

* ``U_rh^2 = u_r^2 + v_r^2``, ``U_rv^2 = u_r^2 + w_r^2`` (``remus100.m``
  234-235);
* drag ``X_r = -1/2 rho U_rh^2 A_r CL_r delta_r^2``, ``X_s = -1/2 rho U_rv^2
  A_s CL_s delta_s^2`` (``remus100.m`` 238-239);
* lift ``Y_r = -1/2 rho U_rh^2 A_r CL_r delta_r``, ``Z_s = -1/2 rho U_rv^2 A_s
  CL_s delta_s`` (``remus100.m`` 242, 245) for the default convention; the
  sign of each lift follows ``convention``;
* moments of point forces on the centre line at ``x_r``, ``x_s``: ``tau =
  [f; r x f]`` (Fossen 2011, eq. 12.226, p. 400), i.e. ``M = -x_s Z_s``,
  ``N = x_r Y_r`` (``remus100.m`` 249-254).

A rudder or fin force ``F = k u`` with a quadratic dependence on speed is
Fossen 2011, eq. 12.225 and Table 12.3, p. 398.

``convention`` (a selector) names which deflection is positive for the
rudder (starboard or port trailing edge) and the stern planes (down or up);
it sets the signs of ``Y_r`` and ``Z_s`` (a construction of this module).
The default ``"starboard_down_positive"`` equals MSS ``remus100.m`` 113-114,
234-254. Deviation from that file: the two areas are separate parameters
(MSS sets both to ``2 S_fin``, ``remus100.m`` 183, 188).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 12,
    eq. 12.225, Table 12.3, p. 398; eq. 12.226, p. 400.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 113-114,
    179-189, 234-254.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from ._common import saturate, wrench

# convention -> (sign of Y_r, sign of Z_s); the default is remus100.m 242, 245
CONVENTIONS = {
    "starboard_down_positive": (-1.0, -1.0),
    "starboard_up_positive": (-1.0, 1.0),
    "port_down_positive": (1.0, -1.0),
    "port_up_positive": (1.0, 1.0),
}


def _positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


_DECLARED = (
    _positive("rudder_area", "m^2", "rudder area A_r (MSS: 2 S_fin)"),
    _positive("stern_plane_area", "m^2", "stern-plane area A_s (MSS: 2 S_fin)"),
    Parameter("rudder_lift_coefficient", (1, 1), "1/rad", "CL_delta_r"),
    Parameter("stern_plane_lift_coefficient", (1, 1), "1/rad", "CL_delta_s"),
    Parameter("rudder_position", (1, 1), "m", "x_r of the rudder from the CO"),
    Parameter("stern_plane_position", (1, 1), "m", "x_s of the stern planes from the CO"),
    _positive("max_deflection", "rad", "deflection limit delta_max"),
    _positive("water_density", "kg/m^3", "water density rho"),
)


def fins_parameters():
    """The declared parameter set: a tuple of ``Parameter`` (name, shape, SI
    unit, meaning, admissible range), in the order of the block's inputs."""
    return _DECLARED


def fins_casadi(*, convention="starboard_down_positive"):
    """The block: ``(delta, nu_r, <parameters by name>) -> tau``.

    Contract
    --------
    ``convention`` one of ``CONVENTIONS`` (unknown: ``ValueError``). Inputs:
    ``delta`` (2x1, rudder and stern-plane deflection in rad, saturated
    inside), ``nu_r`` (6x1, relative to the water) and the parameters of
    ``fins_parameters()`` by name. Output: ``tau`` (6x1, N and N m, BODY,
    about the CO).
    """
    if convention not in CONVENTIONS:
        raise ValueError(f"convention must be one of {tuple(CONVENTIONS)}, got {convention!r}")
    rudder_sign, stern_plane_sign = CONVENTIONS[convention]
    p = symbols(_DECLARED)
    delta = ca.SX.sym("delta", 2)
    nu_r = ca.SX.sym("nu_r", 6)
    limit = p["max_deflection"]
    delta_r = saturate(delta[0], -limit, limit)  # (remus100.m 113)
    delta_s = saturate(delta[1], -limit, limit)  # (remus100.m 114)
    u, v, w = nu_r[0], nu_r[1], nu_r[2]
    rudder_gain = 0.5 * p["water_density"] * p["rudder_area"] * p["rudder_lift_coefficient"]  # 1/2 rho A_r CL_r (remus100.m 242)
    stern_plane_gain = 0.5 * p["water_density"] * p["stern_plane_area"] * p["stern_plane_lift_coefficient"]  # (remus100.m 245)
    horizontal = rudder_gain * (u**2 + v**2)  # 1/2 rho U_rh^2 A_r CL_r (remus100.m 234, 242)
    vertical = stern_plane_gain * (u**2 + w**2)  # 1/2 rho U_rv^2 A_s CL_s (remus100.m 235, 245)

    x_r = -horizontal * delta_r**2  # (remus100.m 238)
    x_s = -vertical * delta_s**2  # (remus100.m 239)
    y_r = rudder_sign * horizontal * delta_r  # (remus100.m 242; sign by convention)
    z_s = stern_plane_sign * vertical * delta_s  # (remus100.m 245; sign by convention)
    # Point forces on the centreline: moments S(r) f give K = 0, M = -x_s Z_s, N = x_r Y_r
    # (Fossen 2011, eq. 12.226, p. 400; remus100.m 249-254)
    tau = wrench(ca.vertcat(x_r, y_r, 0.0), ca.vertcat(p["rudder_position"], 0.0, 0.0)) + wrench(
        ca.vertcat(x_s, 0.0, z_s), ca.vertcat(p["stern_plane_position"], 0.0, 0.0)
    )
    names = [d.name for d in _DECLARED]
    return ca.Function("fins", [delta, nu_r, *[p[name] for name in names]], [tau],
                       ["delta", "nu_r", *names], ["tau"])
