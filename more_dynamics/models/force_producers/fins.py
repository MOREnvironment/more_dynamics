"""Rudder and stern planes: quadratic lift from deflection.

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

The default convention equals MSS ``remus100.m`` 113-114, 234-254. Deviation
from that file: the two areas are separate parameters (MSS and the numpy
source set both to ``2 S_fin``, ``remus100.m`` 183, 188).

Ported from the numpy source [MGM]
``dynamics/propulsion/fins/fins_auv_physical/fins_auv_physical.py``
(``FinsAUVPhysical``, 21-137) with ``fins_base.py``; the parameter dataclass is
``config/dataclass/fins/fins_auv_physical_params.py``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 12, eq. 12.225, Table 12.3, p. 398; eq. 12.226, p. 400.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``CRAFT/AUV/models/remus100.m`` 113-114, 179-189, 234-254.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``:
    the files and lines listed above.
"""

from dataclasses import dataclass

import casadi as ca

from ._common import finite, positive, saturate, wrench

# convention -> (sign of Y_r, sign of Z_s); fins_auv_physical.py lines 87-121
CONVENTIONS = {
    "starboard_down_positive": (-1.0, -1.0),
    "starboard_up_positive": (-1.0, 1.0),
    "port_down_positive": (1.0, -1.0),
    "port_up_positive": (1.0, 1.0),
}


@dataclass(frozen=True)
class FinsConstants:
    """``*_gain = 0.5 rho A CL``; forces are ``gain U^2 delta`` (and ``delta^2``)."""

    rudder_gain: float
    stern_plane_gain: float
    rudder_position: float
    stern_plane_position: float
    max_deflection: float
    rudder_sign: float
    stern_plane_sign: float


def preprocess_fins(
    rudder_area: float,
    stern_plane_area: float,
    rudder_lift_coefficient: float,
    stern_plane_lift_coefficient: float,
    rudder_position: float,
    stern_plane_position: float,
    max_deflection: float,
    water_density: float,
    convention: str = "starboard_down_positive",
) -> FinsConstants:
    """Source names: ``Ar``, ``Ae``, ``CL_delta_r``, ``CL_delta_e``, ``x_r``,
    ``x_e`` (m), ``delta_max`` (here in rad), ``rho``, ``convention``."""
    if convention not in CONVENTIONS:
        raise ValueError(f"convention must be one of {tuple(CONVENTIONS)}")
    rho = positive("water_density", water_density)
    rudder_area = positive("rudder_area", rudder_area)
    stern_plane_area = positive("stern_plane_area", stern_plane_area)
    rudder_sign, stern_plane_sign = CONVENTIONS[convention]
    return FinsConstants(
        rudder_gain=0.5 * rho * rudder_area * finite("rudder_lift_coefficient", rudder_lift_coefficient),
        stern_plane_gain=0.5 * rho * stern_plane_area
        * finite("stern_plane_lift_coefficient", stern_plane_lift_coefficient),
        rudder_position=finite("rudder_position", rudder_position),
        stern_plane_position=finite("stern_plane_position", stern_plane_position),
        max_deflection=positive("max_deflection", max_deflection),
        rudder_sign=rudder_sign,
        stern_plane_sign=stern_plane_sign,
    )


def fins_casadi(constants: FinsConstants) -> ca.Function:
    """``(delta, nu_r) -> tau``; ``delta = (rudder, stern plane)`` rad, saturated."""
    delta = ca.SX.sym("delta", 2)
    nu_r = ca.SX.sym("nu_r", 6)
    limit = constants.max_deflection
    delta_r = saturate(delta[0], -limit, limit)
    delta_s = saturate(delta[1], -limit, limit)
    u, v, w = nu_r[0], nu_r[1], nu_r[2]
    horizontal = constants.rudder_gain * (u**2 + v**2)   # 0.5 rho U_rh^2 A_r CL_r (remus100.m 234)
    vertical = constants.stern_plane_gain * (u**2 + w**2)  # 0.5 rho U_rv^2 A_s CL_s (remus100.m 235)

    x_r = -horizontal * delta_r**2  # (remus100.m 238)
    x_s = -vertical * delta_s**2  # (remus100.m 239)
    y_r = constants.rudder_sign * horizontal * delta_r  # (remus100.m 242)
    z_s = constants.stern_plane_sign * vertical * delta_s  # (remus100.m 245)
    # Point forces on the centreline: moments S(r) f give K = 0, M = -x_s Z_s, N = x_r Y_r
    # (Fossen 2011, eq. 12.226, p. 400; remus100.m 249-254)
    tau = wrench(ca.vertcat(x_r, y_r, 0.0), [constants.rudder_position, 0.0, 0.0]) + wrench(
        ca.vertcat(x_s, 0.0, z_s), [constants.stern_plane_position, 0.0, 0.0]
    )
    return ca.Function("fins", [delta, nu_r], [tau], ["delta", "nu_r"], ["tau"])
