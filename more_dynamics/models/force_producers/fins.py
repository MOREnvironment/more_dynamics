"""Rudder and stern planes: quadratic lift from deflection.

Ported from ``more_generic_models``
``dynamics/propulsion/fins/fins_auv_physical/fins_auv_physical.py``
(``FinsAUVPhysical``, lines 21-137) with ``fins_base.py``; the parameter
dataclass is ``config/dataclass/fins/fins_auv_physical_params.py``. The default
convention equals MSS ``remus100.m`` 114-115, 232-252. Unlike the source, the
two areas are separate parameters (the source sets both to ``2 S_fin``).
"""

from dataclasses import dataclass

import casadi as ca

from ._common import finite, positive, saturate

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
    horizontal = constants.rudder_gain * (u**2 + v**2)   # 0.5 rho U_rh^2 A_r CL_r
    vertical = constants.stern_plane_gain * (u**2 + w**2)  # 0.5 rho U_rv^2 A_s CL_s

    x_r = -horizontal * delta_r**2
    x_s = -vertical * delta_s**2
    y_r = constants.rudder_sign * horizontal * delta_r
    z_s = constants.stern_plane_sign * vertical * delta_s
    tau = ca.vertcat(
        x_r + x_s,
        y_r,
        z_s,
        0.0,
        -constants.stern_plane_position * z_s,
        constants.rudder_position * y_r,
    )
    return ca.Function("fins", [delta, nu_r], [tau], ["delta", "nu_r"], ["tau"])
