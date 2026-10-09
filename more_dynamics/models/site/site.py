"""The ``site`` slot: one water density, one gravity and one kinematic
viscosity, fed to every other part by the vehicle (one quantity, one value).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. WGS-84 normal gravity, through
    more_transformations' ``ECEFNEDtransform.gravity`` (cited there).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: INS/functions/gravity.m 11-12
    (latitude form; remus100.m 96-97); CRAFT/USV/models/otter.m 90-91 (given
    form).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import ECEFNEDtransform, Parameter

from more_dynamics.models.shared.wiring import function_from

_WATER_DENSITY = Parameter("water_density", (1, 1), "kg/m^3", "water density rho", 0.0, minimum_exclusive=True)
_KINEMATIC_VISCOSITY = Parameter("kinematic_viscosity", (1, 1), "m^2/s", "kinematic viscosity nu",
                                 0.0, minimum_exclusive=True)

SITE_AT_LATITUDE_PARAMETERS = (
    Parameter("latitude", (1, 1), "rad", "geodetic latitude (WGS-84 normal gravity)", -1.5708, 1.5708),
    _WATER_DENSITY, _KINEMATIC_VISCOSITY,
)
SITE_GIVEN_PARAMETERS = (
    Parameter("gravity", (1, 1), "m/s^2", "acceleration of gravity g", 0.0, minimum_exclusive=True),
    _WATER_DENSITY, _KINEMATIC_VISCOSITY,
)


def site_at_latitude():
    """``(latitude, water_density, kinematic_viscosity) -> (water_density,
    gravity, kinematic_viscosity)``; gravity from WGS-84 normal gravity at
    the latitude (gravity.m 11-12; remus100.m 96-97)."""
    latitude = ca.SX.sym("latitude")
    water_density = ca.SX.sym("water_density")
    kinematic_viscosity = ca.SX.sym("kinematic_viscosity")
    gravity = ECEFNEDtransform.gravity(latitude, unit="rad")
    return function_from(
        "site_at_latitude",
        {"latitude": latitude, "water_density": water_density, "kinematic_viscosity": kinematic_viscosity},
        {"water_density": water_density, "gravity": gravity, "kinematic_viscosity": kinematic_viscosity},
    )


def site_given():
    """``(gravity, water_density, kinematic_viscosity) -> (water_density,
    gravity, kinematic_viscosity)``; every value given (otter.m 90-91)."""
    gravity = ca.SX.sym("gravity")
    water_density = ca.SX.sym("water_density")
    kinematic_viscosity = ca.SX.sym("kinematic_viscosity")
    return function_from(
        "site_given",
        {"gravity": gravity, "water_density": water_density, "kinematic_viscosity": kinematic_viscosity},
        {"water_density": water_density, "gravity": gravity, "kinematic_viscosity": kinematic_viscosity},
    )
