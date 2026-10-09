"""2-D drag law of a circular section, ``C_D(Re) kappa(L/D)``, the section law
of the cross-flow strip on a torpedo hull; the kinematic viscosity comes from
the ``site`` slot.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/cylinderDrag.m
    78-110 (circular section, Reynolds-dependent drag coefficient and
    aspect-ratio correction).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.cross_flow import cross_flow as cross_flow_module
from more_dynamics.models.shared.wiring import function_from


def circular_cylinder_reynolds():
    """``C_D(Re) kappa(L/D)``, viscosity taken from the ``site`` slot
    (cylinderDrag.m 78-110, through the module's own tables)."""
    nu_r = ca.SX.sym("nu_r", 6)
    length = ca.SX.sym("length")
    beam = ca.SX.sym("beam")
    kinematic_viscosity = ca.SX.sym("kinematic_viscosity")
    interp = cross_flow_module._clamped_interp_casadi
    speed = ca.sqrt(nu_r[1] ** 2 + nu_r[2] ** 2)  # (cylinderDrag.m 78)
    reynolds = speed * beam / kinematic_viscosity  # Re = v D / nu (cylinderDrag.m 79-80)
    kappa = ca.if_else(reynolds < cross_flow_module.CRITICAL_REYNOLDS_NUMBER,
                       interp(length / beam, cross_flow_module.KAPPA_SUBCRITICAL_DATA),
                       interp(length / beam, cross_flow_module.KAPPA_SUPERCRITICAL_DATA))  # (cylinderDrag.m 92-107)
    cd = interp(reynolds, cross_flow_module.CYLINDER_DRAG_DATA) * kappa  # (cylinderDrag.m 83-89, 110)
    return function_from("circular_cylinder_reynolds",
                         {"nu_r": nu_r, "length": length, "beam": beam, "kinematic_viscosity": kinematic_viscosity},
                         {"section_drag_coefficient": cd})


# Every input (nu_r, length, beam, kinematic_viscosity) is a kinematic input or a
# vehicle coupling (hull_form, site); none is a plugin-own parameter.
CIRCULAR_CYLINDER_REYNOLDS_PARAMETERS = ()
