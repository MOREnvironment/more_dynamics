"""The cross-flow strip integral, one member of the ``hydrodynamic_loads``
list, with the 2-D section drag law as its own named sub-slot
(``section_drag``) — the same by-name resolution the vehicle uses,
one level down. One part file serves both the
torpedo (circular section) and the catamaran (Hoerner section).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 6.91-6.92, p. 127 (strip integral).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/crossFlowDrag.m
    54-69 (strip integral); LIBRARY/modeling/cylinderDrag.m 78-110 (circular
    section, Reynolds-dependent drag coefficient and aspect-ratio
    correction); LIBRARY/modeling/Hoerner.m 47-51 (rectangular section,
    Hoerner's table).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow import cross_flow as cross_flow_module

from more_dynamics.models.wiring import function_from


def _strip_positions(length):
    dx = length / cross_flow_module.NUMBER_OF_STRIPS  # (crossFlowDrag.m 38)
    n = cross_flow_module.NUMBER_OF_STRIPS
    return [-length / 2 + (i - 0.5) * dx for i in range(1, n + 1)], dx  # (crossFlowDrag.m 56)


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


def rectangular_section_hoerner():
    """Hoerner's rectangular-section table, ``Cd_2D(B / 2T)`` (Hoerner.m
    47-51)."""
    beam = ca.SX.sym("beam")
    draft = ca.SX.sym("draft")
    cd = cross_flow_module._clamped_interp_casadi(beam / (2 * draft), cross_flow_module.HOERNER_DRAG_DATA)
    return function_from("rectangular_section_hoerner", {"beam": beam, "draft": draft},
                         {"section_drag_coefficient": cd})


def cross_flow_strip_with_section(section):
    """The strip integral (crossFlowDrag.m 54-69), ``section`` resolved as
    a named sub-slot: every one of its inputs not already in the strip's
    own pool becomes an input of the composed part."""
    nu_r = ca.SX.sym("nu_r", 6)
    length = ca.SX.sym("length")
    draft = ca.SX.sym("draft")
    water_density = ca.SX.sym("water_density")
    pool = {"nu_r": nu_r, "length": length, "draft": draft}
    extra = {}
    for name in section.name_in():
        if name not in pool:
            extra[name] = pool[name] = ca.SX.sym(name, section.size_in(name))
    section_drag_coefficient = section(**{n: pool[n] for n in section.name_in()})["section_drag_coefficient"]
    positions, dx = _strip_positions(length)
    factor = -0.5 * water_density * draft * section_drag_coefficient * dx  # (crossFlowDrag.m 63-66)
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    y_sum = z_sum = m_sum = n_sum = ca.SX(0.0)
    for x in positions:
        horizontal = v_r + x * r  # (crossFlowDrag.m 61)
        vertical = w_r + x * q  # (crossFlowDrag.m 62)
        y_sum += ca.fabs(horizontal) * horizontal  # (crossFlowDrag.m 63)
        z_sum += ca.fabs(vertical) * vertical  # (crossFlowDrag.m 64)
        m_sum += x * ca.fabs(vertical) * vertical  # (crossFlowDrag.m 65)
        n_sum += x * ca.fabs(horizontal) * horizontal  # (crossFlowDrag.m 66)
    tau = ca.vertcat(0, factor * y_sum, factor * z_sum, 0, factor * m_sum, factor * n_sum)  # (crossFlowDrag.m 69)
    return function_from(f"cross_flow_strip_{section.name()}",
                         {"nu_r": nu_r, "length": length, "draft": draft, "water_density": water_density, **extra},
                         {"tau": tau})


def cross_flow_strip_circular_cylinder_reynolds():
    """The shared strip integral with the circular (torpedo) section."""
    return cross_flow_strip_with_section(circular_cylinder_reynolds())


def cross_flow_strip_rectangular_section_hoerner():
    """The shared strip integral with the rectangular Hoerner (catamaran)
    section."""
    return cross_flow_strip_with_section(rectangular_section_hoerner())


# Every remaining input of every part here (nu_r, length, beam, draft,
# water_density, kinematic_viscosity) is a kinematic input or a vehicle
# coupling (hull_form, site); none is a plugin-own parameter.
CIRCULAR_CYLINDER_REYNOLDS_PARAMETERS = ()
RECTANGULAR_SECTION_HOERNER_PARAMETERS = ()
CROSS_FLOW_STRIP_PARAMETERS = ()
