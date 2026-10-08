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
    54-69 (strip integral).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow import cross_flow as cross_flow_module
from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.section_drag.circular_cylinder_reynolds import (
    circular_cylinder_reynolds)
from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.section_drag.rectangular_section_hoerner import (
    rectangular_section_hoerner)
from more_dynamics.models.shared.wiring import function_from


def _strip_positions(length):
    dx = length / cross_flow_module.NUMBER_OF_STRIPS  # (crossFlowDrag.m 38)
    n = cross_flow_module.NUMBER_OF_STRIPS
    return [-length / 2 + (i - 0.5) * dx for i in range(1, n + 1)], dx  # (crossFlowDrag.m 56)


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


# Every remaining input (nu_r, length, draft, water_density) is a kinematic input
# or a vehicle coupling (hull_form, site); none is a plugin-own parameter.
CROSS_FLOW_STRIP_PARAMETERS = ()
