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

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.cross_flow import cross_flow as cross_flow_module
from more_dynamics.models.cross_flow.section_drag.circular_cylinder_reynolds import (
    circular_cylinder_reynolds)
from more_dynamics.models.cross_flow.section_drag.rectangular_section_hoerner import (
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


COMBINED_DRAG_CROSS_FLOW_PARAMETERS = (
    Parameter("cross_flow_drag_coefficients", (2, 1), "1", "[Cdy, Cdz] of the two cross-flow planes (npsauv.m 220)"),
    Parameter("cross_flow_section_dimensions", (2, 1), "m",
             "[Hx, Bx] average submerged height and width (npsauv.m 220-221)"),
    Parameter("cross_flow_sections", (1, 1), "1", "number of strips of the integral (npsauv.m 218)",
             1.0, minimum_exclusive=True),
)


def cross_flow_strip_combined_drag(*, sections):
    """The strip integral of a body whose two cross-flow planes (horizontal,
    vertical) share one combined, speed-normalised drag term per strip
    (``npsauv.m`` 217-235), rather than ``cross_flow_strip_with_section``'s
    separate per-plane section-drag sub-slot:

        Ucf = sqrt((v_r + x r)^2 + (w_r - x q)^2) + eps
        drag = Cdy Hx (v_r + x r)^2 + Cdz Bx (w_r - x q)^2
        Cy += dx drag (v_r + x r) / Ucf,  Cz += dx drag (w_r - x q) / Ucf
        Cm += dx drag (w_r + x q) / Ucf x,  Cn += dx drag (v_r + x r) / Ucf x
        tau = (rho / 2) [0, -Cy, -Cz, 0, -Cm, -Cn]

    ``sections`` (an int, not a CasADi symbol: it sets how many strip terms
    the graph unrolls) is read once at build time; the declared parameter
    ``cross_flow_sections`` of the same value is checked against it by the
    plugin, not fed into this graph.

    **Reproduced as written, not fixed (rule 15):** ``Cm``'s strip term uses
    ``(w_r + x q)``, the same sign as the vertical speed term rather than
    ``(w_r - x q)`` that ``Cz`` and the cross-flow kinematics use everywhere
    else in this file (``npsauv.m`` 231); a candidate MSS asymmetry, flagged
    for A-65's register, not corrected here."""
    nu_r = ca.SX.sym("nu_r", 6)
    length = ca.SX.sym("length")
    water_density = ca.SX.sym("water_density")
    cdy_cdz = ca.SX.sym("cross_flow_drag_coefficients", 2)
    hx_bx = ca.SX.sym("cross_flow_section_dimensions", 2)
    cdy, cdz = cdy_cdz[0], cdy_cdz[1]
    hx, bx = hx_bx[0], hx_bx[1]
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    dx = length / sections  # (npsauv.m 218)
    eps = 1e-6  # (npsauv.m 219, avoid dividing by zero)
    y_sum = z_sum = m_sum = n_sum = ca.SX(0.0)
    for i in range(sections + 1):  # (npsauv.m 224: for xL = -L/2 : dxL : L/2, sections + 1 points)
        x = -length / 2 + i * dx
        ucf = ca.sqrt((v_r + x * r) ** 2 + (w_r - x * q) ** 2) + eps  # (npsauv.m 226)
        drag = cdy * hx * (v_r + x * r) ** 2 + cdz * bx * (w_r - x * q) ** 2  # (npsauv.m 227)
        y_sum += dx * drag * (v_r + x * r) / ucf  # (npsauv.m 229)
        z_sum += dx * drag * (w_r - x * q) / ucf  # (npsauv.m 230)
        m_sum += dx * drag * (w_r + x * q) / ucf * x  # (npsauv.m 231, the flagged sign, module docstring)
        n_sum += dx * drag * (v_r + x * r) / ucf * x  # (npsauv.m 232)
    tau = (water_density / 2.0) * ca.vertcat(0, -y_sum, -z_sum, 0, -m_sum, -n_sum)  # (npsauv.m 235)
    return function_from(
        "cross_flow_strip_combined_drag",
        {"nu_r": nu_r, "length": length, "water_density": water_density,
         "cross_flow_drag_coefficients": cdy_cdz, "cross_flow_section_dimensions": hx_bx},
        {"tau": tau})
