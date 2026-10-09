"""2-D drag law of a rectangular section from Hoerner's table, ``Cd_2D(B / 2T)``,
the section law of the cross-flow strip on a pontoon.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/Hoerner.m 47-51
    (rectangular section, Hoerner's table).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.cross_flow import cross_flow as cross_flow_module
from more_dynamics.models.shared.wiring import function_from


def rectangular_section_hoerner():
    """Hoerner's rectangular-section table, ``Cd_2D(B / 2T)`` (Hoerner.m
    47-51)."""
    beam = ca.SX.sym("beam")
    draft = ca.SX.sym("draft")
    cd = cross_flow_module._clamped_interp_casadi(beam / (2 * draft), cross_flow_module.HOERNER_DRAG_DATA)
    return function_from("rectangular_section_hoerner", {"beam": beam, "draft": draft},
                         {"section_drag_coefficient": cd})


# Every input (beam, draft) is a vehicle coupling (hull_form); none is a
# plugin-own parameter.
RECTANGULAR_SECTION_HOERNER_PARAMETERS = ()
