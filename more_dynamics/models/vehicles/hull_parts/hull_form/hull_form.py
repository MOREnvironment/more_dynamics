"""The ``hull_form`` slot: the shared dimensions every other part reads, and
the quantities their defining equations give (rule 16 — a derived quantity
is computed once, not retyped in every consuming part).

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 131-132,
    220-221 (one geometry: length and diameter give beam, draft, span and the
    semi-axes of the equivalent prolate spheroid); CRAFT/USV/models/otter.m
    92-93, 104-107 (twin-pontoon primitives).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.shared.wiring import function_from

PROLATE_SPHEROID_MAIN_DIMENSIONS_PARAMETERS = (
    Parameter("length", (1, 1), "m", "body length L", 0.0, minimum_exclusive=True),
    Parameter("diameter", (1, 1), "m", "body diameter D", 0.0, minimum_exclusive=True),
)
TWIN_PONTOONS_PARAMETERS = (
    Parameter("length", (1, 1), "m", "hull length L", 0.0, minimum_exclusive=True),
    Parameter("beam", (1, 1), "m", "overall beam B", 0.0, minimum_exclusive=True),
    Parameter("pontoon_beam", (1, 1), "m", "beam of one pontoon B_pont", 0.0, minimum_exclusive=True),
    Parameter("pontoon_lateral_offset", (1, 1), "m", "hull centreline y from the CO (each hull at +-y)"),
    Parameter("pontoon_block_coefficient", (1, 1), "1", "block coefficient Cb of one pontoon",
             0.0, 1.0, minimum_exclusive=True),
    Parameter("pontoon_waterplane_coefficient", (1, 1), "1", "waterplane coefficient Cw of one pontoon",
             0.0, 1.0, minimum_exclusive=True),
)


def prolate_spheroid_main_dimensions():
    """``(length, diameter) -> (length, beam, draft, span, section_beam,
    semi_major_axis, semi_minor_axis)``; beam = draft = span = section_beam =
    diameter (remus100.m 220-221: the cross-flow strip's section is the
    circle of the body's diameter); the semi-axes of the equivalent spheroid
    are length/2, diameter/2."""
    length = ca.SX.sym("length")
    diameter = ca.SX.sym("diameter")
    return function_from("prolate_spheroid_main_dimensions", {"length": length, "diameter": diameter}, {
        "length": length, "beam": diameter, "draft": diameter, "span": diameter,
        "section_beam": diameter, "semi_major_axis": length / 2, "semi_minor_axis": diameter / 2,
    })


def twin_pontoons():
    """Primitives only, unchanged (otter.m 92-93, 104-107): length, the
    overall beam, the pontoon beam, its lateral offset, and the pontoons'
    block and waterplane coefficients; ``section_beam``, the beam of one
    cross-flow strip section, is the beam of one pontoon (otter.m 245)."""
    names = ("length", "beam", "pontoon_beam", "pontoon_lateral_offset",
              "pontoon_block_coefficient", "pontoon_waterplane_coefficient")
    s = {n: ca.SX.sym(n) for n in names}
    return function_from("twin_pontoons", s, {**s, "section_beam": s["pontoon_beam"]})
