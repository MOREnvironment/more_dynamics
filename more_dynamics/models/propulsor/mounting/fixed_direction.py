"""Mounting of a fixed thruster given by a direction: the axis is the declared direction
divided by its length (any length is accepted, a zero direction is refused); the point is the
declared position; the steering angle is ignored.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 12, eq. 12.226, p. 400 (unit thrust direction at a point).
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/USV/models/otter.m 242
    @ ed2ce7b.

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import check_values, symbols

from more_dynamics.models.propulsor.declarations import part_function, vector


def fixed_direction_parameters():
    return (vector("position", "m", "CO to the thrust point, BODY (FRD)"),
            vector("direction", "1", "thrust direction, BODY (FRD), normalised inside"))


def check_fixed_direction_values(values):
    """The numbers checked; a zero direction is refused, naming it."""
    numbers = check_values(fixed_direction_parameters(), values)
    if not float(ca.norm_2(numbers["direction"])) > 0.0:
        raise ValueError("parameter 'direction' [1]: must be non-zero")
    return numbers


def fixed_direction_casadi():
    declared = fixed_direction_parameters()
    p = symbols(declared)
    angle = ca.SX.sym("steering_angle")
    axis = p["direction"] / ca.norm_2(p["direction"])  # (Fossen 2011, eq. 12.226, p. 400; otter.m 242)
    return part_function("fixed_direction", {"steering_angle": angle},
                         {"thrust_axis": axis, "thrust_position": p["position"]}, p, declared)
