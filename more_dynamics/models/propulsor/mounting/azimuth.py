"""Mounting of a unit that turns about the body z axis: the axis is ``[cos d, s sin d, 0]`` and
the point the declared position, ``s`` set by ``steering_convention``.

``"port_positive"`` (``s = +1``): a positive angle turns the thrust force to starboard (+y),
so for a unit aft of the CO it turns the bow to port. ``"starboard_positive"`` (``s = -1``):
the name says which way the boat turns for a unit aft of the CO, so a positive angle turns the
bow to starboard. Which way the boat turns depends on the lever arm, so a unit at or ahead of
the CO with the convention is refused by ``check_azimuth_values`` (the meaning of the name
there is not decided).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 12, Table 12.3, p. 398; p. 400 below eq. 12.226 (azimuth thruster).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import check_values, symbols

from more_dynamics.models.propulsor.declarations import part_function, vector

STEERING_CONVENTIONS = {"port_positive": 1.0, "starboard_positive": -1.0}


def _check(steering_convention):
    if steering_convention not in STEERING_CONVENTIONS:
        raise ValueError(f"azimuth: steering_convention must be one of {sorted(STEERING_CONVENTIONS)}, "
                         f"got {steering_convention!r}")


def azimuth_parameters(*, steering_convention="port_positive"):
    _check(steering_convention)
    return (vector("position", "m", "CO to the thrust point, BODY (FRD)"),)


def check_azimuth_values(values, *, steering_convention="port_positive"):
    """The numbers checked; a thrust point at or ahead of the CO is refused (the boat's turn
    named by the convention is defined for a unit aft of it)."""
    numbers = check_values(azimuth_parameters(steering_convention=steering_convention), values)
    if not float(numbers["position"][0]) < 0.0:
        raise ValueError(f"parameter 'position' [m]: steering_convention {steering_convention!r} names the boat's turn "
                         "for a unit aft of the CO; a thrust point at or ahead of it is refused")
    return numbers


def azimuth_casadi(*, steering_convention="port_positive"):
    _check(steering_convention)
    declared = azimuth_parameters(steering_convention=steering_convention)
    p = symbols(declared)
    angle = ca.SX.sym("steering_angle")
    axis = ca.vertcat(ca.cos(angle), STEERING_CONVENTIONS[steering_convention] * ca.sin(angle), 0.0)  # (Fossen 2011, Table 12.3, p. 398)
    return part_function("azimuth", {"steering_angle": angle}, {"thrust_axis": axis, "thrust_position": p["position"]},
                         p, declared)
