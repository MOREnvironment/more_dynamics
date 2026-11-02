"""Housing of an outboard unit that adds the normal force of its leg: a flat-plate-like force on the leg,
proportional to the flow component normal to the motor axis and to the flow speed at the leg, acting at the
thrust point; the thrust-vectoring unit without the leg is the housing ``none``.

The force and its slope ``C_L' = 6.13 AR / (AR + 2.25)`` are computed by
``more_dynamics.models.outboard.outboard_motor.outboard_leg_normal_force_about_normal``; the slope is the rudder lift
slope of Whicker and Fehlner adapted to the leg, a construction (not a measurement of a leg). The
density is the vehicle's, fed to the housing.

The leg lies in the vertical plane through the thrust axis, so its normal is taken from the thrust axis the mounting
gives, ``n = [-a_y, a_x] / sqrt(a_x^2 + a_y^2)`` (the horizontal part of the axis turned by +90 deg), not from the
steering angle: the axis already carries the steering convention of the mounting (``port_positive``: ``a = [cos d,
sin d]``, ``n = [-sin d, cos d]``; ``starboard_positive``: ``a = [cos d, -sin d]``, ``n = [sin d, cos d]``), so the leg
and the thrust axis agree in both. The steering angle input is not read. A vertical axis (no horizontal part) gives no
leg force.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control.
    Wiley. Ch. 12, eq. 12.226, p. 400.
[Krizman] Krizman, E. Model file ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)``, MSS fork branch ``grethe-vehicle`` @ 679a6e8,
    CRAFT/USV/models/grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 269-271, 400-408 (read; the file may still be revised).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.outboard.outboard_motor import outboard_leg_normal_force_about_normal
from more_dynamics.models.propulsor.declarations import part_function, positive

AXIS_FLOOR = 1e-12  # below this the axis has no horizontal part and the leg gives no force


def leg_normal_force_parameters():
    return (positive("leg_area", "m^2", "leg and fin area A_leg"),
            positive("leg_aspect_ratio", "1", "leg aspect ratio AR"))


def leg_normal_force_casadi():
    declared = leg_normal_force_parameters()
    p = symbols(declared)
    inputs = {"nu_r": ca.SX.sym("nu_r", 6), "thrust_axis": ca.SX.sym("thrust_axis", 3),
              "thrust_position": ca.SX.sym("thrust_position", 3), "steering_angle": ca.SX.sym("steering_angle"),
              "water_density": ca.SX.sym("water_density")}
    axis = inputs["thrust_axis"]
    horizontal = ca.sqrt(axis[0]**2 + axis[1]**2)  # |horizontal part of the axis|
    has_axis = horizontal > AXIS_FLOOR
    scale = ca.if_else(has_axis, 1.0 / ca.if_else(has_axis, horizontal, 1.0), 0.0)
    normal = (-axis[1] * scale, axis[0] * scale)  # the horizontal axis turned by +90 deg, unit length
    tau = outboard_leg_normal_force_about_normal(inputs["nu_r"], normal, inputs["thrust_position"],
                                                 inputs["water_density"], p["leg_area"], p["leg_aspect_ratio"])
    return part_function("leg_normal_force", inputs, {"tau": tau}, p, declared)
