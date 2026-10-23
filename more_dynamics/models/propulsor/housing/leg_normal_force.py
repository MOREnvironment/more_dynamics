"""Housing of an outboard unit that adds the normal force of its leg: a flat-plate-like force on the leg,
proportional to the flow component normal to the motor axis and to the flow speed at the leg, acting at the
thrust point; the thrust-vectoring unit without the leg is the housing ``none``.

The force and its slope ``C_L' = 6.13 AR / (AR + 2.25)`` are computed by
``more_dynamics.models.outboard.outboard_motor.outboard_leg_normal_force``; the slope is the rudder lift
slope of Whicker and Fehlner adapted to the leg, a construction (not a measurement of a leg). The
density is the vehicle's, fed to the housing; the steering angle is the unit's actual deflection.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control.
    Wiley. Ch. 12, eq. 12.226, p. 400.
[Krizman] Krizman, E. Model file ``grethe.m``, MSS fork branch ``grethe-vehicle`` @ 01cd570,
    CRAFT/USV/models/grethe.m 265-268, 394-404 (read; the file may still be revised).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.outboard.outboard_motor import outboard_leg_normal_force
from more_dynamics.models.propulsor.declarations import part_function, positive


def leg_normal_force_parameters():
    return (positive("leg_area", "m^2", "leg and fin area A_leg"),
            positive("leg_aspect_ratio", "1", "leg aspect ratio AR"))


def leg_normal_force_casadi():
    declared = leg_normal_force_parameters()
    p = symbols(declared)
    inputs = {"nu_r": ca.SX.sym("nu_r", 6), "thrust_axis": ca.SX.sym("thrust_axis", 3),
              "thrust_position": ca.SX.sym("thrust_position", 3), "steering_angle": ca.SX.sym("steering_angle"),
              "water_density": ca.SX.sym("water_density")}
    tau = outboard_leg_normal_force(inputs["nu_r"], inputs["steering_angle"], inputs["thrust_position"],
                                    inputs["water_density"], p["leg_area"], p["leg_aspect_ratio"])
    return part_function("leg_normal_force", inputs, {"tau": tau}, p, declared)
