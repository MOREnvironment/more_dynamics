"""Drive with a shaft-speed command and a cap on the shaft power: above ``eta P_max`` the
thrust is scaled down by ``eta P_max / P``, and the thrust has its own limit.

The speed is as in ``speed_command``. The load torque ``Q`` of the law gives the shaft power
``P = |2 pi n Q|`` (definition of shaft power); ``thrust_scale = eta P_max / P`` when ``P``
exceeds ``eta P_max`` and 1 otherwise, the efficiency floored at 1e-3 (a guard of the
steerable outboard block this form comes from); the thrust limit is the declared
``thrust_saturation``. The power cap and the thrust cap are constructions of that block,
not published laws; they are kept as a labelled comparison form beside a drive that takes
the motor's own envelope.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 12, eq. 12.265, p. 411 (propeller torque and power).

Author:    Enio Krizman
Date:      2026-10-23
"""
import math

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from more_dynamics.models.propulsor.declarations import part_function, positive
from more_dynamics.models.propulsor.drive.speed_command import (COMMAND_UNITS, DRIVE_OUTPUT_UNIT,
                                                                speed_command_parameters)
from more_dynamics.models.shared.force_producer_common import saturate

MIN_EFFICIENCY = 1e-3  # floor of the efficiency in the power cap (guard of the block this form comes from)


def power_limited_speed_command_parameters(*, command_unit):
    return (*speed_command_parameters(command_unit=command_unit),
            positive("max_power", "W", "electrical input power limit P_max"),
            Parameter("efficiency", (1, 1), "1", "drive efficiency eta: the shaft power cap is eta P_max", 0.0, 1.0,
                      minimum_exclusive=True),
            positive("thrust_saturation", "N", "thrust limit of the unit"))


def power_limited_speed_command_casadi(*, command_unit):
    declared = power_limited_speed_command_parameters(command_unit=command_unit)
    p = symbols(declared)
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 0), ca.SX.sym("shaft_torque")
    speed = saturate(command, p["min_command"], p["max_command"]) / COMMAND_UNITS[command_unit]
    power = ca.fabs(2.0 * math.pi * speed * load)  # P = 2 pi n Q (definition of shaft power)
    limit = p["max_power"] * ca.fmax(p["efficiency"], MIN_EFFICIENCY)  # eta P_max, eta floored: construction
    scale = ca.if_else(power > limit, limit / ca.fmax(power, limit), 1.0)  # thrust scaled by limit / P above it: construction
    return part_function("power_limited_speed_command",
                         {"command": command, "drive_state": state, "shaft_torque": load},
                         {"drive_output": speed, "drive_state_dot": ca.SX.zeros(0, 1), "thrust_scale": scale,
                          "thrust_limit": p["thrust_saturation"]}, p, declared)
