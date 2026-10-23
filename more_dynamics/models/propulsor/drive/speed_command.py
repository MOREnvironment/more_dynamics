"""Drive that takes a shaft-speed command, limits it and turns it into revolutions
per second; no state, no thrust scale, no thrust limit, the load torque ignored.

``command_unit`` selects the unit of the command and of its limits: ``"rpm"``
(``n_p = sat(n, n_max) / 60``, MSS ``remus100.m`` 116) or ``"rad/s"``
(``n = sat(omega, omega_max) / (2 pi)``, the shaft-speed saturation of MSS
``otter.m`` 230). The output ``drive_output`` is in rev/s, the unit the shaft-speed
laws accept.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 116
    and CRAFT/USV/models/otter.m 220 @ 99fdf4d.

Author:    Enio Krizman
Date:      2026-10-23
"""
import math

import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, scalar
from more_dynamics.models.shared.force_producer_common import saturate

DRIVE_OUTPUT_UNIT = "rev/s"
COMMAND_UNITS = {"rpm": 60.0, "rad/s": 2.0 * math.pi}  # command units per revolution per second


def _check(command_unit):
    if command_unit not in COMMAND_UNITS:
        raise ValueError(f"speed_command: command_unit must be one of {sorted(COMMAND_UNITS)}, got {command_unit!r}")


def speed_command_parameters(*, command_unit):
    _check(command_unit)
    return (scalar("min_command", command_unit, "lowest shaft-speed command (astern limit)"),
            scalar("max_command", command_unit, "highest shaft-speed command (ahead limit)"))


def speed_command_casadi(*, command_unit):
    _check(command_unit)
    declared = speed_command_parameters(command_unit=command_unit)
    p = symbols(declared)
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 0), ca.SX.sym("shaft_torque")
    speed = saturate(command, p["min_command"], p["max_command"]) / COMMAND_UNITS[command_unit]  # (remus100.m 116; otter.m 220)
    return part_function("speed_command", {"command": command, "drive_state": state, "shaft_torque": load},
                         {"drive_output": speed, "drive_state_dot": ca.SX.zeros(0, 1), "thrust_scale": ca.SX(1.0),
                          "thrust_limit": ca.SX(ca.inf)}, p, declared)
