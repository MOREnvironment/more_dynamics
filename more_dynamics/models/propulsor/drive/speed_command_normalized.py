"""Drive that takes a signed normalised command and maps it to a shaft speed through the
unit's own limits: ``+1`` is the highest speed ahead, ``-1`` the highest astern.

``n = u n_max`` for ``u >= 0`` and ``n = -u n_min`` for ``u < 0`` (``n_min`` negative, so
``u = -1`` gives ``n_min``), ``u`` clamped to [-1, 1]; the inverse of the allocator's
normalisation ``u = n / n_max`` ahead and ``n / |n_min|`` astern. No state, scale 1, no
thrust limit; the output is in rev/s, the unit the shaft-speed laws accept.

References
----------
Mapping defined by the normalisation of signed commands in the actuator table (a convention
of this library, no published source).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from more_dynamics.models.propulsor.declarations import part_function, positive

DRIVE_OUTPUT_UNIT = "rev/s"


def speed_command_normalized_parameters():
    return (positive("max_shaft_speed", "rev/s", "shaft speed at a command of +1 (ahead)"),
            Parameter("min_shaft_speed", (1, 1), "rev/s", "shaft speed at a command of -1 (astern), not above zero",
                      None, 0.0))


def speed_command_normalized_casadi():
    declared = speed_command_normalized_parameters()
    p = symbols(declared)
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 0), ca.SX.sym("shaft_torque")
    u = ca.fmax(ca.fmin(command[0], 1.0), -1.0)
    speed = ca.if_else(u >= 0.0, u * p["max_shaft_speed"], -u * p["min_shaft_speed"])  # inverse of the signed normalisation
    return part_function("speed_command_normalized", {"command": command, "drive_state": state, "shaft_torque": load},
                         {"drive_output": speed, "drive_state_dot": ca.SX.zeros(0, 1), "thrust_scale": ca.SX(1.0),
                          "thrust_limit": ca.SX(ca.inf)}, p, declared)
