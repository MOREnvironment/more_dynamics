"""Drive that takes a thrust command in newtons and hands it on as it is: no state, no thrust scale, no thrust limit,
the load torque ignored. The thrust is limited where the law knows the limit (a law that caps it by the power
delivered and the speed of the water).

The output ``drive_output`` is a thrust in newtons (unit ``"N"``), the unit the thrust-command laws accept. This is
the form of a unit commanded by its thrust, whose actuator dynamics are the caller's, as ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)`` takes its
thrust ``ui(1)`` as the actual thrust (``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)`` 38-48, "the lag ... applied by the caller");
``thrust_command_lag`` is the same with the lag as a state of the unit.

References
----------
[Krizman] Krizman, E. Model file ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)``, MSS fork branch ``grethe-vehicle`` @ 679a6e8,
    CRAFT/USV/models/grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 38-48, 391-393 (the thrust input and its limit; read, the file may still be revised).

Author:    Enio Krizman
Date:      2026-10-26
"""
import casadi as ca

from more_dynamics.models.propulsor.declarations import part_function

DRIVE_OUTPUT_UNIT = "N"


def thrust_command_parameters():
    return ()


def thrust_command_casadi():
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 0), ca.SX.sym("shaft_torque")
    return part_function("thrust_command", {"command": command, "drive_state": state, "shaft_torque": load},
                         {"drive_output": command[0], "drive_state_dot": ca.SX.zeros(0, 1),
                          "thrust_scale": ca.SX(1.0), "thrust_limit": ca.SX(ca.inf)}, {}, ())
