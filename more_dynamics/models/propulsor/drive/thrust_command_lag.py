"""Drive that takes a thrust command in newtons and delivers it through a first-order lag: the actual thrust is a
state, pulled to the command with the time constant ``thrust_time_constant``; no thrust scale, no thrust limit, the
load torque ignored. The output ``drive_output`` is the actual thrust in newtons (unit ``"N"``).

    d(thrust)/dt = (command - thrust) / T_thrust

The lag of the thrust of the outboard of ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)``: its ``lim.T_thrust`` is a time constant the caller applies to
the thrust command (``SIMgrethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)`` 257, ``Th = Th + h/T_T (thrust_c - Th)``); here it is the continuous form of that
state, integrated with the vehicle's own RK4.

References
----------
[Krizman] Krizman, E. Model file ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)`` and its simulator ``SIMgrethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)``, MSS fork branch ``grethe-vehicle``
    @ 679a6e8, CRAFT/USV/models/grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 268 (``T_thrust``), CRAFT/USV/SIMgrethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 151-154, 257-259 (read, the
    files may still be revised).

Author:    Enio Krizman
Date:      2026-10-26
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive

DRIVE_OUTPUT_UNIT = "N"


def thrust_command_lag_parameters():
    return (positive("thrust_time_constant", "s", "time constant of the thrust's first-order lag"),)


def thrust_command_lag_casadi():
    declared = thrust_command_lag_parameters()
    p = symbols(declared)
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 1), ca.SX.sym("shaft_torque")
    return part_function("thrust_command_lag", {"command": command, "drive_state": state, "shaft_torque": load},
                         {"drive_output": state[0], "drive_state_dot": (command - state) / p["thrust_time_constant"],
                          "thrust_scale": ca.SX(1.0), "thrust_limit": ca.SX(ca.inf)}, p, declared)
