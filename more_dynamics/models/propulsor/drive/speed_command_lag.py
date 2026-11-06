"""Drive with a shaft-speed command followed by a first-order lag: the shaft speed is a state pulled to the command with the
time constant ``speed_time_constant``, ranged to ``[min_command, max_command]``; no thrust scale, no thrust limit, the load
torque ignored. The output ``drive_output`` is the ranged state in rev/s, the unit the shaft-speed laws accept.

    d(n)/dt = (n_c - n) / T_n,        n in [n_min, n_max]

``command_unit`` selects the unit of the command, of the state and of its limits: ``"rpm"`` or ``"rad/s"``. The command is
the target as it is given (it is not ranged before the lag); the range is applied to the state after the flow. A state in
a continuous model has no jump, so the range is carried by two things: the output is the ranged state, and the flow
that would push the state outward at a limit is zero. The state then rests on the limit while the command is beyond it
(as a state ranged after each step does, in the limit of a small step) and leaves it as soon as the command returns inside the
range; it is never wound up toward a command beyond the range. A step of the integrator that crosses a limit leaves the state
beyond it by at most one step of the flow, ``h |n_c - n| / T_n``; the output is the ranged state, so it never exceeds the limit.
The initial state must lie inside the range.

The form is the first-order actuator model of Fossen 2011 (eq. 13.154, p. 446, ``tau_dot = A_thr (tau - tau_com)`` with
``A_thr = -diag(1/T)``), written for the shaft speed instead of the generalised force of a thruster group.

Frame: scalar drive, no frame; the shaft turns positive when the thrust along the shaft axis is positive.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley.
    Ch. 13, eq. 13.154, p. 446 (first-order actuator dynamics; page read).

Author:    Enio Krizman
Date:      2026-11-05
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive, scalar
from more_dynamics.models.propulsor.drive.speed_command import COMMAND_UNITS, DRIVE_OUTPUT_UNIT, _check
from more_dynamics.models.shared.force_producer_common import saturate


def speed_command_lag_parameters(*, command_unit):
    _check(command_unit)
    return (scalar("min_command", command_unit, "lowest shaft speed (astern limit), also the limit of the command's effect"),
            scalar("max_command", command_unit, "highest shaft speed (ahead limit), also the limit of the command's effect"),
            positive("speed_time_constant", "s", "time constant T_n of the shaft-speed lag"))


def speed_command_lag_casadi(*, command_unit):
    _check(command_unit)
    declared = speed_command_lag_parameters(command_unit=command_unit)
    p = symbols(declared)
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 1), ca.SX.sym("shaft_torque")
    flow = (command[0] - state[0]) / p["speed_time_constant"]  # n_dot = (n_c - n) / T_n (Fossen 2011, eq. 13.154, p. 446)
    outward = ca.logic_or(ca.logic_and(state[0] >= p["max_command"], flow > 0.0),
                          ca.logic_and(state[0] <= p["min_command"], flow < 0.0))  # the range of the state, applied after the flow
    dot = ca.if_else(outward, 0.0, flow)
    speed = saturate(state[0], p["min_command"], p["max_command"]) / COMMAND_UNITS[command_unit]
    return part_function("speed_command_lag", {"command": command, "drive_state": state, "shaft_torque": load},
                         {"drive_output": speed, "drive_state_dot": dot, "thrust_scale": ca.SX(1.0),
                          "thrust_limit": ca.SX(ca.inf)}, p, declared)
