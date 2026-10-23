"""Drive on a normalised effort: the effort is a state clamped to [-1, 1], pulled to the
clamped command by a first-order lag with separate rise and fall time constants; the output
is the clamped state (unit 1, an effort, not a shaft speed).

``d(effort)/dt = (ref - effort) / T_rise`` when the reference is above the effort and
``/ T_fall`` otherwise, ``ref`` and ``effort`` clamped to [-1, 1]; the form of the thrust
state of the shared repository's jet nozzle plugin (``jet_nozzle.py`` 58-73), reproduced as
a labelled comparison form.

References
----------
[Mandic] Mandic, L. more_dynamics, ``plugins/force_producers/jet_nozzle.py`` 58-73
    (the thrust state, read, not edited).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive

DRIVE_OUTPUT_UNIT = "1"


def normalized_lag_parameters():
    return (positive("thrust_rise_time", "s", "time constant when the effort rises"),
            positive("thrust_fall_time", "s", "time constant when the effort falls"))


def normalized_lag_casadi():
    declared = normalized_lag_parameters()
    p = symbols(declared)
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 1), ca.SX.sym("shaft_torque")
    effort = ca.fmax(ca.fmin(state[0], 1.0), -1.0)  # (jet_nozzle.py 58)
    ref = ca.fmax(ca.fmin(command[0], 1.0), -1.0)  # (jet_nozzle.py 65)
    dot = ca.if_else(ref > effort, (ref - effort) / p["thrust_rise_time"], (ref - effort) / p["thrust_fall_time"])  # (jet_nozzle.py 73)
    return part_function("normalized_lag", {"command": command, "drive_state": state, "shaft_torque": load},
                         {"drive_output": effort, "drive_state_dot": dot, "thrust_scale": ca.SX(1.0),
                          "thrust_limit": ca.SX(ca.inf)}, p, declared)
