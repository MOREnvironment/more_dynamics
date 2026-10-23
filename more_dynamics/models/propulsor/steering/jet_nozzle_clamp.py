"""Nozzle angle of a steered jet as the shared repository's jet nozzle plugin has it, kept as
a labelled comparison form: a normalised command in [-1, 1] mapped to ``max_angle`` (ahead)
or ``min_angle`` (astern), the angle a state clamped to [min, max], its derivative the angle
error clamped to +-``nozzle_velocity``.

Below the clamp this is a first-order lag of 1 s (the error is in rad and used as rad/s),
above it a rate limit; the lag with rate and angle limits of ``servo.servo`` equals it
inside the limits when the limits are symmetric, the time constant is 1 s, the rate limit is
``nozzle_velocity`` and the command is given in radians.

References
----------
[Mandic] Mandic, L. more_dynamics, ``plugins/force_producers/jet_nozzle.py`` 59-79 (read,
    not edited).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive, scalar


def jet_nozzle_clamp_parameters():
    return (scalar("min_angle", "rad", "lowest nozzle angle"), scalar("max_angle", "rad", "highest nozzle angle"),
            positive("nozzle_velocity", "rad/s", "bound of the angle error that drives the nozzle"))


def jet_nozzle_clamp_casadi():
    declared = jet_nozzle_clamp_parameters()
    p = symbols(declared)
    command, state = ca.SX.sym("command", 1), ca.SX.sym("servo_state", 1)
    angle = ca.fmax(ca.fmin(state[0], p["max_angle"]), p["min_angle"])  # (jet_nozzle.py 62)
    c = ca.fmax(ca.fmin(command[0], 1.0), -1.0)  # (jet_nozzle.py 66)
    ref = ca.if_else(c >= 0.0, c * p["max_angle"], -c * p["min_angle"])  # (jet_nozzle.py 67-71)
    v = p["nozzle_velocity"]
    dot = ca.if_else(ref > angle, ca.fmin(ref - angle, v), ca.fmax(ref - angle, -v))  # (jet_nozzle.py 75-79)
    return part_function("jet_nozzle_clamp", {"command": command, "servo_state": state},
                         {"deflection": angle, "servo_state_dot": dot}, p, declared)
