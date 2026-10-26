"""Inflow of a shaft that sees the surge component of the relative velocity: ``Va = c_Va (1 - w) u_r``
(signed: astern flow gives a negative advance speed), thrust factor ``1 - t``, torque factor 1.

``c_Va`` and the use of the surge component alone are constructions of the steerable outboard
block this form comes from; ``1 - t`` is the thrust deduction of MSS.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 250
    @ ed2ce7b (thrust ``(1 - t) X``).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import fraction, part_function, scalar


def surge_component_parameters():
    return (fraction("thrust_deduction", "thrust deduction t: thrust (1 - t) T"),
            scalar("advance_speed_factor", "1", "c_Va of Va = c_Va (1 - w) u_r"),
            fraction("wake_fraction", "w of Va = c_Va (1 - w) u_r"))


def surge_component_casadi():
    declared = surge_component_parameters()
    p = symbols(declared)
    inputs = {"nu_r": ca.SX.sym("nu_r", 6), "thrust_axis": ca.SX.sym("thrust_axis", 3),
              "thrust_position": ca.SX.sym("thrust_position", 3)}
    advance = p["advance_speed_factor"] * (1.0 - p["wake_fraction"]) * inputs["nu_r"][0]  # Va = c_Va (1 - w) u_r: construction
    return part_function("surge_component", inputs, {"advance_speed": advance, "thrust_factor": 1.0 - p["thrust_deduction"],
                                                     "torque_factor": ca.SX(1.0)}, p, declared)  # (1 - t) (remus100.m 250)
