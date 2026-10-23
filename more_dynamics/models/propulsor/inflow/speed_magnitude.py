"""Inflow of a shaft that sees the vehicle's speed through the water: ``Va = (1 - w) |nu_r[0:3]|``,
thrust factor ``1 - t``, torque factor 1. The speed has a zero (not NaN) derivative at standstill.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 151,
    250 @ 99fdf4d (``Va = (1 - w) U_r``, thrust ``(1 - t) X``).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import fraction, part_function
from more_dynamics.models.shared.force_producer_common import safe_norm


def speed_magnitude_parameters():
    return (fraction("thrust_deduction", "thrust deduction t: thrust (1 - t) T"),
            fraction("wake_fraction", "wake fraction w: Va = (1 - w) U_r"))


def speed_magnitude_casadi():
    declared = speed_magnitude_parameters()
    p = symbols(declared)
    inputs = {"nu_r": ca.SX.sym("nu_r", 6), "thrust_axis": ca.SX.sym("thrust_axis", 3),
              "thrust_position": ca.SX.sym("thrust_position", 3)}
    advance = (1.0 - p["wake_fraction"]) * safe_norm(inputs["nu_r"][0:3])  # (remus100.m 151)
    return part_function("speed_magnitude", inputs, {"advance_speed": advance, "thrust_factor": 1.0 - p["thrust_deduction"],
                                                     "torque_factor": ca.SX(1.0)}, p, declared)  # (1 - t) (remus100.m 250)
