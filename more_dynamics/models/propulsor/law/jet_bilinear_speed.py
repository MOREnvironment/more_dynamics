"""Waterjet thrust bilinear in the shaft speed and the advance speed: ``T = b1 n^2 + b2 n Va`` for ``n > 0``,
zero for a stopped or reversed shaft; no torque.

Anchor: the thrust ratio ``1 + (b2 / b1) u`` of an identified waterjet craft at full throttle (Han et al.
2017, eq. 7 and Table III).

References
----------
[Han 2017] Han, J., et al. (2017). Waterjet craft model identification, eq. 7 and Table III.

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive, scalar

LAW_INPUT_UNIT = "rev/s"


def jet_bilinear_speed_parameters():
    return (positive("speed_squared_coefficient", "N*s^2/rev^2", "b1 of T = b1 n^2 + b2 n Va"),
            scalar("speed_advance_coefficient", "N*s^2/(rev*m)", "b2 of T = b1 n^2 + b2 n Va"))


def jet_bilinear_speed_casadi():
    declared = jet_bilinear_speed_parameters()
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    thrust = ca.if_else(n > 0.0, p["speed_squared_coefficient"] * n**2 + p["speed_advance_coefficient"] * n * advance, 0.0)  # (Han 2017, eq. 7)
    return part_function("jet_bilinear_speed", {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": thrust, "torque": ca.SX(0.0)}, p, declared)
