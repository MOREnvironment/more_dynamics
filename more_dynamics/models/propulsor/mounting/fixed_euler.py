"""Mounting of a fixed shaft given by Euler angles: the axis is the first column of
``R_zyx(orientation)`` in BODY, the point the declared position; the steering angle is ignored.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 2, eq. 2.18, p. 22.
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 250
    @ 6e78e98 (thrust on the first body axis).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms, symbols

from more_dynamics.models.propulsor.declarations import part_function, vector


def fixed_euler_parameters():
    return (vector("position", "m", "CO to the thrust point, BODY (FRD)"),
            vector("orientation", "rad", "shaft orientation (roll, pitch, yaw), zyx Euler angles"))


def fixed_euler_casadi():
    declared = fixed_euler_parameters()
    p = symbols(declared)
    angle = ca.SX.sym("steering_angle")
    axis = MatrixTransforms.Rzyx_explicit(p["orientation"])[:, 0]  # (Fossen 2011, eq. 2.18, p. 22)
    return part_function("fixed_euler", {"steering_angle": angle}, {"thrust_axis": axis, "thrust_position": p["position"]},
                         p, declared)
