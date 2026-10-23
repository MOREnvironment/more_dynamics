"""Mounting of a steered jet by rotation: the axis is ``R_z(a + b) R_y(t) e_x`` for steering
angle ``a``, yaw bias ``b`` and trim ``t``, a unit vector at every angle:
``[cos(a + b) cos t, sin(a + b) cos t, -sin t]``.

Equal to the direction of ``jet_nozzle_direction`` at zero trim; with trim the force keeps the
commanded length where that form lengthens it.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 2, eqs. 2.15-2.18, p. 22 (elementary rotations).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms, symbols

from more_dynamics.models.propulsor.declarations import part_function, scalar, vector


def steered_rotation_parameters():
    return (vector("position", "m", "CO to the nozzle, BODY (FRD)"), scalar("trim_angle", "rad", "nozzle trim"),
            scalar("yaw_bias", "rad", "yaw bias of the nozzle"))


def steered_rotation_casadi():
    declared = steered_rotation_parameters()
    p = symbols(declared)
    angle = ca.SX.sym("steering_angle")
    axis = (MatrixTransforms.Rz(angle + p["yaw_bias"]) @ MatrixTransforms.Ry(p["trim_angle"]))[:, 0]  # R_z R_y e_x (Fossen 2011, eq. 2.15, p. 22)
    return part_function("steered_rotation", {"steering_angle": angle},
                         {"thrust_axis": axis, "thrust_position": p["position"]}, p, declared)
