"""Mounting of a steered jet with the direction formula of the shared repository's jet nozzle
plugin, kept as a labelled comparison form: for steering angle ``a``, yaw bias ``b`` and trim
``t`` the axis is ``[cos(a + b) cos t, sin(a + b), -sin t]``.

The axis is **not** a unit vector when the trim is not zero and the nozzle is not straight:
``|axis|^2 = 1 + sin^2 t sin^2(a + b)``, so the thrust vector is longer than the commanded
thrust (up to 0.54 % at 0.2 rad of trim). ``steered_rotation`` is the unit-length form.

References
----------
[Mandic] Mandic, L. more_dynamics, ``plugins/force_producers/jet_nozzle.py`` 91-95 (read,
    not edited).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, scalar, vector


def jet_nozzle_direction_parameters():
    return (vector("position", "m", "CO to the nozzle, BODY (FRD)"), scalar("trim_angle", "rad", "nozzle trim"),
            scalar("yaw_bias", "rad", "yaw bias of the nozzle"))


def jet_nozzle_direction_casadi():
    declared = jet_nozzle_direction_parameters()
    p = symbols(declared)
    angle = ca.SX.sym("steering_angle")
    heading = angle + p["yaw_bias"]
    axis = ca.vertcat(ca.cos(heading) * ca.cos(p["trim_angle"]), ca.sin(heading), -ca.sin(p["trim_angle"]))  # (jet_nozzle.py 91-95)
    return part_function("jet_nozzle_direction", {"steering_angle": angle},
                         {"thrust_axis": axis, "thrust_position": p["position"]}, p, declared)
