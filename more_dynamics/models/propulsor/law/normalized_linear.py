"""Thrust law on a normalised effort: ``T = T_max c u`` with ``u`` the effort (unit 1), no torque,
indifferent to the advance speed and the water density.

Kept as a labelled comparison form: it is the thrust law of the shared repository's jet nozzle
plugin, which carries ``T_max`` and ``c`` as two numbers of which only the product acts.

References
----------
[Mandic] Mandic, L. more_dynamics, ``plugins/force_producers/jet_nozzle.py`` 89-91 (read, not
    edited).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive, scalar

LAW_INPUT_UNIT = "1"


def normalized_linear_parameters():
    return (positive("max_thrust", "N", "thrust at an effort of 1"), scalar("thrust_coefficient", "1", "scale of the thrust"))


def normalized_linear_casadi():
    declared = normalized_linear_parameters()
    p = symbols(declared)
    effort, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    return part_function("normalized_linear", {"drive_output": effort, "advance_speed": advance, "water_density": density},
                         {"thrust": p["max_thrust"] * p["thrust_coefficient"] * effort, "torque": ca.SX(0.0)}, p, declared)
