"""Waterjet thrust at the bollard, quadratic in the shaft speed: ``T = k n^2`` for ``n > 0``, zero for a
stopped or reversed shaft (a waterjet needs a reversing bucket for astern thrust; none is modelled), no
torque, no effect of the boat speed.

Anchor: the published bollard law ``F = 1.89e-4 N^2`` with ``N`` in rpm, i.e. ``k = 1.89e-4 * 3600``
N s^2/rev^2 (Tajima et al. 2015, eq. 7 and Fig. 3, p. 378).

References
----------
[Tajima 2015] Tajima, Y., et al. (2015). Bollard thrust law of a water-jet propulsor, eq. 7 and Fig. 3,
    p. 378.
[Varyani 2009] Varyani, K. S. and Krishnankutty, P. (2009). J. Mar. Eng. Technol. 13, p. 14 (a waterjet
    gives astern thrust only through a reversing bucket).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive

LAW_INPUT_UNIT = "rev/s"


def jet_bollard_quadratic_parameters():
    return (positive("bollard_thrust_coefficient", "N*s^2/rev^2", "k of T = k n^2 at the bollard"),)


def jet_bollard_quadratic_casadi():
    declared = jet_bollard_quadratic_parameters()
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    thrust = ca.if_else(n > 0.0, p["bollard_thrust_coefficient"] * n**2, 0.0)  # (Tajima 2015, eq. 7, p. 378)
    return part_function("jet_bollard_quadratic", {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": thrust, "torque": ca.SX(0.0)}, p, declared)
