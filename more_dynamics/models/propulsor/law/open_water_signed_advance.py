"""Open-water law with a signed advance number: ``J = Va / (n D)`` (signed), ``K_T`` and ``K_Q`` of the
Wageningen polynomial at that ``J`` (no clip), the reverse factors on both for ``n < 0``, zero force when
the shaft is stopped; ``T = rho D^4 K_T n |n|``, ``Q = rho D^5 K_Q n |n|``.

The reverse factors, the stop guard and the signed ``J`` are constructions of the steerable outboard block
this form comes from; the polynomial and the quadratic law are published. The function is
``more_dynamics.models.outboard.outboard_motor.signed_advance_open_water``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control.
    Wiley. Ch. 12, eq. 12.265, p. 411.
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, LIBRARY/modeling/wageningen.m 6, 35-38
    @ 6e78e98.

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.outboard.outboard_motor import signed_advance_open_water
from more_dynamics.models.propeller.wageningen_kt_kq import wageningen_parameters
from more_dynamics.models.propulsor.declarations import part_function, positive

LAW_INPUT_UNIT = "rev/s"


def open_water_signed_advance_parameters():
    return (positive("diameter", "m", "propeller diameter D"), *wageningen_parameters(),
            positive("reverse_thrust_factor", "1", "factor on K_T for n < 0"),
            positive("reverse_torque_factor", "1", "factor on K_Q for n < 0"))


def open_water_signed_advance_casadi():
    declared = open_water_signed_advance_parameters()
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    thrust, torque = signed_advance_open_water(density, p["diameter"], p["pitch_diameter_ratio"], p["blade_area_ratio"],
                                               p["blade_count"], p["reverse_thrust_factor"], p["reverse_torque_factor"],
                                               n, advance)
    return part_function("open_water_signed_advance", {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": thrust, "torque": torque}, p, declared)
