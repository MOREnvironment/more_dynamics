"""Waterjet thrust and shaft torque from the match of the pump to the system: the flow ``Q`` is the
intersection of a linear pump head line with the system head, ``T = rho Q (v_out - Va)`` and the shaft
torque follows from the pump head; for ``n > 0``, zero otherwise.

Pump line: ``psi = psi_0 - s phi`` with ``phi = Q / (Omega D^3)``, ``psi = g H / (Omega D)^2`` and
``Omega = 2 pi n`` (Bulten 2006, eqs. 2.15-2.16, p. 22), so ``g H_pump = psi_0 (Omega D)^2 - s Omega Q / D``.
System head: ``g H_R = v_out^2 (1 + phi_n) / 2 - Va^2 (1 - eps) / 2`` with ``v_out = Q / A_n``, the nozzle loss
``phi_n``, the inlet loss ``eps`` and no static head (eq. 2.56, p. 35). Setting ``H_pump = H_R`` gives
``a Q^2 + b Q - c = 0`` with ``a = (1 + phi_n) / (2 A_n^2)``, ``b = s Omega / D`` and ``c = psi_0 (Omega D)^2 +
Va^2 (1 - eps) / 2``, and the positive root ``Q = 2 c / (b + sqrt(b^2 + 4 a c))`` (closed form of the
intersection). Shaft torque ``Q_s = rho g H_pump Q / (eta_p Omega)`` (eq. 2.20, p. 22).

References
----------
[Bulten 2006] Bulten, N. W. H. (2006). Numerical analysis of a waterjet propulsion system. PhD thesis,
    TU Eindhoven. Ch. 2, eqs. 2.15, 2.16, 2.20, 2.56, pp. 22, 35.

Author:    Enio Krizman
Date:      2026-10-23
"""
import math

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from more_dynamics.models.propulsor.declarations import fraction, part_function, positive, scalar

LAW_INPUT_UNIT = "rev/s"


def jet_pump_system_match_parameters():
    return (positive("impeller_diameter", "m", "impeller diameter D"),
            positive("nozzle_area", "m^2", "nozzle outlet area A_n"),
            positive("pump_head_coefficient_at_zero_flow", "1", "psi_0 of the pump line psi = psi_0 - s phi"),
            Parameter("pump_head_slope", (1, 1), "1", "s of the pump line psi = psi_0 - s phi", 0.0),
            Parameter("nozzle_loss_coefficient", (1, 1), "1", "phi_n: nozzle loss as a fraction of the jet head", 0.0),
            fraction("inlet_loss_coefficient", "eps: share of the intake ram head lost"),
            Parameter("pump_efficiency", (1, 1), "1", "eta_p of the shaft torque", 0.0, 1.0, minimum_exclusive=True))


def jet_pump_system_match_casadi():
    declared = jet_pump_system_match_parameters()
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    running = n > 0.0
    omega = 2.0 * math.pi * ca.if_else(running, n, 1.0)  # shaft speed in rad/s (a guard value when not running)
    d, area, psi0, slope = p["impeller_diameter"], p["nozzle_area"], p["pump_head_coefficient_at_zero_flow"], p["pump_head_slope"]
    a = (1.0 + p["nozzle_loss_coefficient"]) / (2.0 * area**2)
    b = slope * omega / d
    c = psi0 * (omega * d) ** 2 + advance**2 * (1.0 - p["inlet_loss_coefficient"]) / 2.0
    flow = 2.0 * c / (b + ca.sqrt(b**2 + 4.0 * a * c))  # positive root of a Q^2 + b Q - c = 0
    head = psi0 * (omega * d) ** 2 - slope * omega * flow / d  # g H of the pump line (Bulten 2006, eqs. 2.15-2.16, p. 22)
    thrust = density * flow * (flow / area - advance)  # (Bulten 2006, eq. 2.52, p. 33)
    torque = density * head * flow / (p["pump_efficiency"] * omega)  # (Bulten 2006, eq. 2.20, p. 22)
    return part_function("jet_pump_system_match", {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": ca.if_else(running, thrust, 0.0), "torque": ca.if_else(running, torque, 0.0)},
                         p, declared)
