"""Waterjet thrust from the momentum flux: ``T = rho Q (Q / A_n - Va)`` with the flow ``Q = c_Q n`` through
the nozzle of area ``A_n`` proportional to the shaft speed, for ``n > 0``; zero for a stopped or reversed
shaft; no torque.

The jet leaves at ``v_out = Q / A_n`` and the water enters at the advance speed ``v_in = Va`` (Bulten 2006,
eq. 2.52, p. 33). The thrust falls linearly with the intake speed, is zero when the intake runs at the jet
speed and turns into drag above it. With ``Q`` proportional to ``n`` the law is exactly the bilinear law with
``b1 = rho c_Q^2 / A_n``, ``b2 = -rho c_Q`` and the bollard law at ``Va = 0`` with ``k = rho c_Q^2 / A_n``.
Anchor: the ideal bollard identity ``T = (4 rho A_n P^2)^(1/3)`` of a lossless jet taking the shaft power
``P`` (Loghis and Xiros 2022, App. A, eqs. A2-A8, pp. 23-24).

``flow_coefficient_method``: ``"given"`` declares ``c_Q``; ``"from_bollard_point"`` declares one measured
bollard point (thrust ``T_b`` at shaft speed ``n_b``) and takes ``c_Q = sqrt(T_b A_n / rho) / n_b``, the
inverse of ``T_b = rho c_Q^2 n_b^2 / A_n`` (Kandasamy et al. 2010, eq. 3 as quoted in the literature
memo; with the fed density the point is reproduced at any density).

References
----------
[Bulten 2006] Bulten, N. W. H. (2006). Numerical analysis of a waterjet propulsion system. PhD thesis,
    TU Eindhoven. Ch. 2, eq. 2.52, p. 33.
[Loghis 2022] Loghis, A. and Xiros, N. (2022). Waterjet propulsion model, App. A, eqs. A2-A8, pp. 23-24.
[Kandasamy 2010] Kandasamy, M., et al. (2010). Waterjet thrust from flow, eq. 3.

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive

LAW_INPUT_UNIT = "rev/s"
FLOW_COEFFICIENT_METHODS = ("given", "from_bollard_point")


def _check(flow_coefficient_method):
    if flow_coefficient_method not in FLOW_COEFFICIENT_METHODS:
        raise ValueError(f"flow_coefficient_method must be one of {FLOW_COEFFICIENT_METHODS}, "
                         f"got {flow_coefficient_method!r}")


def jet_momentum_flux_parameters(*, flow_coefficient_method="given"):
    _check(flow_coefficient_method)
    if flow_coefficient_method == "given":
        return (positive("nozzle_area", "m^2", "nozzle outlet area A_n"),
                positive("flow_coefficient", "m^3/rev", "c_Q of Q = c_Q n"))
    return (positive("nozzle_area", "m^2", "nozzle outlet area A_n"),
            positive("bollard_thrust", "N", "thrust T_b measured at the bollard"),
            positive("bollard_shaft_speed", "rev/s", "shaft speed n_b of that bollard point"))


def jet_momentum_flux_casadi(*, flow_coefficient_method="given"):
    _check(flow_coefficient_method)
    declared = jet_momentum_flux_parameters(flow_coefficient_method=flow_coefficient_method)
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    if flow_coefficient_method == "given":
        flow_coefficient = p["flow_coefficient"]
    else:
        flow_coefficient = ca.sqrt(p["bollard_thrust"] * p["nozzle_area"] / density) / p["bollard_shaft_speed"]  # inverse of T_b = rho c_Q^2 n_b^2 / A_n
    flow = flow_coefficient * n  # Q = c_Q n
    thrust = ca.if_else(n > 0.0, density * flow * (flow / p["nozzle_area"] - advance), 0.0)  # T = rho Q (v_out - v_in) (Bulten 2006, eq. 2.52, p. 33)
    return part_function("jet_momentum_flux", {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": thrust, "torque": ca.SX(0.0)}, p, declared)
