"""One DUNE VSIM fin as a force-producer leaf.

Its simulator moment law is distinct from a point-force cross product.
The fed water density is unused: max_force already contains the gain.

References
----------
[DUNE] LSTS. DUNE, EUPL v1.1, src/Simulators/VSIM/VSIM/Fin.cpp 56-102 @ 555ef0b.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols
from more_dynamics.models.shared.force_producer_common import saturate


def proportional_force_fin_parameters():
    return (Parameter("max_force", (3, 1), "N", "simulator force gain along BODY axes"),
            Parameter("fin_position", (3, 1), "m", "CO to fin in BODY axes"),
            Parameter("max_deflection", (1, 1), "rad", "deflection angle limit", 0.0, minimum_exclusive=True))


def proportional_force_fin_casadi():
    declared = proportional_force_fin_parameters()
    p = symbols(declared)
    command, state = ca.SX.sym("command"), ca.SX.sym("state", 0)
    nu_r, density = ca.SX.sym("nu_r", 6), ca.SX.sym("water_density")
    normalized = saturate(command, -p["max_deflection"], p["max_deflection"]) / p["max_deflection"]  # (DUNE Fin.cpp 59-66)
    speed_squared = ca.sumsqr(nu_r[:3])  # (DUNE Fin.cpp 86-93)
    force = speed_squared * normalized * p["max_force"] / 2  # (DUNE Fin.cpp 56-67, 92-93)
    r, f = p["fin_position"], p["max_force"]
    moment = speed_squared * normalized * ca.vertcat(
        ca.sqrt(r[1]**2 + r[2]**2) * ca.fabs(f[0]) / 2,
        r[0] * ca.fabs(f[2]) / 2,
        r[0] * ca.fabs(f[1]) / 2)  # (DUNE Fin.cpp 69-83, 95-101; K,M,N ordering)
    tau = ca.vertcat(force, moment)
    tau = ca.if_else(ca.sqrt(speed_squared) <= 1e-6, ca.DM.zeros(6), tau)  # same zero-speed guard as vsim_fins.py
    names = [d.name for d in declared]
    return ca.Function("proportional_force_fin", [command, state, nu_r, density, *[p[n] for n in names]],
                       [tau, ca.SX.zeros(0, 1)], ["command", "state", "nu_r", "water_density", *names],
                       ["tau", "state_dot"])
