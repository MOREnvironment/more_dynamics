"""Thrust law of a unit commanded by its thrust: the thrust is the command limited to ``+-T_lim`` with
``T_lim = min(T_max, P / U)``, the maximum thrust and the power delivered to the water over the advance speed ``U``
(a constant-power unit gives less thrust the faster it goes); no torque, the water density not used.

    T = sat(command, -T_lim, T_lim),   T_lim = min(T_max, P / U)

The command is the drive's output in newtons (unit ``"N"``). The power cap is the vendor's figure of the power
delivered to the water divided by the speed, a first estimate of the thrust the motor can give at speed, not a
measured thrust curve; raise: bollard pulls and runs at known speed and throttle.

References
----------
[Krizman] Krizman, E. Model file ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)``, MSS fork branch ``grethe-vehicle`` @ 679a6e8,
    CRAFT/USV/models/grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 391-393 (``T_lim = min(T_max, P_prop/U_h)``, ``Th = satlim(ui(1), -T_lim, T_lim)``; read,
    the file may still be revised); its ``T_max`` and ``P_prop`` are the Torqeedo data sheet's maximum thrust and power
    delivered to the water (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 262-263, vendor document, not read here).

Author:    Enio Krizman
Date:      2026-10-26
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive
from more_dynamics.models.shared.force_producer_common import saturate

LAW_INPUT_UNIT = "N"


def thrust_power_limited_parameters():
    return (positive("max_thrust", "N", "maximum thrust of the unit T_max"),
            positive("delivered_power", "W", "power delivered to the water P"))


def thrust_power_limited_casadi():
    declared = thrust_power_limited_parameters()
    p = symbols(declared)
    command, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    limit = ca.fmin(p["max_thrust"], p["delivered_power"] / advance)  # T_lim = min(T_max, P / U) (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 392)
    return part_function("thrust_power_limited", {"drive_output": command, "advance_speed": advance,
                                                  "water_density": density},
                         {"thrust": saturate(command, -limit, limit), "torque": ca.SX(0.0)}, p, declared)  # (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 393)
