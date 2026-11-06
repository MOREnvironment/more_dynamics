"""Power law of the shaft speed: the net thrust is a power of the shaft speed, ``T = k |n|^p`` ahead and ``T = -r_a k |n|^p`` astern, the
reaction torque a fixed fraction of the thrust. The law has no advance-speed term.

Frames and units: ``n`` the shaft speed in rev/s (positive ahead, the unit the shaft-speed drives pass), ``T`` the net thrust along the
thrust axis in N (after the thrust deduction, so ``k`` contains it and the open-water coefficient at the operating advance ratio),
``Q`` the torque about the thrust axis in N m that the propulsor passes to the hull (the propulsor applies ``-Q`` as a couple about the
axis, scaled by its ``shaft_torque_scale``). ``k`` has the unit N/(rev/s)^p.

Equations::

    T(n) =  k |n|^p                  n >= 0                      (the quadratic law F = K n |n| for p = 2, Fossen 2011, eq. 12.265, p. 411)
    T(n) = -r_a k |n|^p              n <  0                      (astern ratio r_a: the astern thrust at -n as a fraction of the ahead thrust at n)
    Q    =  c T                      (c the torque about the thrust axis per unit of net thrust that reaches the hull: positive for a propeller that
                                      turns in the positive sense about the axis when it drives ahead, negative for the opposite hand)

``p = 2`` is the quadratic law of a fixed-pitch screw, ``p = 1`` the linear one; the book states that "some propellers show linear
behavior in n while others are quadratic, even combinations are observed" (Fossen 2011, p. 411), a combined law being a power between
the two or above. The exponent is a number of the identification of a vehicle together with the resistance of its hull: along a steady
run whose resistance goes as the speed squared the advance ratio is constant and the speed dependence of the open-water coefficient
cannot be told from the resistance, so ``k`` and ``p`` are fitted together with the hull's drag law and are not properties of the screw
alone. A propeller with a known open-water curve takes an open-water law. The advance speed and the water density are inputs of the
part and are not used.

Numerics: ``n = 0`` returns exactly zero thrust and torque and the derivative with respect to ``n`` there is taken as zero (it is zero for
``p > 1``; the law is not differentiable at ``n = 0`` for ``p <= 1``); the power is evaluated on ``|n|`` replaced by 1 where ``n = 0``, so
no expression or derivative of the part is ``0^p`` or ``log 0``.

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Ch. 12, eq. 12.265, p. 411 (the quadratic
    law and the remark on linear and combined behaviour).

Author:    Enio Krizman
Date:      2026-11-06
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from more_dynamics.models.propulsor.declarations import part_function, positive, scalar

LAW_INPUT_UNIT = "rev/s"


def power_law_parameters():
    return (positive("power_law_coefficient", "N/(rev/s)^p", "k: T = k |n|^p at a shaft speed n of 1 rev/s ahead, after the thrust deduction"),
            positive("power_law_exponent", "1", "p: the exponent of the shaft speed (2: the quadratic law, 1: linear)"),
            Parameter("astern_thrust_ratio", (1, 1), "1", "r_a: the astern thrust at a shaft speed -n as a fraction of the ahead thrust at n", 0.0),
            scalar("reaction_torque_per_thrust", "N*m/N", "c: the torque Q = c T about the thrust axis that reaches the hull per unit of net "
                                                          "thrust (negative for the opposite hand of the propeller)"))


def power_law_casadi():
    declared = power_law_parameters()
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    moving = ca.fabs(n) > 0.0
    magnitude = ca.if_else(moving, p["power_law_coefficient"] * ca.power(ca.if_else(moving, ca.fabs(n), 1.0), p["power_law_exponent"]), 0.0)
    thrust = ca.if_else(n >= 0.0, magnitude, -p["astern_thrust_ratio"] * magnitude)  # k |n|^p ahead, -r_a k |n|^p astern
    torque = p["reaction_torque_per_thrust"] * thrust  # Q = c T
    return part_function("power_law", {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": thrust, "torque": torque}, p, declared)
