"""Two fixed propellers with a quadratic thrust law (differential thrust), CasADi.

``differential_thruster_parameters()`` declares the numbers (name, shape, SI
unit, meaning, admissible range); ``differential_thruster_casadi()`` builds
``(n, nu_r, <parameters by name>) -> (tau, allocation_matrix, max_speed,
min_speed)``. Per-thruster quantities are ordered ``[left, right]``.

Equations (keys in References):

* ``T_i = k_pos n_i |n_i|`` for ``n_i > 0``, else ``k_neg n_i |n_i|``
  (Fossen 2011, eq. 12.265, p. 411, with separate forward and reverse
  coefficients; MSS ``otter.m`` 223-229).
* Limits ``n_max = sqrt(F_fwd / k_pos)``, ``n_min = -sqrt(F_rev / k_neg)``
  (MSS ``otter.m`` 136-137, there with ``F = 0.5 * 24.4 g`` and
  ``0.5 * 13.6 g``), ``n`` saturated to ``[n_min, n_max]`` (``otter.m`` 220).
* ``tau = B T``, column ``i`` of ``B`` = ``[d_i; r_i x d_i]`` for the unit
  direction ``d_i = e_i / |e_i|`` at ``r_i`` (Fossen 2011, eqs. 12.226,
  p. 400, and 12.229, p. 401; MSS ``otter.m`` 232, which writes it for two
  thrusters along x at ``y = -+ y_pont``, ``otter.m`` 132-133).

The shaft-speed lag is a state of the assembly, not part of this map; ``n``
is the actual shaft speed. A zero direction is refused by
``check_differential_thruster_values``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 2,
    eq. 2.9, p. 20; Ch. 12, eqs. 12.226-12.229, pp. 400-401; eq. 12.265,
    p. 411.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/USV/models/otter.m`` 132-137,
    220-232.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter, check_values, symbols

from more_dynamics.models.shared.force_producer_common import saturate


def _positive_pair(name, unit, meaning):
    return Parameter(name, (2, 1), unit, meaning, 0.0, minimum_exclusive=True)


_DECLARED = (
    _positive_pair("positive_thrust_coefficients", "N*s^2/rad^2", "k_pos of [left, right]: T = k_pos n|n| for n > 0"),
    _positive_pair("negative_thrust_coefficients", "N*s^2/rad^2", "k_neg of [left, right]: T = k_neg n|n| for n <= 0"),
    Parameter("thruster_positions", (2, 3), "m", "rows: CO -> thruster [left; right], body axes (FRD)"),
    Parameter("thruster_directions", (2, 3), "1", "rows: thrust direction [left; right], normalised inside"),
    _positive_pair("max_forward_thrust", "N", "forward thrust at n_max of [left, right]"),
    _positive_pair("max_reverse_thrust", "N", "reverse thrust magnitude at n_min of [left, right]"),
)

DIFFERENTIAL_THRUSTER_OUTPUTS = ("tau", "allocation_matrix", "max_speed", "min_speed")


def differential_thruster_parameters():
    """The declared parameter set: a tuple of ``Parameter`` (name, shape, SI
    unit, meaning, admissible range), in the order of the block's inputs."""
    return _DECLARED


def differential_thruster_casadi():
    """The block: ``(n, nu_r, <parameters by name>) -> (tau, allocation_matrix,
    max_speed, min_speed)``.

    Contract
    --------
    Inputs: ``n`` (2x1, actual shaft speeds [left, right] in rad/s, saturated
    inside), ``nu_r`` (6x1, not used by this law, kept for the common
    producer signature) and the parameters of
    ``differential_thruster_parameters()`` by name. Outputs: ``tau`` (6x1,
    N and N m, BODY, about the CO), ``allocation_matrix`` (6x2, ``B``),
    ``max_speed`` and ``min_speed`` (2x1, rad/s). No value is checked inside
    the graph (``check_differential_thruster_values``).
    """
    p = symbols(_DECLARED)
    n = ca.SX.sym("n", 2)
    nu_r = ca.SX.sym("nu_r", 6)
    k_pos, k_neg = p["positive_thrust_coefficients"], p["negative_thrust_coefficients"]

    columns = []
    for i in range(2):
        direction = p["thruster_directions"][i, :].T
        unit = direction / ca.norm_2(direction)  # d_i = e_i / |e_i| (unit direction; Fossen 2011, eq. 12.226, p. 400)
        arm = MatrixTransforms.skew(p["thruster_positions"][i, :].T) @ unit  # r_i x d_i = S(r_i) d_i (Fossen 2011, eq. 12.226, p. 400; eq. 2.9, p. 20)
        columns.append(ca.vertcat(unit, arm))
    allocation = ca.horzcat(*columns)  # B, columns [d_i; r_i x d_i] (Fossen 2011, eq. 12.229, p. 401; otter.m 232)

    max_speed = ca.sqrt(p["max_forward_thrust"] / k_pos)  # (otter.m 136)
    min_speed = -ca.sqrt(p["max_reverse_thrust"] / k_neg)  # (otter.m 137)
    n_sat = saturate(n, min_speed, max_speed)  # (otter.m 220)
    quadratic = n_sat * ca.fabs(n_sat)  # n |n| (Fossen 2011, eq. 12.265, p. 411)
    # (Fossen 2011, eq. 12.265, p. 411; otter.m 223-229)
    thrust = ca.if_else(n_sat > 0.0, k_pos * quadratic, k_neg * quadratic)
    tau = allocation @ thrust  # (Fossen 2011, eq. 12.229, p. 401; otter.m 232)

    names = [d.name for d in _DECLARED]
    outputs = {"tau": tau, "allocation_matrix": allocation, "max_speed": max_speed, "min_speed": min_speed}
    return ca.Function(
        "differential_thruster",
        [n, nu_r, *[p[name] for name in names]],
        [outputs[name] for name in DIFFERENTIAL_THRUSTER_OUTPUTS],
        ["n", "nu_r", *names],
        list(DIFFERENTIAL_THRUSTER_OUTPUTS),
    )


def check_differential_thruster_values(values):
    """Numbers checked: ``{name: ca.DM}`` (``check_values``: names, shapes,
    finite, positive coefficients and thrusts); a zero thrust direction is
    refused, naming it."""
    numbers = check_values(_DECLARED, values)
    directions = numbers["thruster_directions"]
    for i in range(2):
        if not float(ca.norm_2(directions[i, :])) > 0.0:
            raise ValueError(f"parameter 'thruster_directions' [1]: row {i} must be non-zero")
    return numbers
