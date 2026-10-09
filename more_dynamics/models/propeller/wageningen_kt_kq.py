"""Wageningen B-series open-water coefficients ``K_T(J)``, ``K_Q(J)`` (CasADi).

``wageningen_parameters()`` declares the propeller geometry (name, shape, SI
unit, meaning, admissible range); ``wageningen_casadi(clip_advance_ratio=...)``
builds ``(J, pitch_diameter_ratio, blade_area_ratio, blade_count) -> (KT,
KQ)``; ``wageningen_expression`` gives the same two expressions for a
symbolic ``J`` and geometry, for the propeller and outboard blocks.

Equation (keys in References): ``K = sum_k C_k J^s_k (P/D)^t_k (A_E/A_O)^u_k
z^v_k``, rows ``[C, s, t, u, v]`` (MSS ``wageningen.m`` 33-38, which
attributes the regression to Barnitsas, Ray and Kinley (1981), not read
here). The tables below equal MSS ``LIBRARY/modeling/utiles/WageningData.txt``
row by row (thrust lines 6-44, torque lines 52-98; the same numbers as
``WageningData.mat``, which ``wageningen.m`` 33 loads). The geometry factor
of each term is summed per power of ``J`` first (``s`` = 0..3), so the
result is a cubic in ``J`` whose coefficients are expressions in the
geometry (the same sum, regrouped).

Default: MSS ``wageningen.m``, the polynomial at any ``J``, including
``J < 0`` and ``J > 1.3`` (outside the fitted range). Deviation, behind
``clip_advance_ratio=True`` (a selector): ``J`` replaced by ``min(|J|,
1.3)``, the setting of the earlier template reference (a construction, no
published source).

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``LIBRARY/modeling/wageningen.m`` 23-38;
    ``LIBRARY/modeling/utiles/WageningData.txt`` 6-44, 52-98
    (``WageningData.mat``).

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

MAX_ADVANCE_NUMBER = 1.3  # J bound of the clip selector (the regression's fitted range ends here)

# [C, s, t, u, v] of K_T (WageningData.txt 6-44)
_THRUST_TABLE = (
    (0.00880496, 0, 0, 0, 0),
    (-0.204554, 1, 0, 0, 0),
    (0.166351, 0, 1, 0, 0),
    (0.158114, 0, 2, 0, 0),
    (-0.147581, 2, 0, 1, 0),
    (-0.481497, 1, 1, 1, 0),
    (0.415437, 0, 2, 1, 0),
    (0.0144043, 0, 0, 0, 1),
    (-0.0530054, 2, 0, 0, 1),
    (0.0143481, 0, 1, 0, 1),
    (0.0606826, 1, 1, 0, 1),
    (-0.0125894, 0, 0, 1, 1),
    (0.0109689, 1, 0, 1, 1),
    (-0.133698, 0, 3, 0, 0),
    (0.00638407, 0, 6, 0, 0),
    (-0.00132718, 2, 6, 0, 0),
    (0.168496, 3, 0, 1, 0),
    (-0.0507214, 0, 0, 2, 0),
    (0.0854559, 2, 0, 2, 0),
    (-0.0504475, 3, 0, 2, 0),
    (0.010465, 1, 6, 2, 0),
    (-0.00648272, 2, 6, 2, 0),
    (-0.00841728, 0, 3, 0, 1),
    (0.0168424, 1, 3, 0, 1),
    (-0.00102296, 3, 3, 0, 1),
    (-0.0317791, 0, 3, 1, 1),
    (0.018604, 1, 0, 2, 1),
    (-0.00410798, 0, 2, 2, 1),
    (-0.000606848, 0, 0, 0, 2),
    (-0.0049819, 1, 0, 0, 2),
    (0.0025983, 2, 0, 0, 2),
    (-0.000560528, 3, 0, 0, 2),
    (-0.00163652, 1, 2, 0, 2),
    (-0.000328787, 1, 6, 0, 2),
    (0.000116502, 2, 6, 0, 2),
    (0.000690904, 0, 0, 1, 2),
    (0.00421749, 0, 3, 1, 2),
    (0.0000565229, 3, 6, 1, 2),
    (-0.00146564, 0, 3, 2, 2),
)

# [C, s, t, u, v] of K_Q (WageningData.txt 52-98)
_TORQUE_TABLE = (
    (0.00379368, 0, 0, 0, 0),
    (0.00886523, 2, 0, 0, 0),
    (-0.032241, 1, 1, 0, 0),
    (0.00344778, 0, 2, 0, 0),
    (-0.0408811, 0, 1, 1, 0),
    (-0.108009, 1, 1, 1, 0),
    (-0.0885381, 2, 1, 1, 0),
    (0.188561, 0, 2, 1, 0),
    (-0.00370871, 1, 0, 0, 1),
    (0.00513696, 0, 1, 0, 1),
    (0.0209449, 1, 1, 0, 1),
    (0.00474319, 2, 1, 0, 1),
    (-0.00723408, 2, 0, 1, 1),
    (0.00438388, 1, 1, 1, 1),
    (-0.0269403, 0, 2, 1, 1),
    (0.0558082, 3, 0, 1, 0),
    (0.0161886, 0, 3, 1, 0),
    (0.00318086, 1, 3, 1, 0),
    (0.015896, 0, 0, 2, 0),
    (0.0471729, 1, 0, 2, 0),
    (0.0196283, 3, 0, 2, 0),
    (-0.0502782, 0, 1, 2, 0),
    (-0.030055, 3, 1, 2, 0),
    (0.0417122, 2, 2, 2, 0),
    (-0.0397722, 0, 3, 2, 0),
    (-0.00350024, 0, 6, 2, 0),
    (-0.0106854, 3, 0, 0, 1),
    (0.00110903, 3, 3, 0, 1),
    (-0.000313912, 0, 6, 0, 1),
    (0.0035985, 3, 0, 1, 1),
    (-0.00142121, 0, 6, 1, 1),
    (-0.00383637, 1, 0, 2, 1),
    (0.0126803, 0, 2, 2, 1),
    (-0.00318278, 2, 3, 2, 1),
    (0.00334268, 0, 6, 2, 1),
    (-0.00183491, 1, 1, 0, 2),
    (0.000112451, 3, 2, 0, 2),
    (-0.0000297228, 3, 6, 0, 2),
    (0.000269551, 1, 0, 1, 2),
    (0.00083265, 2, 0, 1, 2),
    (0.00155334, 0, 2, 1, 2),
    (0.000302683, 0, 6, 1, 2),
    (-0.0001843, 0, 0, 2, 2),
    (-0.000425399, 0, 3, 2, 2),
    (0.0000869243, 3, 3, 2, 2),
    (-0.0004659, 0, 6, 2, 2),
    (0.0000554194, 1, 6, 2, 2),
)


def _positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


_DECLARED = (
    _positive("pitch_diameter_ratio", "1", "pitch ratio P/D"),
    _positive("blade_area_ratio", "1", "expanded blade-area ratio A_E/A_O"),
    _positive("blade_count", "1", "number of blades z"),
)


def wageningen_parameters():
    """The declared parameter set: a tuple of ``Parameter`` (name, shape, SI
    unit, meaning, admissible range), in the order of the block's inputs."""
    return _DECLARED


def _collapse(table, pd, aeao, z):
    """``[c_0, c_1, c_2, c_3]`` with ``K = sum_s c_s J^s``: the geometry
    factor of each row summed into its power of ``J`` (wageningen.m 35-38,
    regrouped)."""
    polynomial = [0.0, 0.0, 0.0, 0.0]
    for c, s, t, u, v in table:
        polynomial[int(s)] = polynomial[int(s)] + c * pd**t * aeao**u * z**v  # C J^s (P/D)^t (AE/AO)^u z^v (wageningen.m 35-38)
    return polynomial


def _polyval(coefficients, j):
    """Horner's scheme for ``sum_s c_s J^s`` (numerical method)."""
    result = 0.0
    for c in coefficients[::-1]:
        result = result * j + c  # Horner step (numerical method)
    return result


def wageningen_expression(advance_number, pitch_diameter_ratio, blade_area_ratio, blade_count,
                          clip_advance_ratio=False):
    """``(K_T, K_Q)`` for a symbolic or numeric ``J`` and geometry (the
    clip only with ``clip_advance_ratio``)."""
    j = advance_number
    if clip_advance_ratio:
        j = ca.fmin(ca.fabs(j), MAX_ADVANCE_NUMBER)  # min(|J|, 1.3): construction behind the selector
    geometry = (pitch_diameter_ratio, blade_area_ratio, blade_count)
    return (_polyval(_collapse(_THRUST_TABLE, *geometry), j),  # K_T (wageningen.m 35-36)
            _polyval(_collapse(_TORQUE_TABLE, *geometry), j))  # K_Q (wageningen.m 37-38)


def wageningen_casadi(*, clip_advance_ratio=False):
    """The block: ``(J, pitch_diameter_ratio, blade_area_ratio, blade_count)
    -> (KT, KQ)`` (1x1 each).

    Contract
    --------
    ``clip_advance_ratio`` a bool (default ``False`` = MSS). Inputs: ``J``
    (1x1) and the parameters of ``wageningen_parameters()`` by name.
    """
    if not isinstance(clip_advance_ratio, bool):
        raise ValueError("clip_advance_ratio must be True or False")
    p = symbols(_DECLARED)
    j = ca.SX.sym("J", 1)
    kt, kq = wageningen_expression(j, p["pitch_diameter_ratio"], p["blade_area_ratio"],
                                   p["blade_count"], clip_advance_ratio)
    names = [d.name for d in _DECLARED]
    return ca.Function("wageningen", [j, *[p[name] for name in names]], [kt, kq],
                       ["J", *names], ["KT", "KQ"])
