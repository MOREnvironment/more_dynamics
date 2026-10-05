"""Wageningen B-series open-water coefficients ``K_T(J)``, ``K_Q(J)``.

Ported from ``more_generic_models``
``dynamics/propulsion/propeller_models/wagenigen/wagenigen.py`` (lines 3-127):
the regression of Barnitsas, Ray and Kinley (1981) at Re = 2e6, the same
tables as MSS ``LIBRARY/modeling/utiles/WageningData.mat`` (``wageningen.m``).
Rows are ``[C, s, t, u, v]``: ``K = sum C J^s (P/D)^t (A_E/A_O)^u z^v``. The
tables below are copied by script from the numpy source (lines 14-105).

Default (owner decision E-22 Q1 a): MSS ``wageningen.m``, the polynomial at
any ``J``, including ``J < 0`` and ``J > 1.3`` (outside the fitted range).
``clip_advance_ratio=True`` keeps the numpy source's behaviour: ``J`` is
replaced by ``min(|J|, 1.3)`` (line 11).
"""

from dataclasses import dataclass

import casadi as ca
import numpy as np

from ._common import positive

MAX_ADVANCE_NUMBER = 1.3  # numpy source line 11

_THRUST_TABLE = np.array([
    [0.00880496, 0, 0, 0, 0],
    [-0.204554, 1, 0, 0, 0],
    [0.166351, 0, 1, 0, 0],
    [0.158114, 0, 2, 0, 0],
    [-0.147581, 2, 0, 1, 0],
    [-0.481497, 1, 1, 1, 0],
    [0.415437, 0, 2, 1, 0],
    [0.0144043, 0, 0, 0, 1],
    [-0.0530054, 2, 0, 0, 1],
    [0.0143481, 0, 1, 0, 1],
    [0.0606826, 1, 1, 0, 1],
    [-0.0125894, 0, 0, 1, 1],
    [0.0109689, 1, 0, 1, 1],
    [-0.133698, 0, 3, 0, 0],
    [0.00638407, 0, 6, 0, 0],
    [-0.00132718, 2, 6, 0, 0],
    [0.168496, 3, 0, 1, 0],
    [-0.0507214, 0, 0, 2, 0],
    [0.0854559, 2, 0, 2, 0],
    [-0.0504475, 3, 0, 2, 0],
    [0.010465, 1, 6, 2, 0],
    [-0.00648272, 2, 6, 2, 0],
    [-0.00841728, 0, 3, 0, 1],
    [0.0168424, 1, 3, 0, 1],
    [-0.00102296, 3, 3, 0, 1],
    [-0.0317791, 0, 3, 1, 1],
    [0.018604, 1, 0, 2, 1],
    [-0.00410798, 0, 2, 2, 1],
    [-0.000606848, 0, 0, 0, 2],
    [-0.0049819, 1, 0, 0, 2],
    [0.0025983, 2, 0, 0, 2],
    [-0.000560528, 3, 0, 0, 2],
    [-0.00163652, 1, 2, 0, 2],
    [-0.000328787, 1, 6, 0, 2],
    [0.000116502, 2, 6, 0, 2],
    [0.000690904, 0, 0, 1, 2],
    [0.00421749, 0, 3, 1, 2],
    [0.0000565229, 3, 6, 1, 2],
    [-0.00146564, 0, 3, 2, 2],
])

_TORQUE_TABLE = np.array([
    [0.00379368, 0, 0, 0, 0],
    [0.00886523, 2, 0, 0, 0],
    [-0.032241, 1, 1, 0, 0],
    [0.00344778, 0, 2, 0, 0],
    [-0.0408811, 0, 1, 1, 0],
    [-0.108009, 1, 1, 1, 0],
    [-0.0885381, 2, 1, 1, 0],
    [0.188561, 0, 2, 1, 0],
    [-0.00370871, 1, 0, 0, 1],
    [0.00513696, 0, 1, 0, 1],
    [0.0209449, 1, 1, 0, 1],
    [0.00474319, 2, 1, 0, 1],
    [-0.00723408, 2, 0, 1, 1],
    [0.00438388, 1, 1, 1, 1],
    [-0.0269403, 0, 2, 1, 1],
    [0.0558082, 3, 0, 1, 0],
    [0.0161886, 0, 3, 1, 0],
    [0.00318086, 1, 3, 1, 0],
    [0.015896, 0, 0, 2, 0],
    [0.0471729, 1, 0, 2, 0],
    [0.0196283, 3, 0, 2, 0],
    [-0.0502782, 0, 1, 2, 0],
    [-0.030055, 3, 1, 2, 0],
    [0.0417122, 2, 2, 2, 0],
    [-0.0397722, 0, 3, 2, 0],
    [-0.00350024, 0, 6, 2, 0],
    [-0.0106854, 3, 0, 0, 1],
    [0.00110903, 3, 3, 0, 1],
    [-0.000313912, 0, 6, 0, 1],
    [0.0035985, 3, 0, 1, 1],
    [-0.00142121, 0, 6, 1, 1],
    [-0.00383637, 1, 0, 2, 1],
    [0.0126803, 0, 2, 2, 1],
    [-0.00318278, 2, 3, 2, 1],
    [0.00334268, 0, 6, 2, 1],
    [-0.00183491, 1, 1, 0, 2],
    [0.000112451, 3, 2, 0, 2],
    [-0.0000297228, 3, 6, 0, 2],
    [0.000269551, 1, 0, 1, 2],
    [0.00083265, 2, 0, 1, 2],
    [0.00155334, 0, 2, 1, 2],
    [0.000302683, 0, 6, 1, 2],
    [-0.0001843, 0, 0, 2, 2],
    [-0.000425399, 0, 3, 2, 2],
    [0.0000869243, 3, 3, 2, 2],
    [-0.0004659, 0, 6, 2, 2],
    [0.0000554194, 1, 6, 2, 2],
])


@dataclass(frozen=True)
class WageningenConstants:
    """``K_T = sum_s thrust_polynomial[s] J^s``, same for ``K_Q``, ``s = 0..3``."""

    pitch_diameter_ratio: float
    blade_area_ratio: float
    blade_count: float
    clip_advance_ratio: bool
    thrust_polynomial: np.ndarray
    torque_polynomial: np.ndarray


def _collapse(table: np.ndarray, pd: float, aeao: float, z: float) -> np.ndarray:
    """Sum the geometry factors of each power of ``J`` (numpy, preprocess)."""
    powers = table[:, 1].astype(int)
    weights = table[:, 0] * pd ** table[:, 2] * aeao ** table[:, 3] * z ** table[:, 4]
    polynomial = np.zeros(powers.max() + 1)
    np.add.at(polynomial, powers, weights)
    return polynomial


def preprocess_wageningen(
    pitch_diameter_ratio: float,
    blade_area_ratio: float,
    blade_count: float,
    clip_advance_ratio: bool = False,
) -> WageningenConstants:
    """Source names ``PD``, ``AEAO``, ``z``; ``clip_advance_ratio`` selects the
    numpy source's ``min(|J|, 1.3)`` instead of MSS's unclipped ``J``."""
    pd = positive("pitch_diameter_ratio", pitch_diameter_ratio)
    aeao = positive("blade_area_ratio", blade_area_ratio)
    z = positive("blade_count", blade_count)
    return WageningenConstants(
        pitch_diameter_ratio=pd,
        blade_area_ratio=aeao,
        blade_count=z,
        clip_advance_ratio=bool(clip_advance_ratio),
        thrust_polynomial=_collapse(_THRUST_TABLE, pd, aeao, z),
        torque_polynomial=_collapse(_TORQUE_TABLE, pd, aeao, z),
    )


def _polyval(coefficients: np.ndarray, j):
    result = 0.0
    for c in coefficients[::-1]:
        result = result * j + float(c)
    return result


def wageningen_values(constants: WageningenConstants, advance_number: float):
    """Numpy ``(K_T, K_Q)`` for preprocess steps (clipped only on the flag)."""
    j = float(advance_number)
    if constants.clip_advance_ratio:
        j = min(abs(j), MAX_ADVANCE_NUMBER)
    return _polyval(constants.thrust_polynomial, j), _polyval(constants.torque_polynomial, j)


def wageningen_expression(constants: WageningenConstants, advance_number):
    """Symbolic ``(K_T, K_Q)`` of a CasADi expression ``J`` (clipped only on the flag)."""
    j = advance_number
    if constants.clip_advance_ratio:
        j = ca.fmin(ca.fabs(j), MAX_ADVANCE_NUMBER)
    return _polyval(constants.thrust_polynomial, j), _polyval(constants.torque_polynomial, j)


def wageningen_casadi(constants: WageningenConstants) -> ca.Function:
    """``J -> (KT, KQ)``."""
    j = ca.SX.sym("J", 1)
    kt, kq = wageningen_expression(constants, j)
    return ca.Function("wageningen", [j], [kt, kq], ["J"], ["KT", "KQ"])
