"""Submerged hydrostatics: restoring vector ``g(eta)`` (CasADi).

Ported from ``more_generic_models`` ``dynamics/plant/matrices/restoring_forces.py``:
``RestoringForces.g_restoring_R`` (lines 121-167, the MSS ``gRvect.m`` form),
``g_restoring_submerged`` (19-66, the ``gvect.m`` Euler-angle form) and
``g_restoring_neutral`` (71-117, ``W = B``), as used by
``auv_spheroid.py::AUVSpheroid.get_restoring`` (241-261).

The three source variants are one vector: ``gvect.m`` is ``gRvect.m`` with the
third row of ``Rzyx`` written out, and the neutral case is ``buoyancy =
weight``. The block uses the ``gRvect.m`` form (MSS ``LIBRARY/modeling/gRvect.m``
lines 26-32) with ``R(3,:) = [-s(theta), c(theta)s(phi), c(theta)c(phi)]``
(``Rzyx.m`` line 19), so the translational force is ``-(W - B) R(3,:)^T``.

Conventions: NED, z down; ``eta = [x y z phi theta psi]`` with zyx Euler
angles; ``center_of_gravity`` (``r_bg``) and ``center_of_buoyancy``
(``r_bb``) are measured from the CO in BODY. Yaw and position do not enter
``g``.
"""

from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np


@dataclass(frozen=True)
class SubmergedHydrostaticsConstants:
    """Constants computed before building the CasADi graph.

    ``net_weight`` is ``W - B``; ``first_moment`` is ``r_bg W - r_bb B``.
    """

    weight: float
    buoyancy: float
    center_of_gravity: np.ndarray
    center_of_buoyancy: np.ndarray
    net_weight: float
    first_moment: np.ndarray


def _vector(values: Sequence[float], size: int, name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (size,):
        raise ValueError(f"{name} must contain exactly {size} values")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values")
    return vector


def preprocess_submerged_hydrostatics(
    weight: float,
    buoyancy: float,
    center_of_gravity: Sequence[float],
    center_of_buoyancy: Sequence[float],
) -> SubmergedHydrostaticsConstants:
    """Weight and buoyancy in N, centres in m from the CO (body frame).

    Source names: ``W``, ``B``, ``r_bg``, ``r_bb``. A neutrally buoyant body
    (``g_restoring_neutral``) is ``buoyancy = weight``.
    """
    for name, value in {"weight": weight, "buoyancy": buoyancy}.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")

    r_bg = _vector(center_of_gravity, 3, "center_of_gravity")
    r_bb = _vector(center_of_buoyancy, 3, "center_of_buoyancy")

    return SubmergedHydrostaticsConstants(
        weight=float(weight),
        buoyancy=float(buoyancy),
        center_of_gravity=r_bg,
        center_of_buoyancy=r_bb,
        net_weight=float(weight - buoyancy),
        first_moment=r_bg * weight - r_bb * buoyancy,
    )


def submerged_hydrostatics_casadi(
    constants: SubmergedHydrostaticsConstants,
) -> ca.Function:
    """``eta -> g`` (6x1, N and N m), MSS ``gRvect.m`` lines 26-32."""
    eta = ca.SX.sym("eta", 6)
    phi, theta = eta[3], eta[4]

    # Third row of R = Rzyx(phi, theta, psi), body to NED (Rzyx.m line 19)
    r31 = -ca.sin(theta)
    r32 = ca.cos(theta) * ca.sin(phi)
    r33 = ca.cos(theta) * ca.cos(phi)

    net = constants.net_weight
    mx, my, mz = (float(value) for value in constants.first_moment)

    g = ca.vertcat(
        -net * r31,
        -net * r32,
        -net * r33,
        -my * r33 + mz * r32,
        -mz * r31 + mx * r33,
        -mx * r32 + my * r31,
    )
    return ca.Function("submerged_hydrostatics", [eta], [g], ["eta"], ["g"])
