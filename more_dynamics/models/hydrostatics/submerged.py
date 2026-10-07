"""Submerged hydrostatics: restoring vector ``g(eta)`` (CasADi).

Equations (keys in References):

* ``W = m g``, ``B = rho g nabla`` (Fossen 2011, eq. 4.1, p. 59); a
  neutrally buoyant body has ``W = B`` (Fossen 2011, eq. 4.7, p. 61).
* ``g(eta) = -[f_g + f_b; r_bg x f_g + r_bb x f_b]`` with ``f_g = R^T [0, 0,
  W]``, ``f_b = -R^T [0, 0, B]`` (Fossen 2011, eq. 4.5, p. 60). With
  ``R(3,:) = [-s(theta), c(theta) s(phi), c(theta) c(phi)]`` (Fossen 2011,
  eq. 2.18, p. 22; MSS ``Rzyx.m`` 19) this is MSS ``gRvect.m`` 26-32:
  ``g[0:3] = -(W - B) R(3,:)^T``, ``g[3:6] = -m x R(3,:)^T`` with the first
  moment ``m = r_bg W - r_bb B``; expanded in Euler angles it is Fossen 2011,
  eq. 4.6, p. 60 (MSS ``gvect.m`` 22-28).

``g`` is the left-hand-side restoring vector: the equation of motion is
``M nu_dot + ... + g = tau`` (Fossen 2011, eq. 4.5, p. 60; MSS
``remus100.m`` 258, ``... - g``), so the force applied to the vehicle is
``-g``. For ``B > W`` at level attitude ``g_z = B - W > 0`` and the applied
``-g_z`` points up (negative z in NED).

``R(3,:)`` is ``MatrixTransforms.Rzyx_row3`` of
``more_transformations.more_casadi_transformations``. Conventions: NED, z
down; ``eta = [x y z phi theta psi]`` with zyx Euler angles;
``center_of_gravity`` (``r_bg``) and ``center_of_buoyancy`` (``r_bb``) are
measured from the CO in BODY. Yaw and position do not enter ``g``.

Ported from the numpy source [MGM] ``dynamics/plant/matrices/
restoring_forces.py``: ``RestoringForces.g_restoring_R`` (121-167, the
``gRvect.m`` form), ``g_restoring_submerged`` (19-66, the ``gvect.m`` Euler
form) and ``g_restoring_neutral`` (71-117, ``W = B``), as used by
``plant/auv_spheroid/auv_spheroid.py::AUVSpheroid.get_restoring`` (241-261).
The three source variants are one vector here.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 2,
    eq. 2.18, p. 22; Ch. 4, eqs. 4.1-4.7, pp. 59-61.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/gRvect.m`` 26-32, ``gvect.m`` 22-28;
    ``LIBRARY/kinematics/Rzyx.m`` 16-19; ``CRAFT/AUV/models/remus100.m``
    231, 258.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/`` files and
    lines listed above.
"""

from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np
from more_transformations.more_casadi_transformations.matrix_transforms import (
    MatrixTransforms,
)


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
        net_weight=float(weight - buoyancy),  # W - B (gRvect.m 27-29)
        first_moment=r_bg * weight - r_bb * buoyancy,  # r_bg W - r_bb B (gRvect.m 30-32)
    )


def submerged_hydrostatics_casadi(
    constants: SubmergedHydrostaticsConstants,
) -> ca.Function:
    """``eta -> g`` (6x1, N and N m): Fossen 2011, eq. 4.5, p. 60; MSS
    ``gRvect.m`` 26-32.

    The force on the vehicle is ``-g``.
    """
    eta = ca.SX.sym("eta", 6)
    phi, theta = eta[3], eta[4]

    # Third row of R = Rzyx(phi, theta, psi), body to NED (Fossen 2011, eq. 2.18, p. 22; Rzyx.m 19)
    row3 = MatrixTransforms.Rzyx_row3(phi, theta)
    r31, r32, r33 = row3[0], row3[1], row3[2]

    net = constants.net_weight
    mx, my, mz = (float(value) for value in constants.first_moment)

    # (Fossen 2011, eqs. 4.5-4.6, p. 60; gRvect.m 26-32, one line per row)
    g = ca.vertcat(
        -net * r31,  # gRvect.m 27
        -net * r32,  # gRvect.m 28
        -net * r33,  # gRvect.m 29
        -my * r33 + mz * r32,  # gRvect.m 30
        -mz * r31 + mx * r33,  # gRvect.m 31
        -mx * r32 + my * r31,  # gRvect.m 32
    )
    return ca.Function("submerged_hydrostatics", [eta], [g], ["eta"], ["g"])
