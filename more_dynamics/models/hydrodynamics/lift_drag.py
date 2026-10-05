"""Lift and drag of a slender body (numpy constants, CasADi algebra).

Ported from ``more_generic_models`` ``dynamics/plant/matrices/hydrodynamics.py``:
``HydroForces.get_lift_drag_coeff`` (lines 62-75, ``sigma = 0``) and
``force_lift_drag`` (78-88); translated against MSS (MIT, T. I. Fossen)
``LIBRARY/modeling/forceLiftDrag.m`` lines 26-40, ``coeffLiftDrag.m`` lines
54-68 and ``remus100.m`` lines 127-128. MSS hard-codes ``rho = 1026``
(``forceLiftDrag.m`` 26) and ``e = 0.3`` (``coeffLiftDrag.m`` 54); here both
are parameters.

Lift and drag act in the x-z plane; ``alpha = atan2(w_r, u_r)``,
``U_r = |nu_r[0:3]|``, linear lift ``C_L = C_L_alpha alpha`` (no stall
blending), ``C_D = C_D0 + C_L^2 / (pi e AR)``.
"""

from dataclasses import dataclass

import casadi as ca
import numpy as np


@dataclass(frozen=True)
class LiftDragConstants:
    """Lift slope, induced-drag factor and dynamic-pressure area."""

    lift_slope: float
    parasitic_drag_coefficient: float
    induced_drag_factor: float
    pressure_area: float


def preprocess_lift_drag(
    span: float,
    planform_area: float,
    parasitic_drag_coefficient: float,
    oswald_efficiency: float,
    water_density: float,
) -> LiftDragConstants:
    """Coefficients of ``forceLiftDrag(b, S, CD_0, alpha, U_r)``.

    Source names: ``b``, ``S``, ``CD_0``, ``e``, ``rho``.
    ``AR = b^2 / S``; ``C_L_alpha = pi AR / (1 + sqrt(1 + (AR/2)^2))``.
    """
    positive = {
        "span": span,
        "planform_area": planform_area,
        "oswald_efficiency": oswald_efficiency,
        "water_density": water_density,
    }
    for name, value in positive.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")
    if not np.isfinite(parasitic_drag_coefficient) or parasitic_drag_coefficient < 0.0:
        raise ValueError("parasitic_drag_coefficient must be a non-negative finite value")

    aspect_ratio = span * span / planform_area
    return LiftDragConstants(
        lift_slope=np.pi * aspect_ratio / (1.0 + np.sqrt(1.0 + (0.5 * aspect_ratio) ** 2)),
        parasitic_drag_coefficient=float(parasitic_drag_coefficient),
        induced_drag_factor=1.0 / (np.pi * oswald_efficiency * aspect_ratio),
        pressure_area=0.5 * water_density * planform_area,
    )


def lift_drag_casadi(constants: LiftDragConstants) -> ca.Function:
    """``nu_r -> tau = [X 0 Z 0 0 0]`` (``forceLiftDrag.m`` 30-40)."""
    nu_r = ca.SX.sym("nu_r", 6)
    alpha = ca.atan2(nu_r[2], nu_r[0])
    speed_squared = nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2
    lift_coefficient = constants.lift_slope * alpha
    drag_coefficient = (
        constants.parasitic_drag_coefficient
        + lift_coefficient**2 * constants.induced_drag_factor
    )
    pressure = constants.pressure_area * speed_squared
    sin_alpha, cos_alpha = ca.sin(alpha), ca.cos(alpha)
    surge = -pressure * (drag_coefficient * cos_alpha - lift_coefficient * sin_alpha)
    heave = -pressure * (drag_coefficient * sin_alpha + lift_coefficient * cos_alpha)
    tau = ca.vertcat(surge, 0.0, heave, 0.0, 0.0, 0.0)
    return ca.Function("lift_drag", [nu_r], [tau], ["nu_r"], ["tau"])
