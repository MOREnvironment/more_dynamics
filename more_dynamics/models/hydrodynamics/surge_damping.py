"""Surge damping: linear at low speed, quadratic at high speed, blended.

Equations (keys in References)::

    A11   = 2.7 rho (m / rho)^(5/3) / L^2      (addedMassSurge.m 32-33 @ 72656d1)
    Xudot = -A11                               (forceSurgeDamping.m 60 @ ac77394)
    Xu    = -(m - Xudot) / T1 = -(m + A11) / T1
            (Fossen 2011, eqs. 6.71, p. 124, and 6.76, p. 125;
             forceSurgeDamping.m 61 @ ac77394)
    X     = sigma Xu u_r + (1 - sigma) Xuu |u_r| u_r
                                               (forceSurgeDamping.m 79-82 @ ac77394)

MSS attributes the surge added mass to Söding (1982), not read here.

Two calibrations of the quadratic coefficient, one preprocess each:

* ``max_thrust``: ``Xuu = -max_thrust / max_speed^2`` (the steady state at top
  speed, ``forceSurgeDamping.m`` 64-66 at ``ac77394``);
* ``ittc``: ``Xuu = -1/2 rho S (1 + k) C_f`` with the ITTC-1957 line
  ``C_f = 0.075 / (log10 Rn - 2)^2``, ``Rn = L |u_r| / nu``, ``nu = 1e-6``
  m^2/s (Fossen 2011, eqs. 6.82-6.85, p. 125, with ``C_R = 0``;
  ``forceSurgeDamping.m`` 69-74 at ``ac77394``; ``XuuITTC.m`` 32-39 at
  ``72656d1``), evaluated at the current speed.

MSS hard-codes the form factor ``k = 0.1`` (``forceSurgeDamping.m`` 70) and
the crossover speed ``u_cross = 2`` m/s (line 57); here both are parameters
(Fossen 2011, p. 125: ``k`` typically 0.1 in transit, 0.25 in DP).

Revision: ``LIBRARY/modeling/forceSurgeDamping.m`` is cited at MSS
``ac77394``, the last revision that holds it; it was deleted in MSS
``108ceda`` (release 2.0, 2026-10-06) and is absent from release 2.0.2
(``72656d1``). ``addedMassSurge.m`` and ``XuuITTC.m`` are cited at
``72656d1``.

Deviations from ``forceSurgeDamping.m`` (each with a flag that gives MSS's
line back):

* ``Xudot = -A11`` as MSS; the numpy source takes ``Xudot = +A11``, which is
  wrong in sign (the linearised mode ``(m + A11) du/dt = Xu u`` must decay
  with the time constant ``T1``, Fossen 2011, eqs. 6.71 and 6.76).
* ``surge_blend="symmetric"`` (default): ``sigma = 1 - tanh(|u_r| / u_cross)``.
  ``"mss_tanh"``: ``sigma = 1 - tanh(u_r / u_cross)`` (line 79), which exceeds
  1 for ``u_r < 0``, so the quadratic term changes sign and pushes the vehicle
  astern, against the dissipative property of damping (Fossen 2011,
  Property 6.3, p. 123); the file's own header (lines 8-12) describes a
  symmetric blend. The two are equal for ``u_r >= 0``.
* ``ittc_reynolds_floor=1e5`` (default): ``Rn = max(L |u_r| / nu, floor)``,
  the bound and the line of MSS's own ``XuuITTC.m`` 34-36 (Fossen 2011,
  p. 125: a minimum ``Rn`` should be used, ``C_F`` blows up at low speed).
  ``None``: the unbounded ``log10(Rn + 1e-10)`` of ``forceSurgeDamping.m``
  71-73, which has a pole at ``Rn = 100``. The two are equal for
  ``Rn >= floor`` up to the ``1e-10``.
* Inputs outside their domain raise ``ValueError`` naming the input.

Conventions: ``nu_r = [u v w p q r]`` relative to the water; only ``u_r``
enters. ``tau`` is the force on the vehicle (added to the right-hand side).
A vehicle that uses this block sets the surge term of its linear damping
matrix to zero (MSS ``osv.m`` 188 at ``72656d1``), or surge damping is
counted twice.

Ported from the numpy source [MGM] ``dynamics/plant/matrices/linear_damping.py``
``DampingMatrix.force_surge_damping`` (230-294) and ``added_mass.py``
``AddedMass.added_mass_surge_static`` (338-357).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 6,
    Property 6.3, p. 123; eqs. 6.71-6.85, pp. 124-125.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/addedMassSurge.m`` 32-33, ``XuuITTC.m`` 32-39;
    ``CRAFT/SHIP/models/osv.m`` 188.
[MSS ac77394] Fossen, T. I. *Marine Systems Simulator (MSS)*,
    https://github.com/cybergalactic/MSS, MIT licence, revision ``ac77394``
    (2026-10-05): ``LIBRARY/modeling/forceSurgeDamping.m`` 57-82 (revision
    of 2025-09-23; the file is deleted in later revisions).
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/``
    files and lines listed above.
"""

from dataclasses import dataclass
from typing import Optional

import casadi as ca
import numpy as np

# addedMassSurge.m 33 (72656d1).
SURGE_ADDED_MASS_FACTOR = 2.7

# forceSurgeDamping.m 69 (ac77394), XuuITTC.m 32 (72656d1): nu (m^2/s); Fossen 2011, p. 125.
KINEMATIC_VISCOSITY = 1e-6

# C_f = 0.075 / (log10 Rn - 2)^2 (Fossen 2011, eq. 6.83, p. 125;
# forceSurgeDamping.m 73 at ac77394; XuuITTC.m 36 at 72656d1).
ITTC_FRICTION_FACTOR = 0.075
ITTC_LOG10_OFFSET = 2.0

# forceSurgeDamping.m 71 (ac77394): added to Rn inside the logarithm (MSS line only).
MSS_REYNOLDS_OFFSET = 1e-10

QUADRATIC_MODELS = ("ittc", "max_thrust")
SURGE_BLENDS = ("symmetric", "mss_tanh")


@dataclass(frozen=True)
class SurgeDampingConstants:
    """Coefficients of ``X = sigma Xu u_r + (1 - sigma) Xuu |u_r| u_r``.

    ``quadratic_coefficient`` (``Xuu``) is a number for ``max_thrust`` and NaN
    for ``ittc``, where ``Xuu = ittc_pressure_factor C_f(Rn)`` and
    ``Rn = reynolds_per_speed |u_r|`` are built in the graph.
    """

    added_mass: float
    linear_coefficient: float
    crossover_speed: float
    surge_blend: str
    quadratic_model: str
    quadratic_coefficient: float
    ittc_pressure_factor: float
    reynolds_per_speed: float
    ittc_reynolds_floor: Optional[float]


def _require_positive(**values: float) -> None:
    for name, value in values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")


def _linear_part(
    mass: float,
    length: float,
    water_density: float,
    time_constant: float,
    crossover_speed: float,
    surge_blend: str,
) -> tuple:
    """``(A11, Xu)``: ``addedMassSurge.m`` 32-33 (``72656d1``),
    ``forceSurgeDamping.m`` 60-61 (``ac77394``)."""
    _require_positive(
        mass=mass,
        length=length,
        water_density=water_density,
        time_constant=time_constant,
        crossover_speed=crossover_speed,
    )
    if surge_blend not in SURGE_BLENDS:
        raise ValueError(f"surge_blend must be one of {SURGE_BLENDS}, got {surge_blend!r}")
    displaced_volume = mass / water_density  # (addedMassSurge.m 32)
    added_mass = SURGE_ADDED_MASS_FACTOR * water_density * displaced_volume ** (5 / 3) / length**2  # (addedMassSurge.m 33)
    surge_acceleration_derivative = -added_mass  # (forceSurgeDamping.m 60 @ ac77394)
    # (Fossen 2011, eqs. 6.71, p. 124, and 6.76, p. 125; forceSurgeDamping.m 61 @ ac77394)
    linear_coefficient = -(mass - surge_acceleration_derivative) / time_constant
    return float(added_mass), float(linear_coefficient)


def preprocess_surge_damping_ittc(
    mass: float,
    length: float,
    wetted_surface: float,
    water_density: float,
    time_constant: float,
    form_factor: float,
    crossover_speed: float,
    surge_blend: str = "symmetric",
    ittc_reynolds_floor: Optional[float] = 1e5,
) -> SurgeDampingConstants:
    """``forceSurgeDamping(flag, u_r, m, S, L, T1, rho, u_max)``, the ITTC line.

    Source names: ``m``, ``L``, ``S``, ``rho``, ``T1``, ``k``, ``u_cross``.
    ``ittc_reynolds_floor``: lower bound of the Reynolds number (default: the
    ``Re_min`` of ``XuuITTC.m`` 34 at ``72656d1``); it must exceed 100, the
    pole of the line.
    ``None`` gives ``forceSurgeDamping.m``'s unbounded line.
    """
    added_mass, linear_coefficient = _linear_part(
        mass, length, water_density, time_constant, crossover_speed, surge_blend
    )
    _require_positive(wetted_surface=wetted_surface)
    if not np.isfinite(form_factor) or form_factor < 0.0:
        raise ValueError("form_factor must be a non-negative finite value")
    if ittc_reynolds_floor is not None:
        pole = 10.0**ITTC_LOG10_OFFSET
        if not np.isfinite(ittc_reynolds_floor) or ittc_reynolds_floor <= pole:
            raise ValueError(
                f"ittc_reynolds_floor must be a finite value above {pole:g} "
                "(the pole of the ITTC line), or None"
            )
        ittc_reynolds_floor = float(ittc_reynolds_floor)
    return SurgeDampingConstants(
        added_mass=added_mass,
        linear_coefficient=linear_coefficient,
        crossover_speed=float(crossover_speed),
        surge_blend=surge_blend,
        quadratic_model="ittc",
        quadratic_coefficient=float("nan"),
        # -1/2 rho S (1 + k) (Fossen 2011, eq. 6.82, p. 125; forceSurgeDamping.m 74 @ ac77394)
        ittc_pressure_factor=float(-0.5 * water_density * wetted_surface * (1 + form_factor)),
        reynolds_per_speed=float(length / KINEMATIC_VISCOSITY),  # Rn / |u_r| (Fossen 2011, eq. 6.85, p. 125)
        ittc_reynolds_floor=ittc_reynolds_floor,
    )


def preprocess_surge_damping_max_thrust(
    mass: float,
    length: float,
    water_density: float,
    time_constant: float,
    max_speed: float,
    max_thrust: float,
    crossover_speed: float,
    surge_blend: str = "symmetric",
) -> SurgeDampingConstants:
    """``forceSurgeDamping(flag, u_r, m, S, L, T1, rho, u_max, thrust_max)``.

    Source names: ``m``, ``L``, ``rho``, ``T1``, ``u_max``, ``thrust_max``,
    ``u_cross``. ``Xuu = -max_thrust / max_speed^2`` (``forceSurgeDamping.m``
    66 at ``ac77394``); the wetted surface of the MSS call is not used on this
    path.
    """
    added_mass, linear_coefficient = _linear_part(
        mass, length, water_density, time_constant, crossover_speed, surge_blend
    )
    _require_positive(max_speed=max_speed)
    if not np.isfinite(max_thrust) or max_thrust < 0.0:
        raise ValueError("max_thrust must be a non-negative finite value")
    return SurgeDampingConstants(
        added_mass=added_mass,
        linear_coefficient=linear_coefficient,
        crossover_speed=float(crossover_speed),
        surge_blend=surge_blend,
        quadratic_model="max_thrust",
        quadratic_coefficient=float(-max_thrust / max_speed**2),  # (forceSurgeDamping.m 66 @ ac77394)
        ittc_pressure_factor=float("nan"),
        reynolds_per_speed=float("nan"),
        ittc_reynolds_floor=None,
    )


def _ittc_quadratic_coefficient(constants: SurgeDampingConstants, u_r: ca.SX) -> ca.SX:
    """``Xuu(u_r)``: ``XuuITTC.m`` 35-36 (``72656d1``) with the floor, else
    ``forceSurgeDamping.m`` 72-74 (``ac77394``); Fossen 2011, eqs. 6.82-6.85,
    p. 125."""
    reynolds = constants.reynolds_per_speed * ca.fabs(u_r)  # (forceSurgeDamping.m 72 @ ac77394)
    if constants.ittc_reynolds_floor is None:
        reynolds = reynolds + MSS_REYNOLDS_OFFSET  # (forceSurgeDamping.m 73 @ ac77394)
    else:
        reynolds = ca.fmax(reynolds, constants.ittc_reynolds_floor)  # (XuuITTC.m 35)
    friction = ITTC_FRICTION_FACTOR / (ca.log10(reynolds) - ITTC_LOG10_OFFSET) ** 2  # (Fossen 2011, eq. 6.83, p. 125)
    return constants.ittc_pressure_factor * friction  # (Fossen 2011, eq. 6.82, p. 125; XuuITTC.m 39)


def surge_damping_casadi(constants: SurgeDampingConstants) -> ca.Function:
    """``nu_r -> tau = [X 0 0 0 0 0]`` (``forceSurgeDamping.m`` 79-82 at
    ``ac77394``)."""
    nu_r = ca.SX.sym("nu_r", 6)
    u_r = nu_r[0]
    if constants.quadratic_model == "ittc":
        quadratic_coefficient = _ittc_quadratic_coefficient(constants, u_r)
    else:
        quadratic_coefficient = ca.SX(constants.quadratic_coefficient)
    # symmetric |u_r| by default; u_r as forceSurgeDamping.m 79 (ac77394) with "mss_tanh"
    blend_speed = ca.fabs(u_r) if constants.surge_blend == "symmetric" else u_r
    sigma = 1 - ca.tanh(blend_speed / constants.crossover_speed)
    surge = (  # (forceSurgeDamping.m 82 @ ac77394)
        sigma * constants.linear_coefficient * u_r
        + (1 - sigma) * quadratic_coefficient * ca.fabs(u_r) * u_r
    )
    tau = ca.vertcat(surge, 0.0, 0.0, 0.0, 0.0, 0.0)
    return ca.Function("surge_damping", [nu_r], [tau], ["nu_r"], ["tau"])
