"""Surface hydrostatics: linear restoring ``g = G eta`` of a floating craft (CasADi).

Our own module; the other ``linear_surface.py`` beside it is a different
block and is not used here. The craft is described by its geometry, never by
a vehicle name: ``hull_count`` identical hulls (1 or 2) of length ``L``, beam
``B`` and waterplane coefficient ``Cw``, placed at ``+-hull_lateral_offset``
``y``; ``draft`` ``T`` and ``displacement_volume`` ``nabla`` of the whole
craft.

Equations (keys in References)::

    A_hull = Cw L B                                  (otter.m 173; exShip 32)
    A_wp   = n A_hull
    I_T    = n (1/12) L B^3 6 Cw^3 / ((1 + Cw)(1 + 2 Cw)) + n A_hull y^2
             (Munro-Smith coefficient, exShipHydrostatics.m 38, 47;
              otter.m 174-175; box bound Fossen 2011, eq. 4.35, p. 66)
    I_L    = n c_L (1/12) B L^3                      (otter.m 176; exShip 48;
                                                      Fossen 2011, eq. 4.35, p. 66)
    KB     = (1/3)(5 T / 2 - (nabla / n) / A_KB)     (Morrish; exShip 35,
                                                      osv.m 84, otter.m 177)
    BM = I / nabla                                   (Fossen 2011, eq. 4.33, p. 66)
    KG = T - z_g,  GM = KB + BM - KG                 (Fossen 2011, eq. 4.32, p. 65;
                                                      exShip 41, 55-56; otter.m 178-184)
    G_CF   = diag(0, 0, rho g A_wp, rho g nabla GM_T, rho g nabla GM_L, 0)
                                                     (Fossen 2011, eq. 4.24, p. 65;
                                                      Gmtrx.m 32-35)
    G_CO   = H(r_bF)^T G_CF H(r_bF),  r_bF = [LCF 0 0]
                                                     (Fossen 2011, eq. 7.250, p. 181;
                                                      Gmtrx.m 29, 36)
    G      = H(r_bP)^T G_CO H(r_bP)                  (Gmtrx.m 37)

MSS attributes the Morrish and Munro-Smith formulas to Fossen (2027), 3rd
ed., eqs. 4.37-4.38 (``exShipHydrostatics.m`` 35, 38); that edition is not on
disk and they are cited here from the MSS lines. The numpy source's three
scale factors ``s`` are, in this geometry, ``[n, n c_L, 1/n]``.

``center_of_buoyancy_area`` chooses ``A_KB``:

* ``"waterplane"`` (default): ``A_KB = A_hull``, Morrish with the waterplane
  area, as MSS ``exShipHydrostatics.m`` 35 and ``osv.m`` 84. A wall-sided
  hull then has its centre of buoyancy at exactly ``T/2``.
* ``"length_times_beam"``: ``A_KB = L B``, the form of MSS ``otter.m`` 177,
  kept to reproduce MSS's twin-hull craft. Deviation from ``otter.m`` 177 by
  default: that line divides by ``L B`` of a pontoon, which drops ``Cw``;
  for a wall-sided (prism) hull it puts the centre of buoyancy at
  ``0.583 T`` when ``Cw = 0.75`` instead of ``T/2``, while MSS's own
  ``exShipHydrostatics.m`` 35 and ``osv.m`` 84 use the waterplane area.

For a single hull the default equals MSS.

Gravity is the WGS-84 normal gravity at ``latitude`` (degrees) from
``more_transformations``; ``gravity=`` overrides it with a number (9.81 to
compare with MSS, ``Gmtrx.m`` 26). ``H`` is ``MatrixTransforms.H_matrix``
from ``more_transformations`` (Fossen 2011, eq. 3.24, p. 49; ``Hmtrx.m``
16-18).

``g`` is the left-hand-side restoring vector: the equation of motion is
``M nu_dot + ... + G eta = tau`` (Fossen 2011, eq. 4.25, p. 65; MSS
``otter.m`` 261, ``... - G * eta``), so the force applied to the craft is
``-g``.

Conventions: NED, z down; the CO on the waterline (``KG = T - z_g``, MSS
``exShipHydrostatics.m`` 41); ``center_of_gravity``,
``longitudinal_center_of_flotation`` and ``reference_point`` are measured from
the CO in BODY; small roll, pitch and heave (Fossen 2011, eq. 4.22, p. 64).
Deviation from MSS: inputs outside their domain raise ``ValueError`` naming
the input, where MSS returns complex or non-finite numbers (a negative
stiffness under a square root, a zero dimension).

Ported from the numpy source [MGM] ``dynamics/plant/matrices/
restoring_forces.py``: ``RestoringForces.G_restoring_surface_vessel``
(170-228), ``get_coeff_surface_vessel`` (230-288) and
``get_metacenter_heights_coeff`` (290-307).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    eq. 3.24, p. 49; Ch. 4, eqs. 4.22-4.35, pp. 64-66; Ch. 7, eq. 7.250,
    p. 181.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/Gmtrx.m`` 25-37; ``LIBRARY/kinematics/Hmtrx.m``
    16-18; ``CRAFT/USV/models/otter.m`` 173-191, 261;
    ``mssExamples/exShipHydrostatics.m`` 32-56; ``CRAFT/SHIP/models/osv.m``
    84, 92.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/
    restoring_forces.py``, lines listed above.
"""

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import casadi as ca
import numpy as np
from more_transformations.ecef_ned_transforms import ECEFNEDtransform
from more_transformations.matrix_transforms import MatrixTransforms

CENTER_OF_BUOYANCY_AREAS = ("waterplane", "length_times_beam")


@dataclass(frozen=True)
class SurfaceHydrostaticsConstants:
    """Constants computed before building the CasADi graph.

    Lengths in m, areas in m^2, inertias in m^4; ``stiffness_at_flotation``
    (G_CF), ``stiffness_at_origin`` (G_CO) and ``stiffness_matrix`` (G, at the
    reference point) are 6x6 in N/m, N, N m/rad.
    """

    gravity: float
    water_density: float
    displacement_volume: float
    center_of_buoyancy_area: str
    waterplane_area: float
    transverse_waterplane_inertia: float
    longitudinal_waterplane_inertia: float
    center_of_buoyancy_above_keel: float
    transverse_metacentric_height: float
    longitudinal_metacentric_height: float
    longitudinal_center_of_flotation: float
    reference_point: np.ndarray
    stiffness_at_flotation: np.ndarray
    stiffness_at_origin: np.ndarray
    stiffness_matrix: np.ndarray


def _vector(values: Sequence[float], size: int, name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (size,):
        raise ValueError(f"{name} must contain exactly {size} values")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values")
    return vector


def _positive(values: dict) -> None:
    for name, value in values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")


def _finite(values: dict) -> None:
    for name, value in values.items():
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite")


def surface_restoring_matrices(
    water_density: float,
    gravity: float,
    displacement_volume: float,
    waterplane_area: float,
    transverse_metacentric_height: float,
    longitudinal_metacentric_height: float,
    longitudinal_center_of_flotation: float,
    reference_point: Sequence[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(G, G_CO, G_CF)``: Fossen 2011, eqs. 4.24, p. 65, and 7.250, p. 181;
    MSS ``Gmtrx.m`` 29-37 with ``rho`` and ``g`` as inputs (``Gmtrx.m`` 25-26
    fixes them).

    Source row ``RestoringForces.G_restoring_surface_vessel``. ``G_CF`` is
    diagonal at the centre of flotation ``[LCF 0 0]``, ``G_CO`` at the CO and
    ``G`` at ``reference_point`` (from the CO, BODY). A negative metacentric
    height is accepted (an unstable craft).
    """
    _positive({
        "water_density": water_density,
        "gravity": gravity,
        "displacement_volume": displacement_volume,
        "waterplane_area": waterplane_area,
    })
    _finite({
        "transverse_metacentric_height": transverse_metacentric_height,
        "longitudinal_metacentric_height": longitudinal_metacentric_height,
        "longitudinal_center_of_flotation": longitudinal_center_of_flotation,
    })
    r_bp = _vector(reference_point, 3, "reference_point")

    rho_g = water_density * gravity
    # G_CF (Fossen 2011, eq. 4.24, p. 65; Gmtrx.m 32-35)
    g_cf = np.diag([
        0.0,
        0.0,
        rho_g * waterplane_area,
        rho_g * displacement_volume * transverse_metacentric_height,
        rho_g * displacement_volume * longitudinal_metacentric_height,
        0.0,
    ])
    h_f = MatrixTransforms.H_matrix(
        np.array([longitudinal_center_of_flotation, 0.0, 0.0])
    )
    g_co = h_f.T @ g_cf @ h_f  # (Fossen 2011, eq. 7.250, p. 181; Gmtrx.m 29, 36)
    h_p = MatrixTransforms.H_matrix(r_bp)
    g_p = h_p.T @ g_co @ h_p  # (Gmtrx.m 37)
    return g_p, g_co, g_cf


def metacentric_heights(
    displacement_volume: float,
    center_of_buoyancy_above_keel: float,
    transverse_waterplane_inertia: float,
    longitudinal_waterplane_inertia: float,
    center_of_gravity: Sequence[float],
    draft: float,
) -> np.ndarray:
    """``[BM_T, BM_L, KM_T, KM_L, GM_T, GM_L]`` in m.

    Source row ``RestoringForces.get_metacenter_heights_coeff``:
    ``BM = I / nabla`` (Fossen 2011, eq. 4.33, p. 66), ``KM = KB + BM``,
    ``KG = draft - z_g`` (CO on the waterline), ``GM = KM - KG`` (Fossen 2011,
    eq. 4.32, p. 65; MSS ``otter.m`` 178-184).
    """
    _positive({
        "displacement_volume": displacement_volume,
        "transverse_waterplane_inertia": transverse_waterplane_inertia,
        "longitudinal_waterplane_inertia": longitudinal_waterplane_inertia,
        "draft": draft,
    })
    _finite({"center_of_buoyancy_above_keel": center_of_buoyancy_above_keel})
    r_bg = _vector(center_of_gravity, 3, "center_of_gravity")

    transverse_radius = transverse_waterplane_inertia / displacement_volume  # (Fossen 2011, eq. 4.33, p. 66; otter.m 178)
    longitudinal_radius = longitudinal_waterplane_inertia / displacement_volume  # (Fossen 2011, eq. 4.33, p. 66; otter.m 179)
    center_of_gravity_above_keel = draft - r_bg[2]  # (otter.m 182)
    transverse_metacentre = center_of_buoyancy_above_keel + transverse_radius  # (otter.m 180)
    longitudinal_metacentre = center_of_buoyancy_above_keel + longitudinal_radius  # (otter.m 181)
    # GM = KM - KG (Fossen 2011, eq. 4.32, p. 65; otter.m 183-184)
    return np.array([
        transverse_radius,
        longitudinal_radius,
        transverse_metacentre,
        longitudinal_metacentre,
        transverse_metacentre - center_of_gravity_above_keel,
        longitudinal_metacentre - center_of_gravity_above_keel,
    ])


def preprocess_surface_hydrostatics(
    length: float,
    hull_beam: float,
    draft: float,
    displacement_volume: float,
    waterplane_area_coefficient: float,
    center_of_gravity: Sequence[float],
    hull_count: int,
    hull_lateral_offset: float,
    longitudinal_inertia_factor: float,
    longitudinal_center_of_flotation: float,
    reference_point: Sequence[float],
    water_density: float,
    latitude: float,
    gravity: Optional[float] = None,
    center_of_buoyancy_area: str = "waterplane",
) -> SurfaceHydrostaticsConstants:
    """Geometry to ``G`` (source row ``RestoringForces.get_coeff_surface_vessel``).

    ``length``, ``hull_beam`` (one hull), ``draft``, ``hull_lateral_offset``
    (hull centreline from the CO) in m; ``displacement_volume`` of the whole
    craft in m^3; ``waterplane_area_coefficient`` ``Cw`` of one hull, in
    (0, 1]; ``hull_count`` 1 (on the centreline, offset 0) or 2;
    ``longitudinal_inertia_factor`` ``c_L`` scales the box value of ``I_L``;
    ``water_density`` in kg/m^3; ``latitude`` in degrees; ``gravity`` in
    m/s^2 overrides the latitude gravity when given.
    """
    if center_of_buoyancy_area not in CENTER_OF_BUOYANCY_AREAS:
        raise ValueError(
            f"center_of_buoyancy_area must be one of {CENTER_OF_BUOYANCY_AREAS}"
        )
    if hull_count not in (1, 2):
        raise ValueError("hull_count must be 1 or 2")
    _positive({
        "length": length,
        "hull_beam": hull_beam,
        "draft": draft,
        "displacement_volume": displacement_volume,
        "waterplane_area_coefficient": waterplane_area_coefficient,
        "longitudinal_inertia_factor": longitudinal_inertia_factor,
        "water_density": water_density,
    })
    if waterplane_area_coefficient > 1.0:
        raise ValueError("waterplane_area_coefficient must not exceed 1")
    _finite({"hull_lateral_offset": hull_lateral_offset})
    if hull_count == 1 and hull_lateral_offset != 0.0:
        raise ValueError("hull_lateral_offset must be 0 for a single hull")
    if hull_count == 2 and hull_lateral_offset <= 0.0:
        raise ValueError("hull_lateral_offset must be positive for two hulls")
    if not np.isfinite(latitude) or abs(latitude) > 90.0:
        raise ValueError("latitude must be finite degrees in [-90, 90]")
    if gravity is None:
        gravity = ECEFNEDtransform.gravity(latitude)
    else:
        _positive({"gravity": gravity})
    r_bg = _vector(center_of_gravity, 3, "center_of_gravity")

    n = hull_count
    cw = waterplane_area_coefficient
    hull_area = cw * length * hull_beam  # (otter.m 173; exShipHydrostatics.m 32)
    waterplane_area = n * hull_area
    # Munro-Smith, per hull, plus the parallel-axis term of each hull at +-y
    munro_smith = 6.0 * cw**3 / ((1.0 + cw) * (1.0 + 2.0 * cw))  # (exShipHydrostatics.m 38)
    # (exShipHydrostatics.m 47; otter.m 174-175)
    transverse_inertia = (
        n * length * hull_beam**3 / 12.0 * munro_smith
        + n * hull_area * hull_lateral_offset**2
    )
    # (Fossen 2011, eq. 4.35, p. 66; exShipHydrostatics.m 48; otter.m 176)
    longitudinal_inertia = (
        n * longitudinal_inertia_factor * hull_beam * length**3 / 12.0
    )
    # Morrish, per hull; the area is the waterplane area by default
    # (exShipHydrostatics.m 35, osv.m 84), L B with "length_times_beam" (otter.m 177)
    area_for_kb = (
        hull_area if center_of_buoyancy_area == "waterplane"
        else length * hull_beam
    )
    center_of_buoyancy_above_keel = (
        2.5 * draft - displacement_volume / n / area_for_kb
    ) / 3.0

    heights = metacentric_heights(
        displacement_volume=displacement_volume,
        center_of_buoyancy_above_keel=center_of_buoyancy_above_keel,
        transverse_waterplane_inertia=transverse_inertia,
        longitudinal_waterplane_inertia=longitudinal_inertia,
        center_of_gravity=r_bg,
        draft=draft,
    )
    transverse_metacentric_height = float(heights[4])
    longitudinal_metacentric_height = float(heights[5])

    stiffness_matrix, stiffness_at_origin, stiffness_at_flotation = (
        surface_restoring_matrices(
            water_density=water_density,
            gravity=gravity,
            displacement_volume=displacement_volume,
            waterplane_area=waterplane_area,
            transverse_metacentric_height=transverse_metacentric_height,
            longitudinal_metacentric_height=longitudinal_metacentric_height,
            longitudinal_center_of_flotation=longitudinal_center_of_flotation,
            reference_point=reference_point,
        )
    )

    return SurfaceHydrostaticsConstants(
        gravity=gravity,
        water_density=float(water_density),
        displacement_volume=float(displacement_volume),
        center_of_buoyancy_area=center_of_buoyancy_area,
        waterplane_area=float(waterplane_area),
        transverse_waterplane_inertia=float(transverse_inertia),
        longitudinal_waterplane_inertia=float(longitudinal_inertia),
        center_of_buoyancy_above_keel=float(center_of_buoyancy_above_keel),
        transverse_metacentric_height=transverse_metacentric_height,
        longitudinal_metacentric_height=longitudinal_metacentric_height,
        longitudinal_center_of_flotation=float(longitudinal_center_of_flotation),
        reference_point=_vector(reference_point, 3, "reference_point"),
        stiffness_at_flotation=stiffness_at_flotation,
        stiffness_at_origin=stiffness_at_origin,
        stiffness_matrix=stiffness_matrix,
    )


def surface_hydrostatics_casadi(
    constants: SurfaceHydrostaticsConstants,
) -> ca.Function:
    """``eta -> (g, G)``: ``g = G eta`` (6x1, N and N m) and the constant ``G``.

    ``eta = [x y z phi theta psi]`` (NED, m and rad). ``g`` enters the
    equation of motion on the left, so the force on the craft is ``-g``.
    """
    eta = ca.SX.sym("eta", 6)
    stiffness = ca.SX(ca.DM(constants.stiffness_matrix))
    g = ca.mtimes(stiffness, eta)
    return ca.Function(
        "surface_hydrostatics",
        [eta],
        [g, stiffness],
        ["eta"],
        ["g", "G"],
    )
