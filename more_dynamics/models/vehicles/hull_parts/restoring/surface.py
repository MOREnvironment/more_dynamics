"""Surface hydrostatics: linear restoring ``g = G eta`` of a floating craft (CasADi).

Our own module; the other ``linear_surface.py`` beside it is a different
block and is not used here. The craft is described by its geometry, never by
a vehicle name: ``hull_count`` identical hulls (1 or 2) of length ``L``, beam
``B`` and waterplane coefficient ``Cw``, placed at ``+-hull_lateral_offset``
``y``; ``draft`` ``T`` and ``displacement_volume`` ``nabla`` of the whole
craft.

``surface_hydrostatics_parameters(...)`` declares the parameters of one
choice of the selectors (name, shape, SI unit, meaning, admissible range);
``surface_hydrostatics_casadi(...)`` builds the block ``(eta, <each declared
parameter by name>) -> (g, G, G_CO, G_CF, <coefficients>)``. Selectors
choose the graph: ``hull_count`` (1: one hull on the centreline, no
``hull_lateral_offset``; 2: twin hulls), ``gravity_source`` (``"latitude"``:
WGS-84 normal gravity at ``latitude`` in degrees from
``more_transformations``; ``"value"``: ``gravity`` given, 9.81 to compare
with MSS, ``Gmtrx.m`` 26) and ``center_of_buoyancy_area`` (below).

Equations (keys in References)::

    A_hull = Cw L B                                  (otter.m 174; exShip 32)
    A_wp   = n A_hull
    I_T    = n (1/12) L B^3 6 Cw^3 / ((1 + Cw)(1 + 2 Cw)) + n A_hull y^2
             (Munro-Smith coefficient, exShipHydrostatics.m 38, 47;
              otter.m 175-176; box bound Fossen 2011, eq. 4.35, p. 66)
    I_L    = n c_L (1/12) B L^3                      (otter.m 177; exShip 48;
                                                      Fossen 2011, eq. 4.35, p. 66)
    KB     = (1/3)(5 T / 2 - (nabla / n) / A_KB)     (Morrish; exShip 35,
                                                      osv.m 84, otter.m 178)
    BM = I / nabla                                   (Fossen 2011, eq. 4.33, p. 66)
    KG = T - z_g,  GM = KB + BM - KG                 (Fossen 2011, eq. 4.32, p. 65;
                                                      exShip 41, 55-56; otter.m 179-185)
    G_CF   = diag(0, 0, rho g A_wp, rho g nabla GM_T, rho g nabla GM_L, 0)
                                                     (Fossen 2011, eq. 4.24, p. 65;
                                                      Gmtrx.m 32-35)
    G_CO   = H(r_bF)^T G_CF H(r_bF),  r_bF = [LCF 0 0]
                                                     (Fossen 2011, eq. 7.250, p. 181;
                                                      Gmtrx.m 29, 36)
    G      = H(r_bP)^T G_CO H(r_bP)                  (Gmtrx.m 37)

MSS attributes the Morrish and Munro-Smith formulas to Fossen (2027), 3rd
ed., eqs. 4.37-4.38 (``exShipHydrostatics.m`` 35, 38); that edition is not on
disk and they are cited here from the MSS lines.

``center_of_buoyancy_area`` chooses ``A_KB``:

* ``"waterplane"`` (default): ``A_KB = A_hull``, Morrish with the waterplane
  area, as MSS ``exShipHydrostatics.m`` 35, ``osv.m`` 84 and, since MSS
  ``cf349d4``, ``otter.m`` 178 (revision ``cc07579``). A wall-sided hull then
  has its centre of buoyancy at exactly ``T/2``.
* ``"length_times_beam"``: ``A_KB = L B``, the form of MSS ``otter.m`` 177
  up to revision ``72656d1``, kept only to reproduce that earlier MSS: it
  divides by ``L B`` of a pontoon, which drops ``Cw``; for a wall-sided
  (prism) hull it puts the centre of buoyancy at ``0.583 T`` when
  ``Cw = 0.75`` instead of ``T/2``.

Outputs: ``g = G eta`` (6x1), ``G`` (6x6, at the reference point), ``G_CO``
(at the CO), ``G_CF`` (diagonal, at the centre of flotation), and the
coefficients ``waterplane_area``, ``transverse_waterplane_inertia``,
``longitudinal_waterplane_inertia``, ``center_of_buoyancy_above_keel``,
``transverse_metacentric_height``, ``longitudinal_metacentric_height`` and
``acceleration_of_gravity`` (the ``g`` the stiffness uses; not named
``gravity``, which is the input of ``gravity_source="value"``).

``g`` is the left-hand-side restoring vector: the equation of motion is
``M nu_dot + ... + G eta = tau`` (Fossen 2011, eq. 4.25, p. 65; MSS
``otter.m`` 262, ``... - G * eta``), so the force applied to the craft is
``-g``.

Conventions: NED, z down; the CO on the waterline (``KG = T - z_g``, MSS
``exShipHydrostatics.m`` 41); ``center_of_gravity``,
``longitudinal_center_of_flotation`` and ``reference_point`` are measured from
the CO in BODY; small roll, pitch and heave (Fossen 2011, eq. 4.22, p. 64). A
negative metacentric height is accepted (an unstable craft). The declared
ranges refuse inputs outside their domain, naming the input, where MSS
returns complex or non-finite numbers (a zero dimension).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    eq. 3.24, p. 49; Ch. 4, eqs. 4.22-4.35, pp. 64-66; Ch. 7, eq. 7.250,
    p. 181.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``LIBRARY/modeling/Gmtrx.m`` 25-37;
    ``LIBRARY/kinematics/Hmtrx.m`` 16-18; ``CRAFT/USV/models/otter.m``
    174-192, 262; ``mssExamples/exShipHydrostatics.m`` 32-56;
    ``CRAFT/SHIP/models/osv.m`` 84, 92; revision ``72656d1``:
    ``CRAFT/USV/models/otter.m`` 177 (the earlier KB).

Author:    Enio Krizman
Date:      2026-10-06
"""

import casadi as ca
from more_transformations.more_casadi_transformations import (
    ECEFNEDtransform,
    MatrixTransforms,
    Parameter,
    symbols,
)

CENTER_OF_BUOYANCY_AREAS = ("waterplane", "length_times_beam")
GRAVITY_SOURCES = ("latitude", "value")
HULL_COUNTS = (1, 2)

# Named outputs, in order.
SURFACE_HYDROSTATICS_OUTPUTS = (
    "g",
    "G",
    "G_CO",
    "G_CF",
    "waterplane_area",
    "transverse_waterplane_inertia",
    "longitudinal_waterplane_inertia",
    "center_of_buoyancy_above_keel",
    "transverse_metacentric_height",
    "longitudinal_metacentric_height",
    "acceleration_of_gravity",
)


def _positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


def _check_selectors(hull_count, gravity_source, center_of_buoyancy_area):
    if hull_count not in HULL_COUNTS or isinstance(hull_count, bool):
        raise ValueError(f"hull_count must be one of {HULL_COUNTS}, got {hull_count!r}")
    if gravity_source not in GRAVITY_SOURCES:
        raise ValueError(f"gravity_source must be one of {GRAVITY_SOURCES}, got {gravity_source!r}")
    if center_of_buoyancy_area not in CENTER_OF_BUOYANCY_AREAS:
        raise ValueError(
            f"center_of_buoyancy_area must be one of {CENTER_OF_BUOYANCY_AREAS}, "
            f"got {center_of_buoyancy_area!r}"
        )


def surface_hydrostatics_parameters(*, hull_count, gravity_source="latitude",
                                    center_of_buoyancy_area="waterplane"):
    """The declared parameter set of one choice of the selectors: a tuple of
    ``Parameter`` in the order of the block's inputs. ``hull_lateral_offset``
    exists only for ``hull_count=2``; ``latitude`` only for
    ``gravity_source="latitude"``, ``gravity`` only for ``"value"``. An
    unknown selector value raises ``ValueError``."""
    _check_selectors(hull_count, gravity_source, center_of_buoyancy_area)
    declared = [
        _positive("length", "m", "hull length L"),
        _positive("beam", "m", "beam B of one hull"),
        _positive("draft", "m", "draft T"),
        _positive("displacement_volume", "m^3", "displaced volume nabla of the whole craft"),
        Parameter("waterplane_area_coefficient", (1, 1), "1", "waterplane coefficient Cw = A_hull / (L B) of one hull",
                  0.0, 1.0, minimum_exclusive=True),
        Parameter("center_of_gravity", (3, 1), "m", "CO -> CG r_bg, body axes (FRD); the CO on the waterline"),
    ]
    if hull_count == 2:
        declared.append(_positive("hull_lateral_offset", "m", "hull centreline y from the CO (each hull at +-y)"))
    declared += [
        _positive("longitudinal_inertia_factor", "1", "c_L: I_L = n c_L (1/12) B L^3"),
        Parameter("longitudinal_center_of_flotation", (1, 1), "m", "LCF x_F from the CO, body x"),
        Parameter("reference_point", (3, 1), "m", "CO -> point P where G is given, body axes (FRD)"),
        _positive("water_density", "kg/m^3", "water density rho"),
    ]
    if gravity_source == "latitude":
        declared.append(Parameter("latitude", (1, 1), "deg", "geodetic latitude (WGS-84 normal gravity)",
                                  -90.0, 90.0))
    else:
        declared.append(_positive("gravity", "m/s^2", "acceleration of gravity g"))
    return tuple(declared)


def surface_restoring_matrices(water_density, gravity, displacement_volume, waterplane_area,
                               transverse_metacentric_height, longitudinal_metacentric_height,
                               longitudinal_center_of_flotation, reference_point):
    """``(G, G_CO, G_CF)`` (6x6 each, CasADi; symbols or numbers in):
    Fossen 2011, eqs. 4.24, p. 65, and 7.250, p. 181; MSS ``Gmtrx.m`` 29-37
    with ``rho`` and ``g`` as inputs (``Gmtrx.m`` 25-26 fixes them).

    ``G_CF`` is diagonal at the centre of flotation ``[LCF 0 0]``, ``G_CO``
    at the CO and ``G`` at ``reference_point`` (from the CO, BODY).
    """
    rho_g = water_density * gravity
    # G_CF (Fossen 2011, eq. 4.24, p. 65; Gmtrx.m 32-35)
    g_cf = ca.diag(ca.vertcat(
        0.0,
        0.0,
        rho_g * waterplane_area,  # (Gmtrx.m 32)
        rho_g * displacement_volume * transverse_metacentric_height,  # (Gmtrx.m 33)
        rho_g * displacement_volume * longitudinal_metacentric_height,  # (Gmtrx.m 34)
        0.0,
    ))
    h_f = MatrixTransforms.H_matrix(ca.vertcat(longitudinal_center_of_flotation, 0.0, 0.0))  # r_bF (Gmtrx.m 29)
    g_co = h_f.T @ g_cf @ h_f  # (Fossen 2011, eq. 7.250, p. 181; Gmtrx.m 36)
    h_p = MatrixTransforms.H_matrix(reference_point)
    g_p = h_p.T @ g_co @ h_p  # (Gmtrx.m 37)
    return g_p, g_co, g_cf


def metacentric_heights(displacement_volume, center_of_buoyancy_above_keel,
                        transverse_waterplane_inertia, longitudinal_waterplane_inertia,
                        center_of_gravity, draft):
    """``[BM_T, BM_L, KM_T, KM_L, GM_T, GM_L]`` in m (6x1, CasADi; symbols or
    numbers in): ``BM = I / nabla`` (Fossen 2011, eq. 4.33, p. 66),
    ``KM = KB + BM``, ``KG = draft - z_g`` (CO on the waterline),
    ``GM = KM - KG`` (Fossen 2011, eq. 4.32, p. 65; MSS ``otter.m`` 179-185).
    """
    transverse_radius = transverse_waterplane_inertia / displacement_volume  # (Fossen 2011, eq. 4.33, p. 66; otter.m 179)
    longitudinal_radius = longitudinal_waterplane_inertia / displacement_volume  # (Fossen 2011, eq. 4.33, p. 66; otter.m 180)
    center_of_gravity_above_keel = draft - center_of_gravity[2]  # KG = T - z_g (otter.m 183; exShipHydrostatics.m 41)
    transverse_metacentre = center_of_buoyancy_above_keel + transverse_radius  # (otter.m 181)
    longitudinal_metacentre = center_of_buoyancy_above_keel + longitudinal_radius  # (otter.m 182)
    return ca.vertcat(
        transverse_radius,
        longitudinal_radius,
        transverse_metacentre,
        longitudinal_metacentre,
        transverse_metacentre - center_of_gravity_above_keel,  # GM_T = KM_T - KG (Fossen 2011, eq. 4.32, p. 65; otter.m 184)
        longitudinal_metacentre - center_of_gravity_above_keel,  # GM_L (Fossen 2011, eq. 4.32, p. 65; otter.m 185)
    )


def surface_hydrostatics_casadi(*, hull_count, gravity_source="latitude",
                                center_of_buoyancy_area="waterplane"):
    """The block: ``(eta, <each declared parameter by name>) -> (g, G, G_CO,
    G_CF, <coefficients>)``, ``g = G eta``.

    Contract
    --------
    Keywords choose the graph (``hull_count`` 1 or 2, required;
    ``gravity_source``; ``center_of_buoyancy_area``; module docstring); an
    unknown value raises ``ValueError``. Inputs: ``eta`` (6x1, ``[x y z phi
    theta psi]``, NED, m and rad) and the parameters of
    ``surface_hydrostatics_parameters(...)`` with the same keywords, each
    under its own name and declared shape. Outputs: the names of
    ``SURFACE_HYDROSTATICS_OUTPUTS``. ``g`` enters the equation of motion on
    the left, so the force on the craft is ``-g``. No value is checked
    inside the graph (``check_values`` of the helper does that).
    """
    declared = surface_hydrostatics_parameters(
        hull_count=hull_count, gravity_source=gravity_source,
        center_of_buoyancy_area=center_of_buoyancy_area,
    )
    p = symbols(declared)
    eta = ca.SX.sym("eta", 6)

    n = float(hull_count)
    length, beam, draft = p["length"], p["beam"], p["draft"]
    nabla, cw = p["displacement_volume"], p["waterplane_area_coefficient"]
    offset = p["hull_lateral_offset"] if hull_count == 2 else 0.0  # one hull on the centreline

    if gravity_source == "latitude":
        gravity = ECEFNEDtransform.gravity(p["latitude"])  # WGS-84 normal gravity (gravity.m 11-12, through more_transformations)
    else:
        gravity = p["gravity"]

    hull_area = cw * length * beam  # A_hull = Cw L B (otter.m 174; exShipHydrostatics.m 32)
    waterplane_area = n * hull_area  # A_wp of n hulls (otter.m 187, "2 * Aw_pont")
    munro_smith = 6.0 * cw**3 / ((1.0 + cw) * (1.0 + 2.0 * cw))  # (exShipHydrostatics.m 38)
    # per hull, plus the parallel-axis term of each hull at +-y (exShipHydrostatics.m 47; otter.m 175-176)
    transverse_inertia = n * length * beam**3 / 12.0 * munro_smith + n * hull_area * offset**2
    # (Fossen 2011, eq. 4.35, p. 66; exShipHydrostatics.m 48; otter.m 177)
    longitudinal_inertia = n * p["longitudinal_inertia_factor"] * beam * length**3 / 12.0
    # Morrish, per hull; A_KB = A_hull by default (exShipHydrostatics.m 35; osv.m 84; otter.m 178),
    # L B with "length_times_beam" (otter.m 177 up to revision 72656d1)
    area_for_kb = hull_area if center_of_buoyancy_area == "waterplane" else length * beam
    center_of_buoyancy_above_keel = (2.5 * draft - nabla / n / area_for_kb) / 3.0  # (exShipHydrostatics.m 35; otter.m 178)

    heights = metacentric_heights(nabla, center_of_buoyancy_above_keel, transverse_inertia,
                                  longitudinal_inertia, p["center_of_gravity"], draft)
    stiffness, stiffness_at_origin, stiffness_at_flotation = surface_restoring_matrices(
        p["water_density"], gravity, nabla, waterplane_area, heights[4], heights[5],
        p["longitudinal_center_of_flotation"], p["reference_point"],
    )
    outputs = {
        "g": stiffness @ eta,  # g = G eta (Fossen 2011, eq. 4.25, p. 65; otter.m 262)
        "G": stiffness,
        "G_CO": stiffness_at_origin,
        "G_CF": stiffness_at_flotation,
        "waterplane_area": waterplane_area,
        "transverse_waterplane_inertia": transverse_inertia,
        "longitudinal_waterplane_inertia": longitudinal_inertia,
        "center_of_buoyancy_above_keel": center_of_buoyancy_above_keel,
        "transverse_metacentric_height": heights[4],
        "longitudinal_metacentric_height": heights[5],
        "acceleration_of_gravity": gravity,
    }
    names = [d.name for d in declared]
    return ca.Function(
        "surface_hydrostatics",
        [eta, *[p[name] for name in names]],
        [outputs[name] for name in SURFACE_HYDROSTATICS_OUTPUTS],
        ["eta", *names],
        list(SURFACE_HYDROSTATICS_OUTPUTS),
    )
