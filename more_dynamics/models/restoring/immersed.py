"""Restoring of a body of revolution that can reach the free surface: the buoyancy and its centre follow the part of
the hull under the calm free surface, section by section, as one CasADi block.

Frames: the hull is a body of revolution about a line parallel to BODY x (forward) through ``hull_axis_origin`` (BODY
FRD, from the body origin CO); the section at station ``x`` has its centre at ``p(x) = hull_axis_origin + [x, 0, 0]``
and radius ``r(x)``, linear between the stations of the table ``(hull_profile_x, hull_profile_radius)``. NED, z down;
the free surface is the plane ``z = 0``. ``eta = [x y z phi theta psi]`` with ZYX Euler angles.

The depth below the surface of a BODY point ``q`` is ``D(q) = z + R(3,:) q`` (``R(3,:) = [-s theta, c theta s phi,
c theta c phi]``, Fossen 2011, eq. 2.18, p. 22). In the plane of a section the depth grows along ``n = (s phi, c phi)``
(BODY y, z) at the rate ``cos theta``, so the waterline is a chord at the signed distance ``c = -D(p(x)) / cos theta``
from the centre (``|theta| < 90 deg``). The part of a circle of radius ``r`` below a chord at distance ``c`` (``|c| <=
r``) has the area and centroid of the circular segment::

    A(r, c) = r^2 acos(c / r) - c sqrt(r^2 - c^2),       s(r, c) = (2/3) (r^2 - c^2)^(3/2) / A   (along n)

with ``A = pi r^2, s = 0`` for a section under water and ``A = 0`` for a dry one (Weisstein, Circular Segment). The
displaced volume and the centre of buoyancy are the sums over the sections (trapezoid weights on the refined stations,
``hull_profile_slices`` equal pieces per interval)::

    V(eta) = sum A w,        r_bb(eta) = (1/V) sum A w (p(x) + s [0, n])                (Fossen 2011, eq. 4.1, p. 59)

and ``g(eta)`` is eq. 4.5 (Fossen 2011, p. 60) with ``W = m g``, ``B = rho g V(eta)`` at ``r_bb(eta)``. When every section
is whole under water the law is exactly the submerged one: ``B = rho g displaced_volume`` at ``center_of_buoyancy``.

Configure-time scaling: the radii are scaled by ``sqrt(displaced_volume / V_table)`` so the whole hull displaces
``displaced_volume``; the centroid of the table is carried with a constant shift to ``center_of_buoyancy`` (one truth for
the centre of buoyancy). The output ``profile_centroid_offset`` (m) is the distance a plugin checks against
``hull_profile_cb_tolerance``.

Guarded expressions: at a dry or a whole section the segment formulas are not evaluated with ``|c| = r`` (the derivative
of ``acos`` and of the square root is infinite there); the centroid of an empty section is zero, not ``0/0``. Every
expression and its derivative is finite at every pose with ``|theta| < 90 deg``.

Not in this law: appendages, the dynamic pressure near the surface, waves, ``|theta| >= 90 deg``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Eq. 2.18, p. 22;
    eq. 4.1, p. 59; eq. 4.5, p. 60.
[Weisstein] Weisstein, E. W. Circular Segment. MathWorld (area and first moment of a circular segment).

Author:    Enio Krizman
Date:      2026-11-05
"""
import math

import casadi as ca

from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter

from more_dynamics.models.restoring.submerged import submerged_hydrostatics_casadi
from more_dynamics.models.shared.wiring import function_from


def immersed_restoring_parameters(station_count):
    """Declared parameters of the block for a table of ``station_count`` stations."""
    return (
        Parameter("center_of_buoyancy", (3, 1), "m", "CO -> CB r_bb of the whole hull, body axes (FRD)"),
        Parameter("displaced_volume", (1, 1), "m^3", "displaced volume of the whole hull", 0.0, minimum_exclusive=True),
        Parameter("hull_profile_radius", (int(station_count), 1), "m",
                  "section radius at each station (scaled so the hull displaces displaced_volume)", 0.0),
        Parameter("hull_axis_origin", (3, 1), "m", "axis point at station x = 0, BODY FRD from CO"),
    )


def refine_stations(hull_profile_x, hull_profile_slices):
    """``(x_fine, [(k, s)])``: the stations cut into ``hull_profile_slices`` equal pieces per interval (the last station
    kept), each fine station with its interval ``k`` and fraction ``s`` (the radius is ``(1-s) r_k + s r_k+1``)."""
    x = [float(v) for v in hull_profile_x]
    n = int(hull_profile_slices)
    if len(x) < 2 or any(b <= a for a, b in zip(x[:-1], x[1:])) or n < 1:
        raise ValueError("hull_profile_x must hold at least two increasing stations and hull_profile_slices be >= 1")
    fine, where = [], []
    for k, (a, b) in enumerate(zip(x[:-1], x[1:])):
        step = (b - a) / n
        for i in range(n):
            value = a + i * step
            fine.append(value)
            where.append((k, (value - a) / (b - a)))
    fine.append(x[-1])
    where.append((len(x) - 2, 1.0))
    return fine, where


def _weights(x):
    """Trapezoid weights of an increasing grid."""
    w = [0.0] * len(x)
    for i in range(len(x) - 1):
        dx = x[i + 1] - x[i]
        w[i] += 0.5 * dx
        w[i + 1] += 0.5 * dx
    return w


def circle_segment(radius, centre_depth):
    """``(area, offset)`` of the immersed part of a circle of ``radius`` whose centre lies ``centre_depth`` below the
    waterline in its plane (negative: above it): the circular segment (module docstring), guarded at ``|c| = r``."""
    c = ca.fmin(ca.fmax(-centre_depth, -radius), radius)  # chord distance from the centre, away from the water
    inside = ca.logic_and(c > -radius, c < radius)
    c_in = ca.if_else(inside, c, 0.0)
    r_in = ca.if_else(radius > 0.0, radius, 1.0)
    root = ca.sqrt(ca.if_else(inside, radius ** 2 - c_in ** 2, 0.0) + ca.if_else(inside, 0.0, 1.0))
    segment = radius * radius * ca.acos(c_in / r_in) - c_in * root  # A = r^2 acos(c/r) - c sqrt(r^2 - c^2) (Weisstein)
    area = ca.if_else(c >= radius, 0.0, ca.if_else(c <= -radius, math.pi * radius ** 2, segment))
    wet = ca.logic_and(inside, segment > 1e-300)
    offset = ca.if_else(wet, (2.0 / 3.0) * root ** 3 / ca.if_else(wet, segment, 1.0), 0.0)  # s = (2/3) root^3 / A
    return area, offset


def immersed_restoring(hull_profile_x, hull_profile_slices):
    """The block ``(eta, mass, gravity, water_density, center_of_gravity, center_of_buoyancy, displaced_volume,
    hull_profile_radius, hull_axis_origin) -> (g, weight, buoyancy, center_of_buoyancy, immersed_volume,
    immersed_center_of_buoyancy, profile_centroid_offset)`` for the station positions ``hull_profile_x`` (m, increasing,
    from the axis point) and ``hull_profile_slices`` pieces per interval (module docstring)."""
    stations = [float(v) for v in hull_profile_x]
    fine, where = refine_stations(stations, hull_profile_slices)
    weights = _weights(fine)
    s = {name: ca.SX.sym(name, *shape) for name, shape in (
        ("eta", (6, 1)), ("mass", (1, 1)), ("gravity", (1, 1)), ("water_density", (1, 1)),
        ("center_of_gravity", (3, 1)), ("center_of_buoyancy", (3, 1)), ("displaced_volume", (1, 1)),
        ("hull_profile_radius", (len(stations), 1)), ("hull_axis_origin", (3, 1)))}
    table, a0 = s["hull_profile_radius"], s["hull_axis_origin"]
    radius = [(1.0 - f) * table[k] + f * table[k + 1] for k, f in where]  # linear between the stations
    # whole hull: volume and centroid of the table, then the radii scaled to the displaced volume
    whole = [math.pi * r ** 2 * w for r, w in zip(radius, weights)]
    table_volume = sum(whole, ca.SX(0.0))
    centroid = ca.vertcat(a0[0] + sum((x * a for x, a in zip(fine, whole)), ca.SX(0.0)) / table_volume, a0[1], a0[2])
    scale = ca.sqrt(s["displaced_volume"] / table_volume)
    radius = [scale * r for r in radius]
    shift = s["center_of_buoyancy"] - centroid  # the table's centroid carried to r_bb (one centre of buoyancy)

    eta = s["eta"]
    row3 = MatrixTransforms.Rzyx_row3(eta[3], eta[4])  # R(3,:) (Fossen 2011, eq. 2.18, p. 22)
    cos_theta = ca.sqrt(row3[1] ** 2 + row3[2] ** 2)
    n_y, n_z = row3[1] / cos_theta, row3[2] / cos_theta  # depth direction in the section plane
    depth_at = [(eta[2] + row3[0] * (a0[0] + x) + row3[1] * a0[1] + row3[2] * a0[2]) / cos_theta for x in fine]
    areas, offsets = zip(*(circle_segment(r, d) for r, d in zip(radius, depth_at)))
    wa = [a * w for a, w in zip(areas, weights)]
    volume = sum(wa, ca.SX(0.0))  # V(eta) = sum A w
    safe_volume = ca.if_else(volume > 0.0, volume, 1.0)
    sb = sum((o * a for o, a in zip(offsets, wa)), ca.SX(0.0)) / safe_volume
    r_bb = ca.vertcat(a0[0] + sum((x * a for x, a in zip(fine, wa)), ca.SX(0.0)) / safe_volume,
                      a0[1] + sb * n_y, a0[2] + sb * n_z)
    r_bb = ca.if_else(volume > 0.0, r_bb + shift, ca.SX.zeros(3))
    submerged = ca.mmin(ca.vertcat(*[d - r for d, r in zip(depth_at, radius)])) >= 0.0  # every section whole

    weight = s["mass"] * s["gravity"]  # W = m g (Fossen 2011, eq. 4.1, p. 59)
    full_buoyancy = s["water_density"] * s["gravity"] * s["displaced_volume"]  # B = rho g nabla (eq. 4.1, p. 59)
    buoyancy = ca.if_else(submerged, full_buoyancy, s["water_density"] * s["gravity"] * volume)
    centre = ca.if_else(submerged, s["center_of_buoyancy"], r_bb)
    g = submerged_hydrostatics_casadi()(eta=eta, weight=weight, buoyancy=buoyancy,
                                         center_of_gravity=s["center_of_gravity"], center_of_buoyancy=centre)["g"]
    return function_from("immersed_restoring", s, {
        "g": g, "weight": weight, "buoyancy": buoyancy, "center_of_buoyancy": s["center_of_buoyancy"],
        "immersed_volume": ca.if_else(submerged, s["displaced_volume"], volume),
        "immersed_center_of_buoyancy": centre, "profile_centroid_offset": ca.norm_2(centroid - s["center_of_buoyancy"])})
