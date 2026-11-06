"""Generator of the Marie references: the second geometry of the generality check (``scaled_torpedo_parameters.json``) and the
floating equilibria of both sets (``*_floating_reference.csv``), from the named part's parameter file
(``marie_parameters.json``, values and their kinds).

Pure Python with numpy and scipy, no CasADi and no ``more_dynamics`` import: a second implementation of the cited equations.
Run ``python3 generate_marie_references.py`` in this folder.

Equations:

* Froude similitude of the second geometry at the scale ``s`` (the shaft speed ``s^-0.5``, the power-law thrust coefficient ``k``
  ``s^(3 + p/2)`` so that ``k n^p`` scales as a force, the torque per thrust as a length; lengths ``s``, masses ``s^3``, inertias ``s^5``, times and
  time constants ``s^0.5``, frequencies ``s^-0.5``, the added mass ``A_ij`` as ``s^(3 + a_i + a_j)`` and the potential
  damping as ``s^(2.5 + a_i + a_j)`` with ``a = 0`` for a translation and ``1`` for a rotation, the stiffness ``G_ij`` as
  ``s^(2 + a_i + a_j)``, the speed factors ``kappa`` (s/rad) ``s^0.5`` and ``k_u`` (s/m) ``s^-0.5``; dimensionless
  coefficients unchanged): the scaled body has the same dimensionless dynamics at equal Froude number.
* The smooth body of revolution the station table samples (class ``Body``): a cylinder with equal ellipsoidal ends, ``r(x) = R
  sqrt(1 - ((|x| - c) / a)^2)`` on an end, scaled so the hull displaces ``displaced_volume``; the station table of the part is
  not read, so the library's table and station sums answer for their discretisation against the smooth hull. At a pose with
  pitch ``theta`` (roll zero) the section at ``x`` has its centre at the depth ``d(x) = (z - sin(theta) (a_x + x) + cos(theta)
  a_z) / cos(theta)`` below the waterline in its plane; its immersed part is the circular segment ``A = r^2 acos(c/r) - c
  sqrt(r^2 - c^2)``, ``c = -d``, with the first moment ``(2/3) (r^2 - c^2)^(3/2)`` (Weisstein, Circular Segment, MathWorld). Volume
  and centre of buoyancy by adaptive quadrature (``scipy.integrate.quad``, absolute tolerance 1e-13 per piece, the joints of the
  cylinder and the ends and the kinks of the segment as break points), not by the station sums of the library.
* Restoring of eq. 4.5 (Fossen 2011, p. 60) with ``B = rho g V(eta)``; the equilibrium is the root of the heave force and
  the pitch moment (``scipy.optimize.least_squares``, the pitch bounded to 1.3 rad: the law holds for
  ``|theta| < 90 deg``).
"""
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

HERE = Path(__file__).resolve().parent
SCALE = 1.6
END_AXIS_PER_LENGTH = 0.35 / 2.15              # idealised hull: equal ellipsoidal ends of 0.35 m on the 2.15 m length (cylinder between)
RESERVES = (0.01, 0.02, 0.03, 0.04)
OFFSETS = (0.0, -0.005, -0.010, -0.020)        # m at scale 1: x_G - x_B
DECLARED_OFFSETS = (-0.010, -0.020)


def values(data, group):
    return {n: e["value"] for n, e in data["parameters"][group].items()}


# ------------------------------------------------------------------ similitude
def _matrix_power(base):
    a = [0, 0, 0, 1, 1, 1]
    return [[base + a[i] + a[j] for j in range(6)] for i in range(6)]


def scale_set(data, s):
    out = json.loads(json.dumps(data))
    out["part"] = "Scaled torpedo (generality check)"
    out["description"] = (f"Not a vehicle: the named part's body scaled by {s} with Froude similitude (module docstring of the "
                          "generator); test data of the generality check only.")
    p = out["parameters"]
    v, h, d = p["vehicle"], p["hydrostatics"], p["hydrodynamics"]
    n = len(v["frequencies"]["value"])
    power_a, power_b = _matrix_power(3), _matrix_power(2.5)

    def table(values_, power):
        flat = np.asarray(values_, float).reshape(n, 6, 6) * (s ** np.asarray(power, float))[None]
        return flat.reshape(-1).tolist()

    def mul(entry, factor):
        entry["value"] = (np.asarray(entry["value"], float) * factor).tolist() if isinstance(entry["value"], list) \
            else entry["value"] * factor

    v["frequencies"]["value"] = [w * s ** -0.5 if w < 10.0 else w for w in v["frequencies"]["value"]]  # the sentinel stays
    v["added_mass"]["value"] = table(v["added_mass"]["value"], power_a)
    v["potential_damping"]["value"] = table(v["potential_damping"]["value"], power_b)
    g = np.asarray(v["hydrostatic_stiffness_345"]["value"], float).reshape(3, 3)
    a = [0, 1, 1]
    v["hydrostatic_stiffness_345"]["value"] = (g * np.array([[s ** (2 + a[i] + a[j]) for j in range(3)] for i in range(3)])
                                               ).reshape(-1).tolist()
    for name, factor in (("database_reference", s), ("length", s), ("beam", s), ("height", s), ("body_mass", s ** 3),
                         ("inertia_diagonal", s ** 5), ("inertia_products", s ** 5), ("body_center_of_gravity", s),
                         ("omega_p", s ** -0.5)):
        mul(v[name], factor)
    for name, factor in (("center_of_buoyancy", s), ("displaced_volume", s ** 3), ("hull_profile_x", s),
                         ("hull_profile_radius", s), ("hull_axis_origin", s), ("hull_profile_cb_tolerance", s)):
        mul(h[name], factor)
    for name, factor in (("T_1236", s ** 0.5), ("kappa_456", s ** 0.5), ("surge_damping_blend", s ** -0.5)):
        mul(d[name], factor)
    q = p["propulsor"]
    exponent = q["power_law_exponent"]["value"]
    for name, factor in (("min_command", s ** -0.5), ("max_command", s ** -0.5), ("position", s),
                         ("power_law_coefficient", s ** (3.0 + 0.5 * exponent)), ("reaction_torque_per_thrust", s)):
        mul(q[name], factor)                      # force s^3 = k n^p with n = s^-0.5; the torque per thrust is a length
    for fin in p["fins"]:
        mul(fin["fin_position"], s)
    for name, factor in (("fin_area", s ** 2), ("fin_span", s)):
        mul(p["fin_common"][name], factor)
    for group in p.values():
        for e in (group.values() if isinstance(group, dict) else []):
            if e["kind"] != "design choice":
                e["kind"] = "estimate"
    return out


# ------------------------------------------------------------------ the floating body, second implementation
class Body:
    """The smooth hull the station table is cut from: a cylinder of radius ``R`` = height / 2 and the length ``L`` of the part, with
    equal ellipsoidal ends of semi-axis ``a = END_AXIS_PER_LENGTH L`` along the hull; ``r(x) = R sqrt(1 - ((|x| - (L/2 - a)) /
    a)^2)`` on an end, the axis point at mid-length. The radius is scaled so the hull displaces ``displaced_volume`` (as the
    library scales its table), the centroid carried to ``center_of_buoyancy``. The table of the part is not read: the
    references are those of the smooth hull, so the library's station sums answer for the discretisation too."""

    def __init__(self, data):
        h, v = values(data, "hydrostatics"), values(data, "vehicle")
        length, radius = float(v["length"]), 0.5 * float(v["height"])
        self.half = 0.5 * length
        self.axis = END_AXIS_PER_LENGTH * length
        self.cylinder_half = self.half - self.axis
        self.a0 = np.asarray(h["hull_axis_origin"], float)
        self.volume_target = float(h["displaced_volume"])
        self.r_cb = np.asarray(h["center_of_buoyancy"], float)
        self.r0 = radius
        self.scale_r = 1.0
        whole = self._integrals(1e6, 0.0)
        self.scale_r = math.sqrt(self.volume_target / whole[0])
        whole = self._integrals(1e6, 0.0)
        self.shift = self.r_cb - np.array([self.a0[0] + whole[1] / whole[0], self.a0[1], self.a0[2]])

    def radius(self, x):
        s = max(abs(x) - self.cylinder_half, 0.0) / self.axis
        return self.scale_r * self.r0 * math.sqrt(max(1.0 - s * s, 0.0))

    def _integrals(self, z, theta):
        """(V, int x A dx, int s A dx) at the depth ``z`` of the origin and the pitch ``theta``."""
        c_t, s_t = math.cos(theta), math.sin(theta)

        def chord(x):
            return -(z - s_t * (self.a0[0] + x) + c_t * self.a0[2]) / c_t

        def area(x):
            r, c = self.radius(x), chord(x)
            if c >= r:
                return 0.0
            if c <= -r:
                return math.pi * r * r
            return r * r * math.acos(c / r) - c * math.sqrt(r * r - c * c)

        def moment(x):
            r, c = self.radius(x), chord(x)
            return 0.0 if abs(c) >= r else (2.0 / 3.0) * (r * r - c * c) ** 1.5

        totals = np.zeros(3)
        edges = [-self.half, -self.cylinder_half, self.cylinder_half, self.half]
        for a, b in zip(edges[:-1], edges[1:]):
            # kinks where the waterline meets the section edge: |chord(x)| = r(x); found by bisection on 64 sub-intervals
            points = []
            grid = np.linspace(a, b, 65)
            for sign in (1.0, -1.0):
                f = [chord(x) - sign * self.radius(x) for x in grid]
                for k in range(64):
                    if f[k] * f[k + 1] < 0.0:
                        points.append(optimize.brentq(lambda x: chord(x) - sign * self.radius(x), grid[k], grid[k + 1], xtol=1e-15))
            for k, f_ in enumerate((area, lambda x: x * area(x), moment)):
                totals[k] += integrate.quad(f_, a, b, points=sorted(points) or None, epsabs=1e-13, epsrel=1e-12, limit=400)[0]
        return totals

    def buoyancy_and_centre(self, z, theta, rho, g):
        volume, first_x, first_s = self._integrals(z, theta)
        centre = np.array([self.a0[0] + first_x / volume, self.a0[1], self.a0[2] + first_s / volume]) + self.shift
        return rho * g * volume, volume, centre


def restoring(eta, weight, buoyancy, r_g, r_b):
    """g(eta) of a body with W at r_g and B at r_b (Fossen 2011, eq. 4.5, p. 60), in Euler angles."""
    phi, theta = eta[3], eta[4]
    st, ct, sp, cp = math.sin(theta), math.cos(theta), math.sin(phi), math.cos(phi)
    return np.array([(weight - buoyancy) * st, -(weight - buoyancy) * ct * sp, -(weight - buoyancy) * ct * cp,
                     -(r_g[1] * weight - r_b[1] * buoyancy) * ct * cp + (r_g[2] * weight - r_b[2] * buoyancy) * ct * sp,
                     (r_g[2] * weight - r_b[2] * buoyancy) * st + (r_g[0] * weight - r_b[0] * buoyancy) * ct * cp,
                     -(r_g[0] * weight - r_b[0] * buoyancy) * ct * sp - (r_g[1] * weight - r_b[1] * buoyancy) * st])


def equilibrium(body, mass, r_g, rho, g, guess, depth_bound):
    weight = mass * g

    def residual(v):
        b, _, centre = body.buoyancy_and_centre(v[0], v[1], rho, g)
        f = restoring(np.array([0, 0, v[0], 0, v[1], 0.0]), weight, b, r_g, centre)
        return [f[2] / weight, f[4] / weight]

    best = None
    for start in (guess, [guess[0], 0.3], [guess[0], 0.8]):
        fit = optimize.least_squares(residual, start, bounds=([-depth_bound, -1.3], [depth_bound, 1.3]), xtol=1e-15, ftol=1e-15, gtol=1e-15)
        size = float(np.max(np.abs(fit.fun)))
        if best is None or size < best[1]:
            best = (fit.x, size)
    (z, theta), size = best
    if size > 1e-11:
        raise RuntimeError(f"no floating equilibrium with |pitch| < 1.3 rad (residual {size:.2e} of the weight)")
    b, volume, centre = body.buoyancy_and_centre(z, theta, rho, g)
    return z, theta, volume, centre


def floating_rows(data, scale):
    body = Body(data)
    site, vehicle = values(data, "site"), values(data, "vehicle")
    rho, g = site["water_density"], site["gravity"]
    mass, z_g = vehicle["body_mass"], vehicle["body_center_of_gravity"][2]
    x_b = body.r_cb[0]
    radius = values(data, "vehicle")["height"] / 2.0
    cases = [("declared", mass, 0.0)] + [(f"declared_dx_{d * 1000:+.0f}mm", mass, d * scale) for d in DECLARED_OFFSETS]
    for reserve in RESERVES:
        reserve_mass = rho * body.volume_target / (1.0 + reserve)
        cases += [(f"reserve_{reserve:.2f}_dx_{d * 1000:+.0f}mm", reserve_mass, d * scale) for d in OFFSETS]
    rows = []
    for name, m, dx in cases:
        r_g = np.array([x_b + dx, 0.0, z_g])
        z, theta, volume, centre = equilibrium(body, m, r_g, rho, g, [0.4 * radius, 0.0], vehicle["length"])
        rows.append([name, m, x_b + dx, rho * body.volume_target / m - 1.0, z, theta, volume, centre[0], centre[2]])
    return rows


def scaled_speed_reference(s):
    """``scaled_torpedo_surge_speed_reference.csv`` from ``marie_surge_speed_reference.csv``: shaft speed ``s^-0.5``, speed ``s^0.5``, thrust
    ``s^3`` (Froude similitude)."""
    lines = (HERE / "marie_surge_speed_reference.csv").read_text().splitlines()
    rows = [[float(c) for c in line.split(",")] for line in lines[1:]]
    out = [lines[0]] + [f"{r[0] * s ** -0.5!r},{r[1] * s ** 0.5!r},{r[2] * s ** 3!r}" for r in rows]
    (HERE / "scaled_torpedo_surge_speed_reference.csv").write_text("\n".join(out) + "\n")


def main():
    marie = json.loads((HERE / "marie_parameters.json").read_text())
    scaled = scale_set(marie, SCALE)
    scaled_speed_reference(SCALE)
    (HERE / "scaled_torpedo_parameters.json").write_text(json.dumps(scaled, indent=1) + "\n")
    header = ["case", "mass_kg", "cg_x_m", "reserve", "origin_depth_m", "pitch_rad", "submerged_volume_m3", "buoyancy_x_m",
              "buoyancy_z_down_m"]
    for stem, data, scale in (("marie", marie, 1.0), ("scaled_torpedo", scaled, SCALE)):
        with open(HERE / f"{stem}_floating_reference.csv", "w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(header)
            for row in floating_rows(data, scale):
                writer.writerow([row[0]] + [repr(float(v)) for v in row[1:]])


if __name__ == "__main__":
    main()
