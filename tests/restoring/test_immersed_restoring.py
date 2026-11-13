"""The restoring of a body of revolution through the free surface (``models/restoring/immersed.py``): the circle segment
against its closed forms, the station sums against exact volumes, the fully submerged limit against the submerged block
(Fossen 2011, eq. 4.5, p. 60), finite derivatives at dry and whole sections, and the plugin's refusal of a table whose
centroid is not the centre of buoyancy.

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Eqs. 4.1, 4.5, pp. 59-60.
Weisstein, E. W. Circular Segment. MathWorld (half disc: area pi r^2 / 2, centroid 4 r / (3 pi) from the centre).

Author:    Enio Krizman
Date:      2026-11-05
"""
import math
import sys
from pathlib import Path

import casadi as ca
import numpy as np
import pytest

from more_dynamics.models.restoring.immersed import circle_segment, immersed_restoring
from more_dynamics.models.restoring.submerged import submerged_hydrostatics_casadi

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "vehicle_models"))
import vehicle_contract as vc  # noqa: E402

R = 0.2


def _segment(depth, radius=R):
    area, offset = circle_segment(ca.DM(radius), ca.DM(depth))
    return float(area), float(offset)


def test_the_circle_segment_has_its_closed_forms():
    assert _segment(0.0) == pytest.approx((math.pi * R ** 2 / 2, 4 * R / (3 * math.pi)), rel=1e-14)   # half disc
    assert _segment(R) == (pytest.approx(math.pi * R ** 2, rel=1e-15), 0.0)                         # whole
    assert _segment(5 * R) == (pytest.approx(math.pi * R ** 2, rel=1e-15), 0.0)
    assert _segment(-R) == (0.0, 0.0) and _segment(-5 * R) == (0.0, 0.0)                             # dry
    assert _segment(0.0, 0.0) == (0.0, 0.0)                                                          # no section
    a1, _ = _segment(0.3 * R)
    a2, _ = _segment(-0.3 * R)
    assert a1 + a2 == pytest.approx(math.pi * R ** 2, rel=1e-14)                                     # complements


def test_the_circle_segment_has_finite_derivatives_at_every_depth():
    d = ca.SX.sym("d")
    area, offset = circle_segment(ca.SX(R), d)
    jacobian = ca.Function("j", [d], [ca.jacobian(ca.vertcat(area, offset), d)])
    for depth in (-2 * R, -R, -R + 1e-12, 0.0, R - 1e-12, R, 2 * R):
        assert np.all(np.isfinite(np.asarray(jacobian(depth))))


def _cylinder(length=2.0, slices=40):
    x = [-length / 2, length / 2]
    return immersed_restoring(x, slices), x


def _call(block, eta, radius=R, volume=None, r_cb=(0.0, 0.0, 0.0), mass=None, r_g=(0.0, 0.0, 0.01)):
    volume = math.pi * radius ** 2 * 2.0 if volume is None else volume
    mass = 1025.0 * volume / 1.01 if mass is None else mass
    return block(eta=eta, mass=mass, gravity=9.81, water_density=1025.0, center_of_gravity=r_g, center_of_buoyancy=r_cb,
                 displaced_volume=volume, hull_profile_radius=[radius, radius], hull_axis_origin=[0.0, 0.0, 0.0])


def test_a_level_cylinder_floating_at_its_axis_displaces_half_its_volume_with_the_half_disc_centroid():
    block, _ = _cylinder()
    out = _call(block, [0, 0, 0.0, 0, 0, 0])
    assert float(out["immersed_volume"]) == pytest.approx(math.pi * R ** 2, rel=1e-13)
    assert np.asarray(out["immersed_center_of_buoyancy"]).ravel() == pytest.approx([0.0, 0.0, 4 * R / (3 * math.pi)],
                                                                                    abs=1e-14)


@pytest.mark.parametrize("eta", [[1, 2, 3.0, 0.3, -0.4, 1.0], [0, 0, 1.5, -1.0, 0.8, -2.0]])
def test_fully_submerged_it_is_the_submerged_block_exactly(eta):
    block, _ = _cylinder()
    volume = math.pi * R ** 2 * 2.0
    out = _call(block, eta)
    weight = 1025.0 * volume / 1.01 * 9.81
    expected = submerged_hydrostatics_casadi()(eta=eta, weight=weight, buoyancy=1025.0 * 9.81 * volume,
                                               center_of_gravity=[0, 0, 0.01], center_of_buoyancy=[0, 0, 0])["g"]
    assert np.asarray(out["g"]).ravel() == pytest.approx(np.asarray(expected).ravel(), abs=1e-12 * weight)


def test_the_radii_are_scaled_to_the_displaced_volume_and_the_centroid_shift_is_carried():
    """A cone (radius linear from 0 at x = -1 m to R at x = +1 m): its centroid is 3/4 of the length from the apex,
    x = +0.5 m; the radii are scaled to any given volume, and the cone floating level at its axis displaces half."""
    block = immersed_restoring([-1.0, 1.0], 2000)
    cone = dict(gravity=9.81, water_density=1025.0, center_of_gravity=[0, 0, 0], mass=1.0, hull_axis_origin=[0, 0, 0],
                hull_profile_radius=[0.0, R])
    exact = math.pi * R ** 2 * 2.0 / 3.0
    out = block(eta=[0, 0, 10.0, 0, 0, 0], center_of_buoyancy=[0.5, 0, 0], displaced_volume=exact, **cone)
    assert float(out["profile_centroid_offset"]) == pytest.approx(0.0, abs=1e-6)    # cone centroid 3/4 from the apex
    half = block(eta=[0, 0, 0.0, 0, 0, 0], center_of_buoyancy=[0.5, 0, 0], displaced_volume=2 * exact, **cone)
    assert float(half["immersed_volume"]) == pytest.approx(exact, rel=1e-5)         # level at the axis: half of 2 V


def test_the_plugin_refuses_a_table_whose_centroid_is_not_the_centre_of_buoyancy():
    p = vc.plugins()
    values = {"buoyancy_method": "immersed_body_of_revolution", "center_of_buoyancy": [0.3, 0.0, 0.0],
              "displaced_volume": 0.05, "hull_profile_x": [-1.0, 1.0], "hull_profile_radius": [0.1, 0.1],
              "hull_axis_origin": [0.0, 0.0, 0.0], "hull_profile_slices": 10, "hull_profile_cb_tolerance": 1e-3}
    context = vc.context_of(p.SubmergedRestoring, values, reference=False)
    with vc.refusing(ValueError, match="centroid of the hull profile"):
        vc.configure(context)
    context = vc.context_of(p.SubmergedRestoring, {**values, "center_of_buoyancy": [0.0, 0.0, 0.0]}, reference=False)
    vc.configure(context)
