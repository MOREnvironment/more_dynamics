"""Origin-invariance and sinkage tests of the payload trim (``payload_trim_form``), PT1 (E-211).

In the linear restoring ``g = G (eta - eta_0)`` the plane ``z = 0`` is the unloaded waterline (MSS ``otter.m`` header:
draft T = 13.4 cm when mP = 0 "corresponds to z = 0"; the payload is "an external force g_0"; Fossen 2011, Handbook of
Marine Craft Hydrodynamics and Motion Control, 1st ed., Wiley, sec. 4.2 and 4.4). A payload mP above the centre of
flotation (CF) therefore sinks the craft level by mP / (rho Awp) and trims it not at all; a payload off the CF adds the
couple of its weight about the CF. The equilibrium (sinkage of the CF, pitch, roll) is a property of the craft and the
physical payload position, so it cannot change when the body origin is moved.

The tests move the body origin by ``shift`` along x (the centre of flotation, the payload and the hull centre of
gravity move by ``-shift`` in the body frame) and read ``eta_0`` of the plugin. Hydrostatic expectations are written out
here from ``rho``, ``g``, ``Awp`` and the lever about the CF, not from the plugin's matrix.

Author:    Enio Krizman
Date:      2026-11-02
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "vehicle_models"))
import vehicle_contract as vc  # noqa: E402

FORMS = ("weight_wrench", "lever_arm_and_heave_shift")
RHO, GRAVITY = 1025.0, 9.81
AWP = 2 * 0.75 * 2.0 * 0.25      # n C_w L B_hull of the Otter pontoons (otter.m 187, 190)
LCF = -0.2                       # centre of flotation of the frozen set, body x (otter.m 192)


def _eta_0(form, mass, position, shift=0.0):
    """``eta_0`` of the Otter with the body origin moved by ``shift`` along x; ``position`` is the payload position in
    the frame of the frozen set (origin unmoved)."""
    p = vc.plugins()
    hydrostatics = (p.SurfaceRestoring, {**vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2,
                                         "payload_trim_form": form, "longitudinal_center_of_flotation": LCF - shift})
    cg = list(vc.OTTER["hull_center_of_gravity"])
    cg[0] -= shift
    position = [position[0] - shift, position[1], position[2]]
    vehicle = vc.otter(hydrostatics=hydrostatics,
                       vehicle_params={"payload_mass": float(mass), "payload_position": position,
                                       "hull_center_of_gravity": cg})
    signals = vehicle.signals(np.zeros(12), [])
    return np.asarray(signals["payload_trim_offset"]).ravel(), shift


def _attitude(form, mass, position, shift=0.0):
    """(sinkage of the centre of flotation, roll, pitch) of the equilibrium, origin-independent quantities."""
    eta_0, shift = _eta_0(form, mass, position, shift)
    lcf = LCF - shift
    return np.array([eta_0[2] - lcf * eta_0[4], eta_0[3], eta_0[4]])   # z_CF = z - x_CF theta (small angles)


@pytest.mark.parametrize("form", FORMS)
def test_a_payload_over_the_centre_of_flotation_sinks_the_craft_level(form):
    """25 kg at the CF (any height): heave of the CF is m_P / (rho Awp) = 3.25 cm, pitch zero, no roll."""
    sinkage, roll, pitch = _attitude(form, 25.0, [LCF, 0.0, -0.35])
    assert sinkage == pytest.approx(25.0 / (RHO * AWP), rel=1e-9)
    assert abs(pitch) <= 1e-12 and abs(roll) <= 1e-12


@pytest.mark.parametrize("form", FORMS)
@pytest.mark.parametrize("position", ([-0.2, 0.0, -0.35], [0.05, 0.0, -0.35], [0.3, 0.2, -0.1], [-0.7, -0.15, -0.35]))
def test_the_equilibrium_does_not_depend_on_where_the_body_origin_is(form, position):
    """The same physical payload, the body origin moved by -0.5 .. +0.8 m: same sinkage of the CF, roll and pitch."""
    reference = _attitude(form, 25.0, position, 0.0)
    for shift in (-0.5, 0.2, 0.8):
        assert np.allclose(_attitude(form, 25.0, position, shift), reference, rtol=1e-9, atol=1e-12), (shift, position)


@pytest.mark.parametrize("form", FORMS)
def test_the_sinkage_of_the_centre_of_flotation_is_the_hydrostatic_one_for_any_payload_position(form):
    """Linear hydrostatics: the CF sinks m_P / (rho Awp) wherever the payload is, and the pitch is the couple of the
    weight about the CF over the pitch stiffness (the weight's lever is x_P - x_CF)."""
    sinkage, _, pitch = _attitude(form, 25.0, [0.3, 0.0, -0.1])
    assert sinkage == pytest.approx(25.0 / (RHO * AWP), rel=1e-9)
    assert pitch < 0.0   # a weight forward of the CF trims the bow down (theta positive = bow up)


def test_the_default_form_is_the_physics_one():
    p = vc.plugins()
    vehicle = vc.otter(hydrostatics=(p.SurfaceRestoring, {**vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2}),
                       vehicle_params={"payload_mass": 25.0, "payload_position": [LCF, 0.0, -0.35]})
    eta_0 = np.asarray(vehicle.signals(np.zeros(12), [])["payload_trim_offset"]).ravel()
    assert eta_0[2] - LCF * eta_0[4] == pytest.approx(25.0 / (RHO * AWP), rel=1e-9) and abs(eta_0[4]) <= 1e-12


def test_G4_the_legacy_form_is_told_by_these_gates():
    """The gates can fail: the legacy form (no heave term, moment about the origin) leaves the craft level with a payload
    over the centre of flotation unsunk and pitched bow up (GC2 section E: +1.2 deg for 25 kg), and its pitch moves
    with the body origin."""
    sinkage, _, pitch = _attitude("lever_arm_moment", 25.0, [LCF, 0.0, -0.35])
    assert sinkage == pytest.approx(0.0, abs=1e-12) and np.degrees(pitch) > 1.0
    here = _attitude("lever_arm_moment", 25.0, [0.05, 0.0, -0.35], 0.0)
    moved = _attitude("lever_arm_moment", 25.0, [0.05, 0.0, -0.35], 0.8)
    assert abs(here[2] - moved[2]) > 1e-2
