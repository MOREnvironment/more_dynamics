"""Otter reference gates for ``Catamaran``.

The hull blocks are read through the vehicle's ``signals`` (every named
quantity of the vehicle and of its children's outputs); no actuator is
attached, so the vehicle takes no input, apart from the three current values
in the rigid-body reference, which was run in a current.

Author:    Enio Krizman
Date:      2026-10-08
"""

from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc

TESTS = Path(__file__).resolve().parents[1]
G1_TOLERANCE = 1e-9
CURRENT_SPEED, CURRENT_DIRECTION = None, None  # read from the rigid-body test, which holds the reference's current


def _table(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: np.atleast_1d(data[name]) for name in data.dtype.names}


def _cols(table, prefix, count):
    return np.column_stack([table[f"{prefix}_{i:02d}"] for i in range(1, count + 1)])


def _reference_current():
    import math
    import re
    text = (TESTS / "rigid_body" / "test_rigid_body_block.py").read_text()
    speed = float(re.search(r"^CURRENT_SPEED = ([0-9.]+)$", text, re.M).group(1))
    direction = math.radians(float(re.search(r"^CURRENT_DIRECTION = np.deg2rad\(([0-9.]+)\)$", text, re.M).group(1)))
    return [speed, direction, 0.0]


def test_catamaran_is_one_plugin_and_uses_the_cross_flow_strip_of_the_torpedo_with_its_own_section():
    """The two surface hulls and the torpedo share the strip integral; the section law differs by plugin."""
    from more_dynamics.models.cross_flow.cross_flow_strip import (
        cross_flow_strip_circular_cylinder_reynolds, cross_flow_strip_rectangular_section_hoerner)
    assert set(vc.plugins().Catamaran.COMPONENTS) == {"hydrostatics", "hydrodynamics", "actuators"}
    assert cross_flow_strip_circular_cylinder_reynolds().name_in() != cross_flow_strip_rectangular_section_hoerner().name_in()
    declared = {d.name for d in vc.plugins().Catamaran.PARAMETERS}
    assert set(vc.OTTER) - {"longitudinal_inertia_factor", "longitudinal_center_of_flotation", "max_forward_thrust",
                            "max_speed", "time_constants", "damping_ratios", "yaw_damping_nonlinearity"} <= declared


def test_every_default_equals_the_frozen_otter_value():
    p = vc.plugins()
    for cls in (p.Catamaran, p.SurfaceRestoring, p.SurfaceHullLoads):
        for d in cls.PARAMETERS:
            if d.name in vc.OTTER:
                assert np.array_equal(np.asarray(d.default_value, dtype=float),
                                      np.asarray(vc.OTTER[d.name], dtype=float)), (cls.__name__, d.name)


def test_G1_catamaran_matches_existing_otter_block_references():
    vehicle = vc.otter()
    current_vehicle = vc.otter(current="full_rotation_rate")  # the rigid-body reference was run in a current
    assert vc.input_names(vehicle) == []
    assert vc.input_names(current_vehicle) == list(vc.CURRENT_NAMES)
    current = np.array(_reference_current())
    rigid_in = np.loadtxt(TESTS / "rigid_body" / "data" / "inputs.csv", delimiter=",", skiprows=1)
    rigid_ref = _table(TESTS / "rigid_body" / "data" / "matlab_reference_mss_current.csv")
    hydro_ref = _table(TESTS / "hydrodynamic_loads" / "data" / "catamaran_matlab_reference_mss_cc07579.csv")
    worst = {}

    for k, state in enumerate(rigid_in[:, :12]):
        out = current_vehicle.signals(vc.mss_to_ours(state), current)
        for key, name in (("M_RB", "rigid_body_mass_matrix"), ("M_A", "added_mass_matrix"),
                          ("C_RB", "rigid_body_coriolis_matrix"), ("C_A", "added_mass_coriolis_matrix")):
            ref = _cols(rigid_ref, key, 36)[k].reshape(6, 6)
            worst[key] = max(worst.get(key, 0.0), np.abs(vc.matrix(out[name]) - ref).max())
        ref = _cols(rigid_ref, "C_total", 36)[k].reshape(6, 6)
        total = vc.matrix(out["rigid_body_coriolis_matrix"]) + vc.matrix(out["added_mass_coriolis_matrix"])
        worst["C_total"] = max(worst.get("C_total", 0.0), np.abs(total - ref).max())

    nothing = np.zeros(0)
    for k, nu_r in enumerate(_cols(hydro_ref, "nu_r", 6)):
        out = vehicle.signals(np.r_[np.zeros(6), nu_r], nothing)
        checks = {
            "M": (out["mass_matrix"], _cols(hydro_ref, "M_total", 36)[k]),
            "G": (out["restoring_matrix"], _cols(hydro_ref, "G", 36)[k]),
            "D": (out["damping_matrix"], -_cols(hydro_ref, "D", 36)[k]),
            "damping_wrench": (out["hydrodynamics.damping_force"], _cols(hydro_ref, "tau_damp", 6)[k]),
            "cross_flow_wrench": (out["hydrodynamics.cross_flow_force"], _cols(hydro_ref, "tau_crossflow", 6)[k]),
        }
        for name, (actual, reference) in checks.items():
            if name in ("M", "G", "D"):  # the references are row-major 6x6; the payload matrix is column-major
                reference = reference.reshape(6, 6).ravel(order="F")
            worst[name] = max(worst.get(name, 0.0), np.abs(np.asarray(actual).ravel() - reference).max())
    assert len(rigid_in) == len(_cols(hydro_ref, "nu_r", 6)) == 50
    assert max(worst.values()) <= G1_TOLERANCE, worst
