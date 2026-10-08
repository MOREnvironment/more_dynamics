"""Otter reference gates for the catamaran tree of the same ``MarineCraft6DOF``.

The hull blocks are read through the vehicle's diagnostic outputs (the
payload entries a tree lists in ``diagnostic_outputs``); no force producer is
attached, so the vehicle takes no input.

Author:    Enio Krizman
Date:      2026-10-08
"""

import json

import numpy as np

import vehicle_contract as vc
from vehicle_contract import DATA, make_trees

ROOT = DATA.parent
G1_TOLERANCE = 1e-9


def _table(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: np.atleast_1d(data[name]) for name in data.dtype.names}


def _cols(table, prefix, count):
    return np.column_stack([table[f"{prefix}_{i:02d}"] for i in range(1, count + 1)])


def _vec(matrix):
    return np.asarray(matrix).reshape(6, 6).ravel(order="F")  # a payload matrix is its column-major vec


def test_catamaran_is_data_and_uses_the_same_cross_flow_strip_part():
    otter = make_trees.catamaran()
    torpedo = make_trees.torpedo()
    otter_part = otter["children"]["hydrodynamic_loads"][1]
    remus_part = torpedo["children"]["hydrodynamic_loads"][2]
    assert otter_part["plugin"] == remus_part["plugin"] == "CrossFlowStrip"
    assert otter_part["children"]["section"]["plugin"] != remus_part["children"]["section"]["plugin"]
    values = json.loads((DATA / "otter" / "otter_mss_cc07579.json").read_text())
    for slot, item in otter["children"].items():
        for entry in item if isinstance(item, list) else [item]:
            assert set(entry["params"]) <= set(values) | {"kinematic_viscosity", "reference_point",
                                                         "integration_max_step", "open_inputs",
                                                         "diagnostic_outputs", "current_speed",
                                                         "current_direction", "current_vertical_speed"}, slot


def test_G1_catamaran_matches_existing_otter_block_references():
    vehicle = vc.build("otter")
    current_vehicle = vc.build("otter_current")  # the rigid-body reference was run in a current
    assert vc.input_names(vehicle) == []
    rigid_in = np.loadtxt(ROOT / "rigid_body" / "inputs.csv", delimiter=",", skiprows=1)
    rigid_ref = _table(ROOT / "rigid_body" / "matlab_reference_mss_current.csv")
    hydro_ref = _table(ROOT / "hydrodynamics" / "catamaran_matlab_reference_mss_cc07579.csv")
    worst = {}
    nothing = np.zeros(0)

    for k, state in enumerate(rigid_in[:, :12]):
        out = vc.output_entries(current_vehicle, vc.mss_to_ours(state), nothing)
        for key, name in (("M_RB", "rigid_body_mass_matrix"), ("M_A", "added_mass_matrix"),
                          ("C_RB", "rigid_body_coriolis_matrix"), ("C_A", "added_mass_coriolis_matrix")):
            ref = _cols(rigid_ref, key, 36)[k].reshape(6, 6)
            worst[key] = max(worst.get(key, 0.0), np.abs(vc.matrix(out[name]) - ref).max())
        ref = _cols(rigid_ref, "C_total", 36)[k].reshape(6, 6)
        total = vc.matrix(out["rigid_body_coriolis_matrix"]) + vc.matrix(out["added_mass_coriolis_matrix"])
        worst["C_total"] = max(worst.get("C_total", 0.0), np.abs(total - ref).max())

    for k, nu_r in enumerate(_cols(hydro_ref, "nu_r", 6)):
        out = vc.output_entries(vehicle, np.r_[np.zeros(6), nu_r], nothing)
        checks = {
            "M": (out["mass_matrix"], _cols(hydro_ref, "M_total", 36)[k]),
            "G": (out["restoring_matrix"], _cols(hydro_ref, "G", 36)[k]),
            "D": (out["hydrodynamic_loads.0.damping_matrix"], -_cols(hydro_ref, "D", 36)[k]),
            "damping_wrench": (out["hydrodynamic_loads.0.hydrodynamic_force"], _cols(hydro_ref, "tau_damp", 6)[k]),
            "cross_flow_wrench": (out["hydrodynamic_loads.1.hydrodynamic_force"],
                                  _cols(hydro_ref, "tau_crossflow", 6)[k]),
        }
        for name, (actual, reference) in checks.items():
            if name in ("M", "G", "D"):  # the references are row-major 6x6; the payload matrix is column-major
                reference = reference.reshape(6, 6).ravel(order="F")
            worst[name] = max(worst.get(name, 0.0), np.abs(np.asarray(actual).ravel() - reference).max())
    assert len(rigid_in) == len(_cols(hydro_ref, "nu_r", 6)) == 50
    assert max(worst.values()) <= G1_TOLERANCE, worst
