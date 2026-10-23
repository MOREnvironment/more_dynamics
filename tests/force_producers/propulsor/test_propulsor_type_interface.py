"""Type-interface checks of the propulsor skeleton and its parts.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

* Every part form has the coupling inputs and outputs of its type, names,
  sizes and order, and then its own declared parameters (names, shapes,
  units, admissible ranges as the contract lists them); nothing else, and no
  part declares ``water_density`` as a parameter. A drive module names the
  unit of its output, a law module the unit it accepts.
* ``propulsor_casadi`` refuses, naming the slot: a missing or unknown slot, a
  part without a coupling name of its type (input or output), a coupling of
  the wrong size, a part outside ``law`` and ``housing`` with an input
  ``water_density``, a drive whose output depends on the load torque, and a
  drive / law pair whose units differ (naming both slots and both units; the
  unit check is a proposal not yet accepted). The rpp context
  builder instantiates a child by its plugin name and does not compare its
  type with the slot's, so the skeleton is where a wrong part is caught.
* The skeleton's signature: the four fed inputs, its one declared number,
  then every open part parameter as ``"<slot>.<name>"``; command and state
  are ``[drive; steering]``.
* ``check_propulsor_values`` refuses a negative shaft-torque scale.

Author:    Enio Krizman
Date:      2026-10-10
"""

import casadi as ca
import numpy as np
import pytest

from propulsor_contract import (
    COUPLINGS,
    DECLARED,
    DENSITY_SLOTS,
    DRIVE_UNITS,
    FORMS,
    LAW_UNITS,
    PRODUCER_INPUTS,
    PRODUCER_OUTPUTS,
    SKELETON_DECLARED,
    SLOTS,
    expected_declared,
    form_module,
    max_diff,
    part,
    part_block,
    part_declared,
    propulsor,
    propulsor_from_forms,
    run,
    skeleton,
    units_of,
)


def _fake(name, inputs, outputs):
    """A stand-in part: inputs {name: rows}, outputs {name: rows}, each output an expression of the first input."""
    syms = {n: ca.SX.sym(n, rows) for n, rows in inputs.items()}
    first = syms[next(iter(syms))][0] if syms and next(iter(syms.values())).numel() else ca.SX(1.0)
    outs = {n: ca.repmat(first, rows, 1) if rows else ca.SX(0, 1) for n, rows in outputs.items()}
    return ca.Function(name, list(syms.values()), list(outs.values()), list(syms), list(outs))


JET_FORMS = {"drive": "speed_command_rpm", "steering": "servo_lag_rate_angle", "mounting": "steered_rotation",
             "inflow": "surge_component", "law": "jet_momentum_flux_given", "housing": "none"}
JET_VALUES = {"drive": {"min_command": -6000.0, "max_command": 6000.0},
              "steering": {"max_deflection": 0.5, "max_rate": 1.0, "time_constant": 0.2},
              "mounting": {"position": [-1.0, 0.0, 0.2], "trim_angle": 0.05, "yaw_bias": 0.0},
              "inflow": {"thrust_deduction": 0.05, "advance_speed_factor": 1.0, "wake_fraction": 0.1},
              "law": {"nozzle_area": 7e-3, "flow_coefficient": 2e-3}}


def _parts(forms, values=None):
    values = values or {}
    return {slot: part(slot, forms[slot], values.get(slot) if part_declared(slot, forms[slot]) else None)
            for slot in SLOTS}


def _good_parts():
    return _parts(JET_FORMS, JET_VALUES)


def _units():
    return units_of(JET_FORMS)


# --------------------------------------------------------------------------
# Parts: couplings of their type, declared parameters
# --------------------------------------------------------------------------
@pytest.mark.parametrize("slot, variant", sorted(FORMS))
def test_part_has_the_couplings_of_its_type_then_its_parameters(slot, variant):
    f = part_block(slot, variant)
    inputs, outputs = COUPLINGS[slot]
    declared = [d.name for d in part_declared(slot, variant)]
    assert declared == list(expected_declared(slot, variant)), declared
    assert f.name_in() == [*inputs, *declared], f.name_in()
    assert f.name_out() == list(outputs), f.name_out()
    for name, rows in inputs.items():
        if rows is not None:
            assert f.size_in(name) == (rows, 1), (name, f.size_in(name))
    for name, rows in outputs.items():
        if rows is not None:
            assert f.size_out(name) == (rows, 1), (name, f.size_out(name))


@pytest.mark.parametrize("slot, variant", sorted(FORMS))
def test_part_declares_its_numbers_as_the_contract_lists_them(slot, variant):
    expected = expected_declared(slot, variant)
    for d in part_declared(slot, variant):
        shape, unit, minimum, exclusive, maximum = expected[d.name]
        assert (d.shape, d.unit, d.minimum, d.minimum_exclusive, d.maximum) == (shape, unit, minimum, exclusive, maximum), d
        assert d.meaning, d.name


@pytest.mark.parametrize("slot, variant", sorted(FORMS))
def test_part_is_one_casadi_function_and_declares_no_density(slot, variant):
    f = part_block(slot, variant)
    assert isinstance(f, ca.Function)
    declared = [d.name for d in part_declared(slot, variant)]
    assert "water_density" not in declared
    if slot not in DENSITY_SLOTS:
        assert "water_density" not in f.name_in()


def test_the_contract_table_and_the_declaration_table_list_the_same_forms():
    assert sorted(DECLARED) == sorted(FORMS)


@pytest.mark.parametrize("variant", sorted(DRIVE_UNITS))
def test_a_drive_module_names_the_unit_of_its_output(variant):
    assert form_module("drive", variant).DRIVE_OUTPUT_UNIT == DRIVE_UNITS[variant]


@pytest.mark.parametrize("slot, variant", sorted(k for k in FORMS if k[0] == "law"))
def test_a_law_module_names_the_unit_it_accepts(slot, variant):
    assert form_module(slot, variant).LAW_INPUT_UNIT == LAW_UNITS.get(variant, "rev/s")


# --------------------------------------------------------------------------
# Skeleton
# --------------------------------------------------------------------------
def test_skeleton_declares_its_one_number():
    declared = skeleton().propulsor_parameters()
    assert [d.name for d in declared] == list(SKELETON_DECLARED)
    for d in declared:
        shape, unit, minimum, exclusive, maximum = SKELETON_DECLARED[d.name]
        assert (d.shape, d.unit, d.minimum, d.minimum_exclusive, d.maximum) == (shape, unit, minimum, exclusive, maximum), d
        assert d.meaning


def test_skeleton_signature():
    f = skeleton().propulsor_casadi(_good_parts(), _units())
    names = [d.name for d in skeleton().propulsor_parameters()]
    assert f.name_in()[:len(PRODUCER_INPUTS) + len(names)] == [*PRODUCER_INPUTS, *names], f.name_in()
    assert f.name_out() == PRODUCER_OUTPUTS
    assert f.size_in("command") == (2, 1)      # drive command, then steering command
    assert f.size_in("state") == (1, 1)        # no drive state, one servo state
    assert f.size_in("nu_r") == (6, 1) and f.size_in("water_density") == (1, 1)
    assert f.size_out("tau") == (6, 1) and f.size_out("state_dot") == (1, 1)
    assert f.size_out("thrust") == (1, 1) and f.size_out("steering_angle") == (1, 1)


def test_open_part_parameters_become_inputs_named_slot_dot_name():
    parts = _good_parts()
    parts["law"] = part("law", "jet_momentum_flux_given")  # parameters left open
    f = skeleton().propulsor_casadi(parts, _units())
    assert f.name_in()[-2:] == ["law.nozzle_area", "law.flow_coefficient"], f.name_in()


def test_command_and_state_are_drive_then_steering():
    forms = {**JET_FORMS, "drive": "normalized_lag", "law": "normalized_linear", "mounting": "steered_rotation"}
    values = {**JET_VALUES, "drive": {"thrust_rise_time": 0.5, "thrust_fall_time": 0.2},
              "law": {"max_thrust": 100.0, "thrust_coefficient": 1.0}}
    f = propulsor_from_forms(forms, values)
    assert f.size_in("command") == (2, 1) and f.size_in("state") == (2, 1)
    # drive state 0.5, servo state 0.1; commands: effort 1, angle 0.3
    out = run(f, [1.0, 0.3], np.zeros(6), 1025.0, state=[0.5, 0.1])
    assert out["state_dot"].shape == (2,)
    assert abs(out["state_dot"][0] - (1.0 - 0.5) / 0.5) <= 1e-12            # the drive's lag (rise time 0.5 s)
    assert abs(out["state_dot"][1] - (0.3 - 0.1) / 0.2) <= 1e-12            # the servo's lag (0.2 s) inside its rate limit


@pytest.mark.parametrize("change, words", [
    (lambda p: p.pop("law"), ["law"]),
    (lambda p: p.update(propeller_hub=p["housing"]), ["propeller_hub"]),
    (lambda p: p.update(inflow=_fake("bad_inflow", {"nu_r": 6, "thrust_axis": 3},
                                     {"advance_speed": 1, "thrust_factor": 1, "torque_factor": 1})),
     ["inflow", "thrust_position"]),
    (lambda p: p.update(inflow=_fake("bad_inflow", {"nu_r": 6, "thrust_axis": 3, "thrust_position": 3},
                                     {"advance_speed": 1, "thrust_factor": 1})), ["inflow", "torque_factor"]),
    (lambda p: p.update(mounting=_fake("bad_mount", {"steering_angle": 2}, {"thrust_axis": 3, "thrust_position": 3})),
     ["mounting", "steering_angle"]),
    (lambda p: p.update(mounting=_fake("bad_mount", {"steering_angle": 1}, {"thrust_axis": 2, "thrust_position": 3})),
     ["mounting", "thrust_axis"]),
    (lambda p: p.update(law=_fake("bad_law", {"drive_output": 1, "advance_speed": 1}, {"thrust": 1, "torque": 1})),
     ["law", "water_density"]),
    (lambda p: p.update(drive=_fake("bad_drive", {"command": 1, "drive_state": 0, "shaft_torque": 1},
                                    {"shaft_speed": 1, "drive_state_dot": 0, "thrust_scale": 1, "thrust_limit": 1})),
     ["drive", "drive_output"]),
    (lambda p: p.update(drive=_fake("dense_drive", {"command": 1, "drive_state": 0, "shaft_torque": 1, "water_density": 1},
                                    {"drive_output": 1, "drive_state_dot": 0, "thrust_scale": 1, "thrust_limit": 1})),
     ["drive", "water_density"]),
    (lambda p: p.update(steering=_fake("bad_servo", {"command": 1, "servo_state": 1},
                                       {"deflection": 1, "servo_state_dot": 0})), ["steering", "state"]),
    (lambda p: p.update(housing=_fake("bad_housing", {"nu_r": 6, "thrust_axis": 3, "thrust_position": 3,
                                                       "steering_angle": 1, "water_density": 1}, {"tau": 3})),
     ["housing", "tau"]),
])
def test_skeleton_refuses_a_part_that_does_not_fit_its_slot(change, words):
    parts = _good_parts()
    change(parts)
    with pytest.raises(ValueError) as error:
        skeleton().propulsor_casadi(parts, _units())
    for word in words:
        assert word in str(error.value), (word, str(error.value))


def test_a_drive_whose_output_depends_on_the_load_is_refused():
    """Load-dependent shaft speed must be a state: a drive that returns the speed as a function of the load torque
    would make the speed, torque and law a loop."""
    command, state, load = ca.SX.sym("command", 1), ca.SX.sym("drive_state", 0), ca.SX.sym("shaft_torque")
    loaded = ca.Function("loaded_drive", [command, state, load],
                         [command - 0.01 * load, ca.SX(0, 1), ca.SX(1.0), ca.SX(ca.inf)],
                         ["command", "drive_state", "shaft_torque"],
                         ["drive_output", "drive_state_dot", "thrust_scale", "thrust_limit"])
    parts = _good_parts()
    parts["drive"] = loaded
    with pytest.raises(ValueError) as error:
        skeleton().propulsor_casadi(parts, _units())
    assert "drive" in str(error.value) and "load" in str(error.value), str(error.value)


def test_a_drive_and_a_law_in_different_units_are_refused_naming_both():  # the unit check is a proposal, not yet accepted
    forms = {**JET_FORMS, "drive": "normalized_lag"}                   # effort, unit "1"
    parts = _parts(forms)
    with pytest.raises(ValueError) as error:
        skeleton().propulsor_casadi(parts, units_of(forms))             # law jet_momentum_flux accepts "rev/s"
    text = str(error.value)
    for word in ("drive", "law", "'1'", "'rev/s'"):
        assert word in text, (word, text)


def test_the_matching_pairs_are_accepted():
    for drive, law in (("speed_command_rpm", "jet_momentum_flux_given"), ("speed_command_rad_s", "quadratic_asymmetric"),
                       ("normalized_lag", "normalized_linear")):
        units = units_of({"drive": drive, "law": law})
        assert units["drive"] == units["law"], (drive, law, units)


@pytest.mark.parametrize("value", [-0.1, -1.0])
def test_a_negative_shaft_torque_scale_is_refused(value):
    with pytest.raises(ValueError, match="shaft_torque_scale"):
        skeleton().check_propulsor_values({"shaft_torque_scale": value})


def test_the_fed_density_scales_a_law_that_reads_it():
    """One density enters, from the vehicle: a law that reads it scales with it (momentum flux: thrust ~ rho)."""
    f = propulsor_from_forms(JET_FORMS, JET_VALUES)
    base = run(f, [3000.0, 0.0], [0, 0, 0, 0, 0, 0], 1000.0, state=[0.0])
    more = run(f, [3000.0, 0.0], [0, 0, 0, 0, 0, 0], 1030.0, state=[0.0])
    assert max_diff(more["tau"], 1.03 * base["tau"]) <= 1e-9 * max(1.0, float(np.max(np.abs(more["tau"]))))
