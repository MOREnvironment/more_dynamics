"""The ``Propulsor`` as a plugin of the shared repository's ``ForceProducer`` type, built with rpp's own
``ComponentContext`` inside a vehicle.

Written 2026-10-10, before the plugin exists (contract: ``propulsor_contract``).

One plugin, parts as selectors (the way the fin is built): the
parameters are one selector per socket (named as the socket, valued by a
variant of the contract's ``FORMS``), the numbers of every part variant, and
the skeleton's ``shaft_torque_scale``. Density is no parameter: the vehicle feeds
it, with the relative velocity, by name.

* Listed in ``plugins.json``, one class per file, a default for every
  parameter, the defaults the fixed jet (our servo and rotation), every
  number of every variant declared and valid, every variant selectable
  alone (with a compatible partner), each default documented in
  ``DEFAULTS.md`` under the plugin's own heading, no density parameter.
* In a vehicle: inputs ``drive_command`` then ``steering_command`` (a size-zero
  command is absent), states ``drive_state`` then ``servo_state`` before the vehicle's,
  first output ``generated_force`` (the wrench, BODY, about the CO), then
  ``thrust`` and ``steering_angle``.
* Equalities: the jet (his forms) on his frozen rows through the vehicle,
  output and dynamics; the fixed propeller equals the landed ``Propeller``
  plugin on the same states.
* Refusals when the vehicle is built: an unknown selector value names the
  socket and the choices; a drive / law pair in different units names both.

Author:    Enio Krizman
Date:      2026-10-10
"""

import importlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from propulsor_contract import (EXACT_TOLERANCE, FORMS, G2_TOLERANCE, SLOTS, DRIVE_UNITS, LAW_UNITS, expected_declared,
                                jet_forms_his, jet_parts, jet_reference, jet_values_his, max_diff, random_nu_r, rng,
                                remus_composition)

LIBRARY = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vehicle_models"))
import vehicle_contract as vc  # noqa: E402

RHO = vc.REMUS["water_density"]


def _propulsor():
    vc.plugins()
    return importlib.import_module("more_dynamics.plugins.force_producers.propulsor").Propulsor


def _flat(forms, part_values=None, scale=0.0):
    """Selectors and numbers as the plugin's flat parameter dictionary."""
    flat = {slot: forms[slot] for slot in SLOTS}
    for numbers in (part_values or {}).values():
        flat.update(numbers)
    flat["shaft_torque_scale"] = scale
    return flat


def _vehicle(actuator_params, extra=()):
    p = vc.plugins()
    return vc.remus(coriolis="munk_couplings_removed", current="none",
                    actuators=[p.PrescribedWrench, (_propulsor(), actuator_params), *extra])


def _evaluate(vessel, x, u):
    graph = vc.graph_of(vessel)
    return np.asarray(graph.step(x, u)).ravel(), vessel.signals(x, u)


def _default_numbers(cls):
    return {d.name: d.default_value for d in cls.PARAMETERS}


# --------------------------------------------------------------------------
# The plugin
# --------------------------------------------------------------------------
def test_the_propulsor_is_listed_has_one_class_in_its_file_and_a_default_for_every_parameter():
    cls = _propulsor()
    registry = json.loads((LIBRARY / "plugins.json").read_text())
    entries = {e["Name"]: e["Path"] for e in registry["Plugins"]}
    assert entries.get("Propulsor") == "more_dynamics/plugins/force_producers/propulsor.py"
    text = (LIBRARY / entries["Propulsor"]).read_text()
    assert text.count("\nclass ") == 1
    for d in cls.PARAMETERS:
        assert d.default_value is not None, d.name


def test_the_parameters_are_the_selectors_every_variants_numbers_and_the_shaft_torque_scale_and_no_density():
    cls = _propulsor()
    names = [d.name for d in cls.PARAMETERS]
    assert len(names) == len(set(names))
    for slot in SLOTS:
        assert slot in names, slot
    assert "shaft_torque_scale" in names and "water_density" not in names
    for (slot, variant) in FORMS:
        for number in expected_declared(slot, variant):
            assert number in names, (slot, variant, number)


def test_the_default_selectors_are_the_fixed_jet():
    defaults = _default_numbers(_propulsor())
    assert (defaults["drive"], defaults["steering"], defaults["mounting"], defaults["inflow"], defaults["law"],
            defaults["housing"]) == ("normalized_lag", "servo_lag_rate_angle", "steered_rotation", "none",
                                     "normalized_linear", "none")


def _partner(slot, variant):
    """A compatible drive/law pair for a variant being switched alone."""
    forms = dict(zip(SLOTS, ("normalized_lag", "servo_lag_rate_angle", "steered_rotation", "none", "normalized_linear",
                             "none")))
    forms[slot] = variant
    if slot == "drive":
        forms["law"] = {"1": "normalized_linear", "N": "thrust_power_limited"}.get(DRIVE_UNITS[variant],
                                                                                   "quadratic_asymmetric")
    if slot == "law":
        forms["drive"] = {"1": "normalized_lag", "N": "thrust_command"}.get(LAW_UNITS.get(variant, "rev/s"),
                                                                            "speed_command_rpm")
    return forms


def _numbers_of(forms):
    """Arbitrary in-range numbers for every part of ``forms`` (the plugin carries none: its numbers come from a named
    part): positive or free scalars 0.5, a fraction 0.5, a bound at zero -0.5, a vector [-1, 0.5, 0.25] (a point aft of
    the CO, a non-zero direction), the four open-water coefficients [0.4, 0.07, 0.18, 0.03]."""
    numbers = {"shaft_torque_scale": 0.0}
    for slot in SLOTS:
        for name, (shape, unit, minimum, exclusive, maximum) in expected_declared(slot, forms[slot]).items():
            if shape == (4, 1):
                numbers[name] = [0.4, 0.07, 0.18, 0.03]
            elif shape == (3, 1):
                numbers[name] = [-1.0, 0.5, 0.25]
            else:
                numbers[name] = -0.5 if maximum == 0.0 else 0.5
    return numbers


@pytest.mark.parametrize("slot, variant", sorted(FORMS))
def test_every_variant_can_be_selected_with_numbers_given_as_a_named_part_gives_them(slot, variant):
    forms = _partner(slot, variant)
    vessel = _vehicle({**{s: forms[s] for s in SLOTS}, **_numbers_of(forms)})
    assert vc.input_names(vessel)[0] == "desired_wrench"
    assert vessel.graph().outputDescription        # the vehicle assembles and exposes its outputs


def test_the_defaults_are_documented_under_the_plugins_own_heading_in_defaults_md():
    text = (LIBRARY / "DEFAULTS.md").read_text()
    match = re.search(r"^## `more_dynamics::Propulsor`.*?(?=^## |\Z)", text, re.S | re.M)
    assert match, "DEFAULTS.md has no section for more_dynamics::Propulsor"
    section = match.group(0)
    for d in _propulsor().PARAMETERS:
        assert f"`{d.name}`" in section, d.name


# --------------------------------------------------------------------------
# In a vehicle
# --------------------------------------------------------------------------
def test_the_jet_in_a_vehicle_orders_inputs_states_and_outputs():
    forms = jet_forms_his()
    vessel = _vehicle(_flat(forms, jet_values_his(jet_parts()["trimmed"])))
    graph = vessel.graph()
    assert vc.input_names(vessel) == ["desired_wrench", "drive_command", "steering_command"]
    assert [(d.name, d.size) for d in graph.stateDescription] == [("drive_state", 1), ("servo_state", 1), ("state", 12)]


def test_a_fixed_thruster_has_no_steering_command_and_no_states():
    forms, values, scale, _ = remus_composition("linearized")
    vessel = _vehicle(_flat(forms, values, scale))
    assert vc.input_names(vessel) == ["desired_wrench", "drive_command"]
    assert [(d.name, d.size) for d in vessel.graph().stateDescription] == [("state", 12)]


@pytest.mark.parametrize("name", sorted(jet_parts()))
def test_his_jet_through_the_plugin_equals_his_frozen_rows_output_and_dynamics(name):
    ref = jet_reference()[name]
    vessel = _vehicle(_flat(jet_forms_his(), jet_values_his(jet_parts()[name])))
    worst_out = worst_dot = 0.0
    for state, command, output, dot in list(zip(ref.state, ref.command, ref.output, ref.dot))[::3]:
        x = np.concatenate([state, np.zeros(12)])               # the vehicle at rest: relative velocity zero
        xdot, out = _evaluate(vessel, x, np.concatenate([np.zeros(6), command]))
        worst_out = max(worst_out, max_diff(out["actuators.1.generated_force"], output))
        worst_dot = max(worst_dot, max_diff(xdot[:2], dot))
    assert worst_out <= EXACT_TOLERANCE and worst_dot <= EXACT_TOLERANCE, (worst_out, worst_dot)


def test_the_fixed_propeller_through_the_plugin_equals_the_landed_propeller_plugin():
    p = vc.plugins()
    forms, values, scale, _ = remus_composition("linearized", open_water_coefficients="given")
    propulsor_vessel = _vehicle(_flat(forms, values, scale))
    propeller_vessel = vc.remus(coriolis="munk_couplings_removed", current="none",
                                actuators=[p.PrescribedWrench, (p.Propeller, vc.pick(p.Propeller, vc.remus_values()))])
    g = rng()
    worst = 0.0
    for nu, rpm in zip(random_nu_r(g, 40, speed=2.5, rate=0.4), g.uniform(-1800, 1800, 40)):
        x = np.concatenate([np.zeros(6), nu[:6]])
        u = np.concatenate([np.zeros(6), [rpm]])
        a = _evaluate(propulsor_vessel, x, u)[1]["actuators.1.generated_force"]
        b = _evaluate(propeller_vessel, x, u)[1]["actuators.1.generated_force"]
        worst = max(worst, max_diff(a, b))
    assert worst <= G2_TOLERANCE, worst


# --------------------------------------------------------------------------
# Refusals
# --------------------------------------------------------------------------
@pytest.mark.parametrize("slot, bad", [("drive", "first_order"), ("steering", "rudder"), ("mounting", "swivel"),
                                       ("inflow", "wake"), ("law", "four_quadrant"), ("housing", "leg")])
def test_an_unknown_selector_value_is_refused_naming_the_socket_and_the_choices(slot, bad):
    with pytest.raises(ValueError, match=f"Propulsor: {slot} must be one of"):
        _vehicle({slot: bad})


def test_a_drive_and_a_law_in_different_units_are_refused_when_the_vehicle_is_built_naming_both():
    with pytest.raises(ValueError) as error:
        _vehicle({"drive": "speed_command_rpm", "law": "normalized_linear"})
    text = str(error.value)
    assert "drive" in text and "law" in text and "rev/s" in text and "'1'" in text, text
