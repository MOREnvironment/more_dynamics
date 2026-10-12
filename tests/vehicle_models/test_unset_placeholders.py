"""Neutral placeholders: the plugins carry no vehicle, device or site number. Every value of a plugin is declared
unset (``unset()``, an empty list) and the composition gives it; a plugin built with nothing set refuses in
``initialize`` with one message that names the plugin and lists every value its selected forms read, with unit and
meaning; and exactly the listed values are what the build needs.

Per plugin (13 of ours): (1) every declared default is unset, a text selector, a switch or one of the three named
numerical settings; (2) built with no value set, it refuses and the message lists the values; (3) given exactly the
listed values (from the reference part or data file) it builds; (4) leaving any one of them out refuses again,
naming that one. A vehicle type's children keep their reference values so the refusal is the vehicle's own.

Author:    Enio Krizman
Date:      2026-10-12
"""
import re

import numpy as np
import pytest

import vehicle_contract as vc
from more_transformations.more_casadi_transformations import is_unset

NUMERICAL_SETTINGS = {"integration_max_step", "smooth_speed_epsilon", "ittc_reynolds_floor"}
LEAVES = ["SubmergedRestoring", "SurfaceRestoring", "AuvHullLoads", "SurfaceHullLoads", "CoefficientHullLoads",
          "Site", "Fin", "FinPairsDeflectionOnly", "Propeller"]
VEHICLES = {"SpheroidAuv": ("SubmergedRestoring", "AuvHullLoads"), "Catamaran": ("SurfaceRestoring", "SurfaceHullLoads"),
            "Monohull": ("SurfaceRestoring", "SurfaceHullLoads"),
            "DerivativeAuv": ("SubmergedRestoring", "CoefficientHullLoads")}
PLUGIN_NAMES = LEAVES + list(VEHICLES)


def _build_leaf(cls, values):
    context = vc.context_of(cls, values, reference=False)
    context.initialize()
    return context


def _build_vehicle(cls, values):
    p = vc.plugins()
    hydrostatics, hydrodynamics = (getattr(p, name) for name in VEHICLES[cls.__name__])
    if cls is p.Monohull:  # one hull: the surface restoring of the Otter part has two
        hydrostatics = (hydrostatics, {"hull_count": 1})
    return vc.build_context(cls, values, hydrostatics=hydrostatics, hydrodynamics=hydrodynamics,
                            actuators=[p.PrescribedWrench], reference=False)


def _build(cls, values):
    return (_build_vehicle if cls.__name__ in VEHICLES else _build_leaf)(cls, values)


def _refusal(cls, values):
    """The listed parameter names of the refusal, and its message."""
    UnsetValues = vc.plugins().UnsetValues
    with pytest.raises(UnsetValues) as caught:
        _build(cls, values)
    message = str(caught.value)
    assert message.startswith(f"{cls.__name__}: not set in the composition -- missing parameter(s): "), message
    return re.findall(r"'(\w+)' \[", message), message


@pytest.mark.parametrize("name", PLUGIN_NAMES)
def test_every_default_is_unset_a_selector_a_switch_or_a_named_numerical_setting(name):
    cls = getattr(vc.plugins(), name)
    for d in cls.PARAMETERS:
        value = d.default_value
        assert (is_unset(value) or isinstance(value, (str, bool)) or d.name in NUMERICAL_SETTINGS), \
            f"{name}.{d.name} defaults to {value!r}: a vehicle number belongs in a named part"


@pytest.mark.parametrize("name", PLUGIN_NAMES)
def test_built_with_nothing_set_it_refuses_naming_every_value_with_unit_and_meaning(name):
    cls = getattr(vc.plugins(), name)
    listed, message = _refusal(cls, {})
    assert listed, message
    assert len(set(listed)) == len(listed), f"a value listed twice: {message}"
    declared = {d.name for d in cls.PARAMETERS if is_unset(d.default_value)}
    assert set(listed) <= declared, set(listed) - declared
    assert re.search(r"'\w+' \[[^\]]*\] \([^)]+\)", message), "unit and meaning are named"
    assert "set them in the part's params/parameters.py" in message


@pytest.mark.parametrize("name", PLUGIN_NAMES)
def test_the_listed_values_are_exactly_what_the_build_needs(name):
    p = vc.plugins()
    cls = getattr(p, name)
    listed, _ = _refusal(cls, {})
    reference = vc.reference_values(cls)
    assert set(listed) <= set(reference), f"no reference value for {sorted(set(listed) - set(reference))}"
    given = {n: reference[n] for n in listed}
    _build(cls, given)  # one listing was enough: nothing else is missing
    for n in listed:  # and each of them is needed
        others = {k: v for k, v in given.items() if k != n}
        relisted, message = _refusal(cls, others)
        assert n in relisted, (n, message)


def test_a_vehicle_whose_child_is_unset_refuses_in_the_child_before_the_vehicle_is_built():
    """rpp initialises children first: the unset child's own refusal names the child, and the vehicle's
    ``graph()`` is never asked for."""
    p = vc.plugins()
    with pytest.raises(p.UnsetValues, match="^SurfaceHullLoads: not set in the composition"):
        vc.build_context(p.Catamaran, hydrostatics=p.SurfaceRestoring,
                         hydrodynamics=vc.context_of(p.SurfaceHullLoads, reference=False))


def test_zero_one_and_false_are_values_and_only_the_empty_list_is_unset():
    assert not is_unset(0) and not is_unset(0.0) and not is_unset([0.0]) and not is_unset(False)
    assert not is_unset("") and not is_unset(None) and is_unset([])
