"""Part plugins as plugin classes in Luka's form: every ``PARAMETERS`` entry has
a default, the defaults are the frozen reference values the trees carry, and
the model layer imports no plugin.

Author:    Enio Krizman
Date:      2026-10-08
"""

import importlib
import json
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc
from vehicle_contract import make_trees

LIBRARY = Path(__file__).resolve().parents[2]
OURS = {"SiteAtLatitude", "SiteGiven", "NoCurrent", "HorizontalCurrentYawRateTerms",
        "HorizontalCurrentFullRotationRate", "ProlateSpheroidMainDimensions", "TwinPontoons",
        "HomogeneousSpheroid", "HullWithPointPayload", "LambSpheroid", "ScaledDerivatives", "KirchhoffFull",
        "MunkCouplingsRemoved", "SubmergedNeutral", "TwinHullMetacentricEquilibriumDraft",
        "TimeConstantDampingSubmerged", "TimeConstantDampingSurface", "HullLiftDrag", "CrossFlowStrip",
        "SurgeResistanceIttc", "CircularCylinderReynolds", "RectangularSectionHoerner", "PrescribedWrench",
        "FinPairsDeflectionOnly", "PropellerLinearizedOpenWater", "MarineCraft6DOF"}


def _classes():
    vc.builder()
    registry = json.loads((LIBRARY / "plugins.json").read_text())
    found = {}
    for entry in registry["Plugins"]:
        if entry["Name"] in OURS:
            module = importlib.import_module(entry["Path"][:-3].replace("/", "."))
            found[entry["Name"]] = getattr(module, entry["Name"])
    return found


def test_every_part_of_ours_is_listed_in_plugins_json_and_has_a_class():
    assert set(_classes()) == OURS


def test_every_parameter_has_a_default():
    for name, cls in _classes().items():
        for description in cls.PARAMETERS:
            assert description.default_value is not None, (name, description.name)


def _nodes(tree):
    yield tree
    for item in tree["children"].values():
        for child in item if isinstance(item, list) else [item]:
            yield from _nodes(child)


@pytest.mark.parametrize("tree_name", ["remus100_full", "remus100_hull", "otter", "otter_current"])
def test_defaults_equal_the_frozen_values_of_the_tree(tree_name):
    """The literal defaults in each part are the numbers of the frozen parameter files (REMUS 100 one-value set,
    Otter inputs); the vehicle-level settings and the current of a tree that sets one are the tree's own."""
    classes = _classes()
    trees = make_trees.gate_trees()
    skip = {"open_parameters", "open_inputs", "diagnostic_outputs"}
    for node in _nodes(trees[tree_name]):
        if node["plugin"] in ("HorizontalCurrentFullRotationRate",) and tree_name == "otter_current":
            continue  # the rigid-body test's current (a fixture), not a default
        defaults = {d.name: d.default_value for d in classes[node["plugin"]].PARAMETERS}
        for name, value in node["params"].items():
            if name in skip:
                continue
            assert name in defaults, (node["plugin"], name)
            assert np.array_equal(np.asarray(defaults[name], dtype=float), np.asarray(value, dtype=float)), \
                (node["plugin"], name, defaults[name], value)


def test_the_model_layer_imports_no_plugin():
    models = LIBRARY / "more_dynamics" / "models"
    offenders = [str(path.relative_to(LIBRARY)) for path in models.rglob("*.py")
                 if "more_dynamics.plugins" in path.read_text() or "from ..plugins" in path.read_text()]
    assert not offenders


@pytest.mark.parametrize("vehicle, commands", [
    ("remus100", ["fin_deflection_command", "shaft_speed_command"]),
    ("otter", []),
])
def test_user_vehicles_build_from_their_script_descriptions(vehicle, commands):
    """The two vehicles of ``scripts/vehicles/.rppws`` build in rpp's builder, close to nothing open, and are the
    gate trees' parts with the physics Coriolis form."""
    script = LIBRARY / "scripts" / "vehicles" / ".rppws" / "script_descriptions" / f"{vehicle}.json"
    vessel = vc.build(vehicle, script)
    assert vc.input_names(vessel) == commands
    assert [d.name for d in vessel.graph().outputDescription] == ["output"]
