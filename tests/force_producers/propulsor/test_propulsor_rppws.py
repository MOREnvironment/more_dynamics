"""The ``Propulsor`` built through rpp: its named parts of the library's ``.rppws`` and vehicles whose actuators are
those parts, built by rpp's own ``ComponentContextBuilder`` (not assembled by hand as the other propulsor tests do).

* Every named part of ``more_dynamics::Propulsor`` loads with rpp's parameter store, names only parameters the plugin
  declares, and builds the plugin.
* The committed REMUS 100 (MSS) vehicle links the named part ``REMUS 100 propulsor`` as its propulsor; its equality
  with the earlier composition (the propeller plugin and the deflection-only fin pairs) and its frozen runs are in
  ``tests/vehicle_models/test_named_vehicles.py`` and ``test_rppws_runs.py``.
* An Otter whose two thrusters are the named parts ``Otter port thruster`` and ``Otter starboard thruster``, linked
  into its ``actuators`` list, is built by rpp's builder and gives the same dynamics as the Otter driven by the
  differential-thruster block's wrench. The two thrusters' inputs are both called ``drive_command``: the standard names
  an input by the plugin, so a command schedule cannot tell them apart by name (positional use is unaffected); this is
  recorded in the test as the platform's behaviour, not worked around.

Author:    Enio Krizman
Date:      2026-10-23
"""

import importlib
import json
import shutil
import sys
import uuid
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vehicle_models"))
import vehicle_contract as vc  # noqa: E402

from propulsor_contract import rng, two_thruster_values  # noqa: E402

LIBRARY = Path(__file__).resolve().parents[3]
RPPWS = LIBRARY / ".rppws"
PARTS = ("REMUS 100 propulsor", "Otter port thruster", "Otter starboard thruster")
OTTER_PART = "Otter (MSS)"
LIBRARY_NAME = "more_dynamics"


def _propulsor_class():
    vc.plugins()
    return importlib.import_module("more_dynamics.plugins.force_producers.propulsor").Propulsor


def _part_dir(name):
    found = [d.parent for d in RPPWS.joinpath("parts").rglob("description.json")
             if json.loads(d.read_text()).get("Name") == name]
    assert len(found) == 1, (name, found)
    return found[0]


def _builder(workspace):
    try:
        from rpp_py.context_builder import ComponentContextBuilder
        from rpp_py.data_manager import DataManager
        from rpp_plugin_registrator.library_manager import LibraryManager
    except Exception as exc:  # pragma: no cover - environment without rpp
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
    vc.plugins()
    data_manager = DataManager(library_manager=LibraryManager(), workspace_path=str(workspace))
    return ComponentContextBuilder(data_manager=data_manager)


def _build(workspace, description_path, configuration):
    context = _builder(workspace).build_script_from_description_path(str(description_path), configuration=configuration)
    context.initialize()
    return context.get_component("vessels")[0]


def _workspace(tmp_path, parts):
    """A ``.rppws`` with the library's script and the part folders ``parts`` (folder names under ``parts/``), each
    copied whole."""
    root = tmp_path / LIBRARY_NAME  # the script runner finds a script description of a workspace named as the library
    (root / "scripts" / "vehicles").mkdir(parents=True)
    shutil.copy(LIBRARY / "scripts" / "vehicles" / "vehicle_simulation.py",
                root / "scripts" / "vehicles" / "vehicle_simulation.py")
    for folder in parts:
        target = root / ".rppws" / "parts" / folder.parent.name / folder.name
        shutil.copytree(folder, target)
    return root


def _description(tmp_path, configuration, vehicle_id, plugin_name, commands=None):
    configuration_entry = {"Description": "test", "Components": {"vessels": [{"Id": vehicle_id, "PluginName": plugin_name}]}}
    if commands:
        configuration_entry["Commands"] = commands
    description = {"ScriptPath": "scripts/vehicles/vehicle_simulation.py", "Language": "python",
                   "Configurations": {configuration: configuration_entry}, "ActiveConfiguration": configuration,
                   "Spec": {"vessels": "List[more_dynamics::VehicleModel3D]"}}
    path = tmp_path / ".rppws" / "script_descriptions" / "vehicle_simulation.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(description, indent=4))
    return path


# --------------------------------------------------------------------------
# The named parts
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", PARTS)
def test_a_named_part_names_only_declared_parameters_and_builds_the_plugin(name):
    cls = _propulsor_class()
    values = vc.reference_part(name)
    declared = {d.name for d in cls.PARAMETERS}
    assert set(values) <= declared, sorted(set(values) - declared)
    context = vc.context_of(cls, values, reference=False)
    context.initialize()
    assert context.get_instance().graph().outputDescription[0].name == "generated_force"


@pytest.mark.parametrize("name", PARTS)
def test_a_named_part_describes_a_propulsor_force_producer_with_no_vehicle(name):
    description = json.loads((_part_dir(name) / "description.json").read_text())
    assert description["PluginName"] == "more_dynamics::Propulsor"
    assert description["PluginType"] == "more_dynamics::ForceProducer"
    assert description["ParentComponentInfo"] is None
    assert (_part_dir(name) / "SOURCE.md").is_file()


# --------------------------------------------------------------------------
# The Otter with two named thrusters linked into its actuator list
# --------------------------------------------------------------------------
def _otter_with_thrusters(tmp_path):
    vehicle_dir = _part_dir(OTTER_PART)
    thrusters = [_part_dir("Otter port thruster"), _part_dir("Otter starboard thruster")]
    workspace = _workspace(tmp_path, [vehicle_dir, *thrusters])
    copied = workspace / ".rppws" / "parts" / vehicle_dir.parent.name / vehicle_dir.name
    parent = json.loads((copied / "description.json").read_text())
    entries = []
    for thruster in thrusters:
        thruster_id = json.loads((thruster / "description.json").read_text())["Id"]
        link_id = str(uuid.uuid4())
        folder = copied / "subcomponents" / link_id
        folder.mkdir(parents=True)
        (folder / "description.json").write_text(json.dumps({
            "Id": link_id, "Name": json.loads((thruster / "description.json").read_text())["Name"],
            "LinkedComponentId": thruster_id, "LinkedComponentWorkspace": LIBRARY_NAME,
            "ParentComponentInfo": {"Id": parent["Id"], "PluginType": "more_dynamics::VehicleModel3D",
                                    "PluginName": "more_dynamics::Catamaran", "SlotName": "actuators",
                                    "Library": LIBRARY_NAME, "IsLinked": False}}, indent=4))
        entries.append({"Id": link_id, "PluginType": "more_dynamics::ForceProducer", "PluginName": "more_dynamics::Propulsor",
                        "SlotName": "actuators", "Library": LIBRARY_NAME, "IsLinked": True})
    parent["Subcomponents"]["actuators"] = entries
    (copied / "description.json").write_text(json.dumps(parent, indent=4))
    return workspace, parent["Id"]


def test_the_otter_with_two_linked_named_thrusters_equals_the_otter_driven_by_the_thruster_blocks_wrench(tmp_path):
    from more_dynamics.models.thruster.differential_thruster import (check_differential_thruster_values,
                                                                      differential_thruster_casadi,
                                                                      differential_thruster_parameters)
    from more_transformations.more_casadi_transformations import freeze
    workspace, vehicle_id = _otter_with_thrusters(tmp_path)
    description = _description(workspace, "otter_thrusters", vehicle_id, "more_dynamics::Catamaran")
    built = _build(workspace, description, "otter_thrusters")
    assert type(built).__name__ == "Catamaran"
    names = vc.input_names(built)
    assert names == ["drive_command", "drive_command"]            # one name for both: the standard names an input by its plugin
    block_values = {**two_thruster_values()[0]}
    port, starboard = (vc.reference_part(n) for n in ("Otter port thruster", "Otter starboard thruster"))
    block_values.update({"positive_thrust_coefficients": [port["positive_thrust_coefficient"], starboard["positive_thrust_coefficient"]],
                         "negative_thrust_coefficients": [port["negative_thrust_coefficient"], starboard["negative_thrust_coefficient"]],
                         "thruster_positions": [port["position"], starboard["position"]],
                         "thruster_directions": [port["direction"], starboard["direction"]]})
    # the limits of the part follow the thrust limits through n = sqrt(F / k): give the block the same ones
    block_values["max_forward_thrust"] = [port["positive_thrust_coefficient"] * port["max_command"] ** 2,
                                          starboard["positive_thrust_coefficient"] * starboard["max_command"] ** 2]
    block_values["max_reverse_thrust"] = [port["negative_thrust_coefficient"] * port["min_command"] ** 2,
                                          starboard["negative_thrust_coefficient"] * starboard["min_command"] ** 2]
    declared = differential_thruster_parameters()
    thrusters = freeze(differential_thruster_casadi(), declared, check_differential_thruster_values(block_values))
    reference = vc.otter(actuators=[vc.plugins().PrescribedWrench])
    a, b = vc.graph_of(built), vc.graph_of(reference)
    g = rng()
    worst = 0.0
    for _ in range(40):
        x = np.concatenate([g.uniform(-0.3, 0.3, 6), g.uniform(-0.5, 2.0, 3), g.uniform(-0.2, 0.2, 3)])
        n = g.uniform(-1.2, 1.2, 2) * np.array([port["max_command"], starboard["max_command"]])
        tau = np.asarray(thrusters(n=n, nu_r=np.zeros(6))["tau"]).ravel()   # the thrusters take nothing from the water
        xa = np.asarray(a.step(x, n)).ravel()
        xb = np.asarray(b.step(x, tau)).ravel()
        worst = max(worst, float(np.max(np.abs(xa - xb) / np.maximum(1.0, np.abs(xb)))))
    assert worst <= 1e-10, worst


# --------------------------------------------------------------------------
# The editor's round trip
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", PARTS)
def test_a_named_part_survives_the_editors_save_and_gives_the_same_propulsor(name, tmp_path):
    """Saving the part's parameters the way rpp's editor does (it rewrites the file without the comments) and reading
    them back gives the same values, and the plugin built from them has the same wrench and dynamics."""
    from rpp_orchestrator.component_storage import ComponentParameterStore
    folder = tmp_path / "part"
    shutil.copytree(_part_dir(name), folder)
    store = ComponentParameterStore()
    before = store.load(folder)
    store.save(folder, before)
    after = store.load(folder)
    assert after == before
    assert "#" not in (folder / "params" / "parameters.py").read_text().split("class ComponentParameters")[1]
    cls = _propulsor_class()
    original = vc.context_of(cls, before, reference=False)
    original.initialize()
    saved = vc.context_of(cls, after, reference=False)
    saved.initialize()
    from more_common.casadi_graph import RppCasadiGraph
    a, b = RppCasadiGraph(original.get_instance().graph()), RppCasadiGraph(saved.get_instance().graph())
    g = rng()
    n_in = sum(d.size for d in original.get_instance().graph().inputDescription)
    n_state = sum(d.size for d in original.get_instance().graph().stateDescription)
    for _ in range(10):
        x, u = g.uniform(-1, 1, n_state), g.uniform(-1, 1, n_in)
        assert np.array_equal(np.asarray(a.output(x, u)), np.asarray(b.output(x, u)))
