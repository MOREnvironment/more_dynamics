"""The ``site`` slot of the vehicle types: the water density, the gravity and the kinematic viscosity are one
component, ``Site``, in the slot of every vehicle (no vehicle declares them), private to a vehicle by default and
linked into several vehicles by a scenario that puts them in the same water.

The linking is rpp's one sharing mechanism: ``ComponentContextBuilder.resolve_component`` redirects a linked
subcomponent record to the original component's own folder
(``rpp_py/context_builder.py`` 190-206, ``rpp_py/data_manager.py`` 123-127), so every parent that links the
component reads the same ``params/parameters.py``. The scenario tree is written here in a temporary folder in the
shape rpp-orchestrator saves (a linked subcomponent is a record with ``LinkedComponentId``); it is never written beside
the library's own ``.rppws``.

Author:    Enio Krizman
Date:      2026-10-10
"""
import json
import shutil
import uuid
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc

LIBRARY = Path(__file__).resolve().parents[2]
LIBRARY_NAME = "more_dynamics"
REST = np.zeros(12)


def _write(folder, description, params=None):
    (folder / "params").mkdir(parents=True, exist_ok=True)
    (folder / "description.json").write_text(json.dumps(description, indent=4))
    (folder / "callbacks.py").write_text("from __future__ import annotations\n")
    body = "".join(f"    {k} = {v!r}\n" for k, v in (params or {}).items()) or "    pass\n"
    (folder / "params" / "parameters.py").write_text(
        f"from __future__ import annotations\n\n\nclass ComponentParameters:\n{body}")


def _record(part_id, name, plugin_type, plugin_name, spec=None, subcomponents=None):
    return {"Id": part_id, "Name": name, "PluginType": plugin_type, "PluginName": plugin_name, "Library": LIBRARY_NAME,
            "SubcomponentSpec": spec or {}, "Subcomponents": subcomponents or {}, "ParentComponentInfo": None}


def _entry(part_id, plugin_type, plugin_name, slot, linked=False):
    return {"Id": part_id, "PluginType": plugin_type, "PluginName": plugin_name, "SlotName": slot,
            "Library": LIBRARY_NAME, "IsLinked": linked}


def _scenario(tmp_path, site_id, vehicles):
    """A ``.rppws`` tree: the top-level ``Site`` ``site_id`` (parameters ``site_params``), and one ``Catamaran`` per
    entry of ``vehicles`` (``name -> None`` links the shared site, ``name -> params`` gives the vehicle its own
    ``Site`` with those parameters). Returns ``({name: vehicle id}, description path)``."""
    (tmp_path / "scripts" / "vehicles").mkdir(parents=True)
    shutil.copy(LIBRARY / "scripts" / "vehicles" / "vehicle_simulation.py",
                tmp_path / "scripts" / "vehicles" / "vehicle_simulation.py")
    parts = tmp_path / ".rppws" / "parts"
    spec = {"actuators": "List[more_dynamics::ForceProducer]", "site": "more_dynamics::SiteModel",
            "hydrostatics": "more_dynamics::HydrostaticsModel", "hydrodynamics": "more_dynamics::HydrodynamicsModel"}
    ids = {}
    for name, own_site in vehicles.items():
        vehicle_id, restoring_id, loads_id, child_site_id = (str(uuid.uuid4()) for _ in range(4))
        ids[name] = vehicle_id
        folder = parts / "more_dynamics__catamaran" / vehicle_id
        children = folder / "subcomponents"
        _write(children / restoring_id, _record(restoring_id, f"{name} hydrostatics", "more_dynamics::HydrostaticsModel",
                                                "more_dynamics::SurfaceRestoring"), vc.reference_part("Otter hydrostatics"))
        _write(children / loads_id, _record(loads_id, f"{name} hydrodynamics", "more_dynamics::HydrodynamicsModel",
                                            "more_dynamics::SurfaceHullLoads"), vc.reference_part("Otter hull loads"))
        if own_site is None:  # a linked record: the shared site's own folder is read
            (children / child_site_id).mkdir(parents=True)
            (children / child_site_id / "description.json").write_text(json.dumps({
                "Id": child_site_id, "Name": "shared site", "LinkedComponentId": site_id,
                "LinkedComponentWorkspace": LIBRARY_NAME,
                "ParentComponentInfo": {"Id": vehicle_id, "PluginType": "more_dynamics::VehicleModel3D",
                                        "PluginName": "more_dynamics::Catamaran", "SlotName": "site",
                                        "Library": LIBRARY_NAME, "IsLinked": False}}, indent=4))
            site_entry = _entry(child_site_id, "more_dynamics::SiteModel", "more_dynamics::Site", "site", linked=True)
        else:
            _write(children / child_site_id, _record(child_site_id, f"{name} site", "more_dynamics::SiteModel",
                                                     "more_dynamics::Site"), own_site)
            site_entry = _entry(child_site_id, "more_dynamics::SiteModel", "more_dynamics::Site", "site")
        _write(folder, _record(vehicle_id, name, "more_dynamics::VehicleModel3D", "more_dynamics::Catamaran", spec, {
            "site": site_entry,
            "hydrostatics": _entry(restoring_id, "more_dynamics::HydrostaticsModel", "more_dynamics::SurfaceRestoring",
                                   "hydrostatics"),
            "hydrodynamics": _entry(loads_id, "more_dynamics::HydrodynamicsModel", "more_dynamics::SurfaceHullLoads",
                                    "hydrodynamics")}), {**vc.reference_part("Otter (MSS)"), "current_form": "none"})
    description = {
        "ScriptPath": "scripts/vehicles/vehicle_simulation.py", "Language": "python",
        "Configurations": {"scenario": {"Description": "vehicles in one water (test)", "Components": {
            "vessels": [{"Id": i, "PluginName": "more_dynamics::Catamaran"} for i in ids.values()]}}},
        "ActiveConfiguration": "scenario", "Spec": {"vessels": "List[more_dynamics::VehicleModel3D]"}}
    path = tmp_path / ".rppws" / "script_descriptions" / "scenario.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(description, indent=4))
    return ids, path


def _build(tmp_path, description_path):
    try:
        from rpp_py.context_builder import ComponentContextBuilder
        from rpp_py.data_manager import DataManager
        from rpp_plugin_registrator.library_manager import LibraryManager
    except Exception as exc:  # pragma: no cover - environment without rpp
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
    vc.plugins()
    data_manager = DataManager(library_manager=LibraryManager(), workspace_path=str(tmp_path))
    context = ComponentContextBuilder(data_manager=data_manager).build_script_from_description_path(
        str(description_path), configuration="scenario")
    context.initialize()
    return context.get_component("vessels")


def _site_values(vessel):
    signals = vessel.signals(REST, [])
    return tuple(float(signals[name][0]) for name in ("water_density", "gravity", "kinematic_viscosity"))


def test_a_site_linked_into_two_vehicles_is_edited_in_one_place(tmp_path):
    """Two vehicles link one ``Site``, a third has its own: editing the shared site alone moves the two linked
    vehicles identically and leaves the third untouched."""
    site_id = str(uuid.uuid4())
    shared = tmp_path / ".rppws" / "parts" / "more_dynamics__site" / site_id
    _write(shared, _record(site_id, "water of the scenario", "more_dynamics::SiteModel", "more_dynamics::Site"),
           {"water_density": 1025.0, "gravity": 9.81, "kinematic_viscosity": 1e-06})
    ids, description = _scenario(tmp_path, site_id, {"first": None, "second": None, "alone": {
        "water_density": 1025.0, "gravity": 9.81, "kinematic_viscosity": 1e-06}})
    first, second, alone = _build(tmp_path, description)
    assert _site_values(first) == _site_values(second) == _site_values(alone) == (1025.0, 9.81, 1e-06)

    _write(shared, _record(site_id, "water of the scenario", "more_dynamics::SiteModel", "more_dynamics::Site"),
           {"water_density": 1030.0, "gravity": 9.80, "kinematic_viscosity": 1.2e-06})
    first, second, alone = _build(tmp_path, description)
    assert _site_values(first) == _site_values(second) == (1030.0, 9.80, 1.2e-06)
    assert _site_values(alone) == (1025.0, 9.81, 1e-06)


def test_a_vehicle_type_declares_none_of_the_site_quantities():
    p = vc.plugins()
    for cls in (p.SpheroidAuv, p.Catamaran, p.Monohull, p.DerivativeAuv):
        assert "site" in cls.COMPONENTS and cls.COMPONENTS["site"] == "more_dynamics::SiteModel", cls
        assert not {d.name for d in cls.PARAMETERS} & set(vc.SITE_NAMES), cls


def test_the_site_gives_its_three_quantities_in_both_forms():
    p = vc.plugins()
    for form, values in (("given", {"gravity": 9.5, "water_density": 1000.0, "kinematic_viscosity": 1.1e-06}),
                         ("latitude", {"latitude": 0.0, "water_density": 1000.0, "kinematic_viscosity": 1.1e-06})):
        site = vc.child((p.Site, {"site_form": form, **values}))
        site.initialize()
        from more_common.casadi_graph import RppCasadiGraph
        graph = RppCasadiGraph(site.get_instance().graph())
        assert [d.name for d in graph.payload.outputDescription] == ["water_density", "gravity", "kinematic_viscosity"]
        out = np.asarray(graph.output(np.zeros(0), np.zeros(0))).ravel()
        assert out[0] == 1000.0 and out[2] == 1.1e-06
        assert out[1] == (9.5 if form == "given" else pytest.approx(9.7803253359, rel=1e-13))


def test_an_unknown_site_form_is_refused_naming_the_choices():
    p = vc.plugins()
    with pytest.raises(p.CompositionError, match="site_form must be one of"):
        vc.child((p.Site, {"site_form": "geodetic"})).initialize()


def test_the_committed_shared_site_scenario_links_one_site_into_both_boats():
    """``vehicle_simulation.json`` configuration ``shared_site``: two boats of the Otter's parameters link the one
    ``Scenario site``; each builds to the frozen Otter vehicle, and their site is the one component's."""
    try:
        from rpp_py.context_builder import ComponentContextBuilder
        from rpp_py.data_manager import DataManager
        from rpp_plugin_registrator.library_manager import LibraryManager
    except Exception as exc:  # pragma: no cover - environment without rpp
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
    vc.plugins()
    description = LIBRARY / ".rppws" / "script_descriptions" / "vehicle_simulation.json"
    data_manager = DataManager(library_manager=LibraryManager(), workspace_path=str(LIBRARY))
    context = ComponentContextBuilder(data_manager=data_manager).build_script_from_description_path(
        str(description), configuration="shared_site")
    context.initialize()
    boats = context.get_component("vessels")
    assert len(boats) == 2
    linked = [entry for entry in (json.loads(p.read_text()) for p in (LIBRARY / ".rppws" / "parts").rglob("description.json")
                                  if "LinkedComponentId" in p.read_text())
              if entry["ParentComponentInfo"]["SlotName"] == "site"]  # the site links, not the actuator links
    assert len(linked) == 2 and len({entry["LinkedComponentId"] for entry in linked}) == 1
    reference = vc.graph_of(vc.otter()).step
    x = np.r_[np.zeros(6), 1.0, 0.1, 0.0, 0.0, 0.0, 0.2]
    for boat in boats:
        assert _site_values(boat) == (1025.0, 9.81, 1e-06)
        assert np.array_equal(np.asarray(vc.graph_of(boat).step(x, [])), np.asarray(reference(x, [])))
