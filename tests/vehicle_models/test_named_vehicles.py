"""The named vehicles of the library's own ``.rppws`` (REMUS 100 (MSS), Otter
(MSS) and Grethe), built by rpp's own ``ComponentContextBuilder`` from
``.rppws/script_descriptions/vehicle_simulation.json``. REMUS 100 (MSS) and
Otter (MSS) equal the plugins' declared defaults, whose gates are in
``test_spheroid_auv.py`` and ``test_catamaran.py``; Grethe's own physics gates
are in ``test_monohull_grethe.py`` (U8a/U8b) -- here only the part itself is
checked: it builds, and every parameter file of every part names only
parameters its plugin declares, so rpp does not drop a value silently. These
three parts are the reference vehicles: a vehicle with other numbers is made
in the GUI as a copy, not by editing these.

Author:    Enio Krizman
Date:      2026-10-09
"""

import importlib
import json
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc

LIBRARY = Path(__file__).resolve().parents[2]
RPPWS = LIBRARY / ".rppws"
DESCRIPTION = RPPWS / "script_descriptions" / "vehicle_simulation.json"
NAMED = {"remus100": "REMUS 100 (MSS)", "otter": "Otter (MSS)", "grethe": "Grethe"}


def _build(configuration):
    vc.plugins()
    try:
        from rpp_py.data_manager import DataManager
        from rpp_plugin_registrator.library_manager import LibraryManager
        data_manager = DataManager(library_manager=LibraryManager(), workspace_path=str(LIBRARY))
        from rpp_py.context_builder import ComponentContextBuilder
        context = ComponentContextBuilder(data_manager=data_manager).build_script_from_description_path(
            str(DESCRIPTION), configuration=configuration)
        context.initialize()
    except Exception as exc:  # more_dynamics is not registered here
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
    return context.get_component("vessels")[0]


def _states(n, seed):
    g = np.random.default_rng(seed)
    for _ in range(n):
        yield np.concatenate([g.uniform(-0.3, 0.3, 6), g.uniform(-1.0, 2.5, 3), g.uniform(-0.3, 0.3, 3)])


def test_the_script_description_lists_one_configuration_per_named_vehicle():
    description = json.loads(DESCRIPTION.read_text())
    assert set(description["Configurations"]) == set(NAMED) | {"shared_site"}  # the last is a scenario, not a vehicle
    assert description["Spec"] == {"vessels": "List[more_dynamics::VehicleModel3D]"}
    assert (LIBRARY / description["ScriptPath"]).is_file()
    for key, configuration in description["Configurations"].items():
        if key == "shared_site":
            continue  # two boats linking one site: test_site_slot.py
        assert len(configuration["Components"]["vessels"]) == 1
        part_id = configuration["Components"]["vessels"][0]["Id"]
        # the part's own description.json carries the E-106 name, found by id (rule 20, R10)
        matches = [p for p in RPPWS.joinpath("parts").rglob("description.json")
                  if json.loads(p.read_text()).get("Id") == part_id]
        assert len(matches) == 1, (key, part_id)
        assert json.loads(matches[0].read_text())["Name"] == NAMED[key]


def test_remus_100_from_the_workspace_equals_the_plugins_on_the_part_values():
    """The workspace part built by rpp's builder equals the plugins built from the part's values and selectors
    (still water; the part's own site: gravity from the latitude, 1026 kg/m^3; MSS's Munk-removed Coriolis)."""
    p = vc.plugins()
    built = _build("remus100")
    site = vc.reference_part("REMUS 100 (MSS) site")
    reference = vc.build_context(p.SpheroidAuv, {"current_form": "none", "coriolis_form": "munk_couplings_removed"},
                                 site=(p.Site, {"site_form": "latitude", "latitude": site["latitude"],
                                                "water_density": site["water_density"]}),
                                 hydrostatics=p.SubmergedRestoring, hydrodynamics=p.AuvHullLoads,
                                 actuators=[p.FinPairsDeflectionOnly, p.Propeller])
    assert type(built).__name__ == "SpheroidAuv"
    assert vc.input_names(built) == vc.input_names(reference) == ["fin_deflection_command", "shaft_speed_command"]
    a, b = vc.graph_of(built).step, vc.graph_of(reference).step
    rng = np.random.default_rng(5)
    for x in _states(30, 3):
        u = rng.uniform(-0.3, 0.3, 2)
        u = np.r_[u, rng.uniform(0.0, 1500.0)]
        assert np.array_equal(np.asarray(a(x, u)), np.asarray(b(x, u)))


def test_otter_from_the_workspace_equals_the_frozen_otter_set():
    built = _build("otter")
    reference = vc.otter()
    assert type(built).__name__ == "Catamaran"
    assert vc.input_names(built) == []
    a, b = vc.graph_of(built).step, vc.graph_of(reference).step
    for x in _states(30, 4):
        assert np.array_equal(np.asarray(a(x, [])), np.asarray(b(x, [])))


def _own_parameters(part_id):
    """``{name: value}`` of one part's own ``params/parameters.py``, found by id (as
    ``test_every_parameter_file_names_only_parameters_its_plugin_declares`` reads them)."""
    for description in RPPWS.joinpath("parts").rglob("description.json"):
        if json.loads(description.read_text()).get("Id") == part_id:
            namespace = {}
            exec((description.parent / "params" / "parameters.py").read_text(), namespace)
            cp = namespace["ComponentParameters"]
            return {k: getattr(cp, k) for k in vars(cp) if not k.startswith("_")}
    raise AssertionError(f"no part with id {part_id!r}")


def test_grethe_from_the_workspace_builds_with_its_own_selectors():
    """Structural only (physics gates are U8a/U8b's ``test_monohull_grethe.py``): the part
    builds as a ``Monohull`` with no open actuator input, and its two children' own parameter
    files (not the plugins' defaults) carry E-101's selectors."""
    built = _build("grethe")
    assert type(built).__name__ == "Monohull"
    assert vc.input_names(built) == []
    description = json.loads(DESCRIPTION.read_text())
    vessel = description["Configurations"]["grethe"]["Components"]["vessels"][0]
    root = next(p for p in RPPWS.joinpath("parts").rglob("description.json")
               if json.loads(p.read_text()).get("Id") == vessel["Id"])
    record = json.loads(root.read_text())
    hydrostatics = _own_parameters(record["Subcomponents"]["hydrostatics"]["Id"])
    hydrodynamics = _own_parameters(record["Subcomponents"]["hydrodynamics"]["Id"])
    root_params = _own_parameters(vessel["Id"])
    assert root_params["added_mass_mass_basis"] == "displaced"
    assert hydrostatics["wetted_surface_method"] == "regression_table"
    assert hydrodynamics["manoeuvring_damping"] == "linear_coupled"
    assert hydrodynamics["surge_resistance"] == "ittc_residual"
    assert hydrodynamics["sway_yaw_derivatives"] == [-251.875, 0.0, 0.0, -229.86112500000002]


def test_every_parameter_file_names_only_parameters_its_plugin_declares():
    vc.plugins()
    registry = json.loads((LIBRARY / "plugins.json").read_text())
    checked = 0
    for description in sorted(RPPWS.joinpath("parts").rglob("description.json")):
        record = json.loads(description.read_text())
        if "LinkedComponentId" in record:
            continue  # a link to a component whose own folder is checked
        name = record["PluginName"].split("::")[1]
        if name in ("HullVessel", "JetNozzle", "LinearSurfaceHydrostatics", "LinearSurfaceHydrodynamics"):
            continue  # Luka's parts keep his files
        path = next((e["Path"] for e in registry["Plugins"] if e["Name"] == name), None)
        if path is None or record["Library"] != "more_dynamics":
            continue  # a part of another library in Luka's tree
        cls = getattr(importlib.import_module(path[:-3].replace("/", ".")), name)
        namespace = {}
        exec((description.parent / "params" / "parameters.py").read_text(), namespace)
        values = {k for k in vars(namespace["ComponentParameters"]) if not k.startswith("_")}
        assert values <= {d.name for d in cls.PARAMETERS}, (name, values - {d.name for d in cls.PARAMETERS})
        checked += 1
    assert checked >= 8
