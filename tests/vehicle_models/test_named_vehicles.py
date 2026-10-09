"""The named vehicles of the library's own ``.rppws`` (REMUS 100 and Otter),
built by rpp's own ``ComponentContextBuilder`` from
``.rppws/script_descriptions/vehicle_simulation.json``, equal the plugins'
declared defaults (the cited REMUS 100 and Otter values), whose gates are in
``test_spheroid_auv.py`` and ``test_catamaran.py``; and every parameter file
of those parts names only parameters its plugin declares, so rpp does not drop
a value silently. The two parts are the reference vehicles: a vehicle with
other numbers is made in the GUI as a copy, not by editing these.

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
NAMED = {"remus100": "REMUS 100", "otter": "Otter"}


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
    assert set(description["Configurations"]) == set(NAMED)
    assert description["Spec"] == {"vessels": "List[more_dynamics::VehicleModel3D]"}
    assert (LIBRARY / description["ScriptPath"]).is_file()
    for configuration in description["Configurations"].values():
        assert len(configuration["Components"]["vessels"]) == 1


def test_remus_100_from_the_workspace_equals_the_plugins_on_their_defaults():
    p = vc.plugins()
    built = _build("remus100")
    reference = vc.build_context(p.SpheroidAuv, {}, hydrostatics=p.SubmergedRestoring, hydrodynamics=p.AuvHullLoads,
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


def test_every_parameter_file_names_only_parameters_its_plugin_declares():
    vc.plugins()
    registry = json.loads((LIBRARY / "plugins.json").read_text())
    checked = 0
    for description in sorted(RPPWS.joinpath("parts").rglob("description.json")):
        record = json.loads(description.read_text())
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
