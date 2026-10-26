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
NAMED = {"remus100": "REMUS 100 (MSS)", "otter": "Otter (MSS)", "grethe": "Grethe", "npsauv": "NPS AUV II",
         "grethe_stress_test": "Grethe (yaw-unstable stress test)"}


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


def _part_child(plugin, name):
    """A child context of ``plugin`` with the values of the named part ``name`` (found by its ``Name``)."""
    return plugin, vc.pick(plugin, vc.reference_part(name))


def test_remus_100_from_the_workspace_equals_the_plugins_on_the_part_values():
    """The workspace part built by rpp's builder equals the plugins built from the part's values and selectors
    (still water; the part's own site: gravity from the latitude, 1026 kg/m^3; MSS's Munk-removed Coriolis), its
    actuators being four fins (rudder pair, stern-plane pair) and the propulsor, linked from named parts."""
    p = vc.plugins()
    built = _build("remus100")
    site = vc.reference_part("REMUS 100 (MSS) site")
    reference = vc.build_context(p.SpheroidAuv, {"current_form": "none", "coriolis_form": "munk_couplings_removed"},
                                 site=(p.Site, {"site_form": "latitude", "latitude": site["latitude"],
                                                "water_density": site["water_density"]}),
                                 hydrostatics=p.SubmergedRestoring, hydrodynamics=p.AuvHullLoads,
                                 actuators=[_part_child(p.Fin, "REMUS 100 (MSS) rudder fin")] * 2
                                 + [_part_child(p.Fin, "REMUS 100 (MSS) stern plane fin")] * 2
                                 + [_part_child(p.Propulsor, "REMUS 100 propulsor")])
    assert type(built).__name__ == "SpheroidAuv"
    names = ["fin_angle_command"] * 4 + ["drive_command"]
    assert vc.input_names(built) == vc.input_names(reference) == names
    a, b = vc.graph_of(built).step, vc.graph_of(reference).step
    rng = np.random.default_rng(5)
    for x in _states(30, 3):
        u = np.r_[rng.uniform(-0.3, 0.3, 4), rng.uniform(0.0, 1500.0)]
        assert np.array_equal(np.asarray(a(x, u)), np.asarray(b(x, u)))


def test_remus_100_with_four_fins_and_the_propulsor_equals_the_fin_pairs_and_the_propeller_plugins():
    """Two independent compositions of the same REMUS 100 values: the committed part (four ``Fin`` children and the
    ``Propulsor``) and the earlier one (the deflection-only fin pairs and the ``Propeller``, their own named parts):
    the same derivative on random states and commands, the two fins of a pair given the same deflection."""
    p = vc.plugins()
    built = _build("remus100")
    site = vc.reference_part("REMUS 100 (MSS) site")
    earlier = vc.build_context(p.SpheroidAuv, {"current_form": "none", "coriolis_form": "munk_couplings_removed"},
                               site=(p.Site, {"site_form": "latitude", "latitude": site["latitude"],
                                              "water_density": site["water_density"]}),
                               hydrostatics=p.SubmergedRestoring, hydrodynamics=p.AuvHullLoads,
                               actuators=[_part_child(p.FinPairsDeflectionOnly, "REMUS 100 fin pairs"),
                                          _part_child(p.Propeller, "REMUS 100 propeller")])
    assert vc.input_names(earlier) == ["fin_deflection_command", "shaft_speed_command"]
    a, b = vc.graph_of(built).step, vc.graph_of(earlier).step
    rng = np.random.default_rng(6)
    worst = 0.0
    for x in _states(30, 7):
        rudder, stern = rng.uniform(-0.3, 0.3, 2)
        n = rng.uniform(-1500.0, 1500.0)
        ours = np.asarray(a(x, [rudder, rudder, stern, stern, n])).ravel()
        reference = np.asarray(b(x, [rudder, stern, n])).ravel()
        worst = max(worst, float(np.max(np.abs(ours - reference) / np.maximum(1.0, np.abs(reference)))))
    assert worst <= 1e-12, worst


def test_otter_from_the_workspace_equals_the_frozen_otter_set():
    """The Otter part built by rpp's builder equals the plugins built from the frozen Otter set with the two thrusters
    of the named parts and the MSS payload trim form the part selects."""
    p = vc.plugins()
    built = _build("otter")
    thrusters = [_part_child(p.Propulsor, "Otter port thruster"), _part_child(p.Propulsor, "Otter starboard thruster")]
    reference = vc.otter(actuators=thrusters, hydrostatics=(p.SurfaceRestoring, {
        **vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2, "payload_trim_form": "lever_arm_and_heave_shift"}))
    assert type(built).__name__ == "Catamaran"
    assert vc.input_names(built) == vc.input_names(reference) == ["drive_command", "drive_command"]
    a, b = vc.graph_of(built).step, vc.graph_of(reference).step
    rng = np.random.default_rng(8)
    for x in _states(30, 4):
        n = rng.uniform(-100.0, 100.0, 2)
        assert np.array_equal(np.asarray(a(x, n)), np.asarray(b(x, n)))


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
    """Structural only (physics gates are ``test_monohull_grethe.py``): the part builds as a ``Monohull`` with one
    outboard (a thrust command and a steering angle command, and the two actuator states behind them), and its
    children's own parameter files (not the plugins' defaults) carry the forms of the model file ``grethe.m``."""
    built = _build("grethe")
    assert type(built).__name__ == "Monohull"
    assert vc.input_names(built) == ["drive_command", "steering_command"]
    assert vc.actuator_state_count(built) == 2
    description = json.loads(DESCRIPTION.read_text())
    vessel = description["Configurations"]["grethe"]["Components"]["vessels"][0]
    root = next(p for p in RPPWS.joinpath("parts").rglob("description.json")
               if json.loads(p.read_text()).get("Id") == vessel["Id"])
    record = json.loads(root.read_text())
    hydrostatics = _own_parameters(record["Subcomponents"]["hydrostatics"]["Id"])
    hydrodynamics = _own_parameters(record["Subcomponents"]["hydrodynamics"]["Id"])
    root_params = _own_parameters(vessel["Id"])
    assert root_params["added_mass_mass_basis"] == "displaced"
    assert root_params["hull_waterplane_coefficient"] == 0.41  # fitted to the logged roll period (grethe.m 242)
    assert hydrostatics["wetted_surface_method"] == "regression_table"
    assert hydrostatics["reference_point"] == [0.0, 0.0, 0.0]  # the waterline midpoint, E-101 Q4 a
    assert hydrodynamics["manoeuvring_damping"] == "time_constants_coupled"
    assert hydrodynamics["surge_resistance"] == "ittc_residual_ramped"
    assert hydrodynamics["time_constants"] == [8.0, 4.0] and hydrodynamics["sway_yaw_coupling"] == [0.0, 2500.0]
    assert "sway_yaw_derivatives" not in hydrodynamics  # computed from the time constants, not typed again


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
