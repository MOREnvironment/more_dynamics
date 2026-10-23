"""Vehicle gate helpers: build a vehicle type of ``more_dynamics`` with rpp's
own ``ComponentContext`` (as Luka's and Matko's tests build a plugin), its
hydrostatics, hydrodynamics and actuators as child contexts, and read its
payload. A vehicle is one plugin with its own parameters; REMUS 100 and the
Otter differ only in the parameter values and the children handed to it,
the ``site`` child (the water density, the gravity and the viscosity) among
them. The
reference numbers live under ``tests/vehicle_models/data``.

The plugins subclass rpp's generated plugin types, so the tests need rpp
installed and the generated types available (``rpp library register
./more_dynamics --link``, ``workspace/`` as the working folder); without it
they skip with that message and the model-layer gates run as usual.

Author:    Enio Krizman
Date:      2026-10-09
"""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DATA = Path(__file__).resolve().parent / "data"
REGISTER_MESSAGE = ("rpp's generated plugin types for more_dynamics are not available: run "
                    "`rpp library register ./more_dynamics --link` (workspace/ as the working folder)")
KINEMATIC_VISCOSITY = 1e-6  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s
CURRENT_NAMES = ("current_speed", "current_direction", "current_vertical_speed")
SITE_NAMES = ("site_form", "latitude", "gravity", "water_density", "kinematic_viscosity")  # the site's, not the vehicle's

REMUS_FILE = json.loads((DATA / "remus100" / "remus100_parameters_consistent.json").read_text())
REMUS = {name: entry["value"] for name, entry in REMUS_FILE["parameters"].items()}
OTTER = json.loads((DATA / "otter" / "otter_mss_99fdf4d.json").read_text())

_PLUGINS = {}

RPPWS = Path(__file__).resolve().parents[2] / ".rppws"
# the named part whose values a plugin is built from when a test gives none (the plugins carry no vehicle
# numbers: every vehicle, device and site value comes from a part); a plugin with no part reads its frozen data file
REFERENCE_PARTS = {
    "SpheroidAuv": "REMUS 100 (MSS)", "Monohull": "Grethe", "Catamaran": "Otter (MSS)",
    "SubmergedRestoring": "REMUS 100 hydrostatics", "SurfaceRestoring": "Otter hydrostatics",
    "AuvHullLoads": "REMUS 100 hull loads", "SurfaceHullLoads": "Otter hull loads",
    "FinPairsDeflectionOnly": "REMUS 100 fin pairs", "Propeller": "REMUS 100 propeller",
    "Site": "Otter (MSS) site", "Fin": "REMUS 100 fin",
}
REFERENCE_FILES = {"DerivativeAuv": "npsauv/nps_auv_ii_parameters.json",
                   "CoefficientHullLoads": "npsauv/nps_auv_ii_parameters.json"}


def reference_part(name):
    """``{parameter: value}`` of the named part ``name`` of the library's workspace (its own ``params``, read with
    rpp's own parameter loader, the part found by the ``Name`` of its ``description.json``)."""
    from rpp_orchestrator.component_storage import ComponentParameterStore
    found = [d.parent for d in RPPWS.joinpath("parts").rglob("description.json")
             if json.loads(d.read_text()).get("Name") == name]
    assert len(found) == 1, f"named part {name!r}: found {len(found)} under {RPPWS}"
    return ComponentParameterStore().load(found[0])


def reference_values(cls):
    """The values a plugin class takes where its default is unset and a test gives none: its reference part's, or
    its frozen data file's."""
    from more_transformations.more_casadi_transformations import is_unset
    if cls.__name__ in REFERENCE_PARTS:
        source = reference_part(REFERENCE_PARTS[cls.__name__])
    elif cls.__name__ in REFERENCE_FILES:
        entries = json.loads((DATA / REFERENCE_FILES[cls.__name__]).read_text())["parameters"]
        source = {name: entry["value"] for name, entry in entries.items()}
    else:
        return {}
    return {d.name: source[d.name] for d in cls.PARAMETERS if is_unset(d.default_value) and d.name in source}


def plugins():
    """The plugin classes of the library, importable once rpp's generated types are; skips otherwise."""
    if not _PLUGINS:
        try:
            from rpp_py.data_manager import DataManager
            DataManager()  # before the generated modules are imported
            from more_dynamics.plugins.vehicle_models.spheroid_auv import SpheroidAuv
            from more_dynamics.plugins.vehicle_models.monohull import Monohull
            from more_dynamics.plugins.vehicle_models.catamaran import Catamaran
            from more_dynamics.plugins.vehicle_models.derivative_auv import DerivativeAuv
            from more_dynamics.plugins.site.site import Site
            from more_dynamics.plugins.hydrostatics.submerged_restoring import SubmergedRestoring
            from more_dynamics.plugins.hydrostatics.surface_restoring import SurfaceRestoring
            from more_dynamics.plugins.hydrodynamics.auv_hull_loads import AuvHullLoads
            from more_dynamics.plugins.hydrodynamics.surface_hull_loads import SurfaceHullLoads
            from more_dynamics.plugins.hydrodynamics.coefficient_hull_loads import CoefficientHullLoads
            from more_dynamics.plugins.force_producers.fin import Fin
            from more_dynamics.plugins.force_producers.fin_pairs_deflection_only import FinPairsDeflectionOnly
            from more_dynamics.plugins.force_producers.propeller import Propeller
            from more_dynamics.plugins.force_producers.prescribed_wrench import PrescribedWrench
            from more_dynamics.plugins.shared.vehicle_graph import CompositionError
            from more_dynamics.plugins.shared.payload_io import UnsetValues
        except Exception as exc:  # no rpp, or the generated types are missing
            pytest.skip(f"{REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
        _PLUGINS.update(locals())
        _PLUGINS.pop("DataManager", None)
        _PLUGINS.pop("exc", None)
    return SimpleNamespace(**_PLUGINS)


def context_of(cls, params=None, children=None, spec=None, reference=True):
    """A ``ComponentContext`` of ``cls`` with its declared defaults, its reference part's values where the default is
    unset (``reference_values``), ``params`` over them, and child contexts by slot (``{slot: [context, ...]}``)."""
    from rpp_py.context import ComponentContext
    from rpp_py.parameter_handler import ParameterHandler
    resolved = ParameterHandler.resolve_params(cls.PARAMETERS, {})
    unknown = set(params or {}) - {d.name for d in cls.PARAMETERS}
    assert not unknown, f"{cls.__name__} declares no parameter {sorted(unknown)}"
    if reference:
        resolved.params.update(reference_values(cls))
    resolved.params.update(params or {})
    return ComponentContext(instance=cls(), params=resolved, subcomponents=children or {},
                            spec=spec if spec is not None else getattr(cls, "COMPONENTS", {}))


def child(item):
    """``(class, parameters)``, a class, or an already made context -> its context."""
    if not isinstance(item, (type, tuple)):
        return item
    cls, params = item if isinstance(item, tuple) else (item, None)
    return context_of(cls, params)


def build_context(vehicle_class, params=None, *, hydrostatics, hydrodynamics, actuators=(), site=None, reference=True):
    """The initialised context of a vehicle: ``site``, ``hydrostatics`` and ``hydrodynamics`` are a class or a
    ``(class, parameters)`` (``site`` default: the ``Site`` plugin with its own defaults), ``actuators`` a list of
    them. Returns the vehicle plugin."""
    children = {"site": [child(plugins().Site if site is None else site)],
                "hydrostatics": [child(hydrostatics)], "hydrodynamics": [child(hydrodynamics)],
                "actuators": [child(a) for a in actuators]}
    context = context_of(vehicle_class, params, children, reference=reference)
    context.initialize()
    return context.get_instance()


def pick(cls, values, names=None):
    """The entries of ``values`` that ``cls`` declares as parameters (only ``names`` when given)."""
    declared = {d.name for d in cls.PARAMETERS}
    return {n: values[n] for n in (names or declared) if n in declared and n in values}


def with_site(site_values, vehicle_params):
    """``(site parameters, vehicle parameters)``: the names of ``SITE_NAMES`` a test passes among a vehicle's
    parameters belong to the site child and are handed to it over ``site_values``."""
    vehicle_params = dict(vehicle_params or {})
    return {**site_values, **{n: vehicle_params.pop(n) for n in SITE_NAMES if n in vehicle_params}}, vehicle_params


def remus_values(values=None):
    """The REMUS 100 set as the plugins' parameter names: the frozen file plus the site's latitude and the one
    kinematic viscosity."""
    v = dict(REMUS if values is None else values)
    v["diameter"] = v["beam"]
    v.setdefault("body_mass", REMUS_FILE["one_value"]["mass"])  # remus100.m:3, the given mass; the density follows
    v["latitude"] = REMUS_FILE["gravity"]["latitude_rad"]
    v["kinematic_viscosity"] = KINEMATIC_VISCOSITY
    return v


def remus_site(values=None):
    """The ``Site`` parameters of the REMUS 100 set: gravity from the latitude, its water density and viscosity."""
    v = remus_values(values)
    return {"site_form": "latitude", "latitude": v["latitude"], "water_density": v["water_density"],
            "kinematic_viscosity": v["kinematic_viscosity"]}


def remus(*, full=False, coriolis="kirchhoff_full", current="yaw_rate_terms", values=None, actuators=None,
          vehicle_params=None):
    """REMUS 100: the hull alone (a prescribed wrench, ``full=False``) or with the fin pairs and the propeller.
    ``values`` replaces the frozen parameter set (a perturbed set in a test)."""
    p = plugins()
    v = remus_values(values)
    site_p, vehicle_params = with_site(remus_site(values), vehicle_params)
    vehicle_p = {**pick(p.SpheroidAuv, v), **({"coriolis_form": coriolis} if coriolis else {}),
                 "current_form": current, **vehicle_params}
    if actuators is None:
        actuators = ([(p.FinPairsDeflectionOnly, pick(p.FinPairsDeflectionOnly, v)), (p.Propeller, pick(p.Propeller, v))]
                     if full else [p.PrescribedWrench])
    return build_context(p.SpheroidAuv, vehicle_p, site=(p.Site, site_p),
                         hydrostatics=(p.SubmergedRestoring, pick(p.SubmergedRestoring, v, ("center_of_buoyancy",))),
                         hydrodynamics=(p.AuvHullLoads, pick(p.AuvHullLoads, v)), actuators=actuators)


def otter(*, current="none", hydrostatics=None, hydrodynamics=None, actuators=(), vehicle_params=None):
    """The Otter as a catamaran: the frozen MSS set, no force producer unless ``actuators`` gives some."""
    p = plugins()
    o = dict(OTTER)
    o["kinematic_viscosity"] = KINEMATIC_VISCOSITY
    o["added_mass_form"] = "scaled_derivatives_loaded_surge"   # otter.m 162 at 99fdf4d: A11 from m + mP
    o["manoeuvring_damping"] = "time_constants_surge_blend"    # otter.m 211-212, 245-246 at 99fdf4d
    site_p, vehicle_params = with_site(
        {"site_form": "given", "gravity": o["gravity"], "water_density": o["water_density"],
         "kinematic_viscosity": o["kinematic_viscosity"]}, vehicle_params)
    vehicle_p = {**pick(p.Catamaran, o), "current_form": current, **vehicle_params}
    return build_context(
        p.Catamaran, vehicle_p, site=(p.Site, site_p),
        hydrostatics=hydrostatics or (p.SurfaceRestoring, {**pick(p.SurfaceRestoring, o), "hull_count": 2}),
        hydrodynamics=hydrodynamics or (p.SurfaceHullLoads, pick(p.SurfaceHullLoads, o)), actuators=actuators)


def graph_of(vessel):
    from more_common.casadi_graph import RppCasadiGraph
    return RppCasadiGraph(vessel.graph())


def input_names(vessel):
    return [d.name for d in vessel.getInputDescriptions()]


def input_size(vessel):
    """Number of stacked input values of the vehicle."""
    return sum(d.size for d in vessel.getInputDescriptions())


def matrix(vector):
    """A payload matrix: its column-major ``vec`` back to 6x6."""
    return np.asarray(vector).reshape(6, 6, order="F")


def dynamics_map(vessel, n):
    """``dynamics`` mapped over ``n`` columns: ``f(x (12 x n), u (m x n)) -> xdot (12 x n)``."""
    return graph_of(vessel).step.map(n)


def actuator_state_count(vessel):
    """Number of states the vehicle's actuators carry: the graph stacks them, in list order, before the 12 vessel
    states (``vehicle_graph.py``), so it is the graph's state size less 12."""
    return graph_of(vessel).num_states - 12


def mss_to_ours(x, vessel=None):
    """MSS state [nu; eta; actuator states] to ours [actuator states; eta; nu]. ``vessel`` names the vehicle whose
    actuators carry states; without it the state is the 12 vessel states. An MSS vector whose actuator columns
    are not the vehicle's states is refused, not truncated: a column with no state behind it in the vehicle
    (``npsauv.m``'s propeller shaft speed, for one) is a difference to see, not to drop."""
    x = np.asarray(x)
    count = 0 if vessel is None else actuator_state_count(vessel)
    if x.shape[-1] - 12 != count:
        raise ValueError(f"the MSS state has {x.shape[-1] - 12} actuator states, the vehicle's actuators carry {count}")
    return np.concatenate([x[..., 12:], x[..., 6:12], x[..., 0:6]], axis=-1)


def ours_to_mss(x, vessel=None):
    """Ours [actuator states; eta; nu] to the MSS order [nu; eta; actuator states] (``mss_to_ours``, reversed)."""
    x = np.asarray(x)
    count = 0 if vessel is None else actuator_state_count(vessel)
    if x.shape[-1] - 12 != count:
        raise ValueError(f"the state has {x.shape[-1] - 12} actuator states, the vehicle's actuators carry {count}")
    return np.concatenate([x[..., count + 6:], x[..., count:count + 6], x[..., :count]], axis=-1)
