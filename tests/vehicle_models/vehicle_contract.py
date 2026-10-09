"""Vehicle gate helpers: build a vehicle type of ``more_dynamics`` with rpp's
own ``ComponentContext`` (as Luka's and Matko's tests build a plugin), its
hydrostatics, hydrodynamics and actuators as child contexts, and read its
payload. A vehicle is one plugin with its own parameters; REMUS 100 and the
Otter differ only in the parameter values and the children handed to it. The
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

REMUS_FILE = json.loads((DATA / "remus100" / "remus100_parameters_consistent.json").read_text())
REMUS = {name: entry["value"] for name, entry in REMUS_FILE["parameters"].items()}
OTTER = json.loads((DATA / "otter" / "otter_mss_cc07579.json").read_text())

_PLUGINS = {}


def plugins():
    """The plugin classes of the library, importable once rpp's generated types are; skips otherwise."""
    if not _PLUGINS:
        try:
            from rpp_py.data_manager import DataManager
            DataManager()  # before the generated modules are imported
            from more_dynamics.plugins.vehicle_models.spheroid_auv import SpheroidAuv
            from more_dynamics.plugins.vehicle_models.monohull import Monohull
            from more_dynamics.plugins.vehicle_models.catamaran import Catamaran
            from more_dynamics.plugins.hydrostatics.submerged_restoring import SubmergedRestoring
            from more_dynamics.plugins.hydrostatics.surface_restoring import SurfaceRestoring
            from more_dynamics.plugins.hydrodynamics.auv_hull_loads import AuvHullLoads
            from more_dynamics.plugins.hydrodynamics.surface_hull_loads import SurfaceHullLoads
            from more_dynamics.plugins.force_producers.fin import Fin
            from more_dynamics.plugins.force_producers.fin_pairs_deflection_only import FinPairsDeflectionOnly
            from more_dynamics.plugins.force_producers.propeller import Propeller
            from more_dynamics.plugins.force_producers.prescribed_wrench import PrescribedWrench
            from more_dynamics.plugins.shared.vehicle_graph import CompositionError
        except Exception as exc:  # no rpp, or the generated types are missing
            pytest.skip(f"{REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
        _PLUGINS.update(locals())
        _PLUGINS.pop("DataManager", None)
        _PLUGINS.pop("exc", None)
    return SimpleNamespace(**_PLUGINS)


def context_of(cls, params=None, children=None, spec=None):
    """A ``ComponentContext`` of ``cls`` with its declared defaults, ``params`` over them, and child contexts by
    slot (``{slot: [context, ...]}``)."""
    from rpp_py.context import ComponentContext
    from rpp_py.parameter_handler import ParameterHandler
    resolved = ParameterHandler.resolve_params(cls.PARAMETERS, {})
    unknown = set(params or {}) - {d.name for d in cls.PARAMETERS}
    assert not unknown, f"{cls.__name__} declares no parameter {sorted(unknown)}"
    resolved.params.update(params or {})
    return ComponentContext(instance=cls(), params=resolved, subcomponents=children or {},
                            spec=spec if spec is not None else getattr(cls, "COMPONENTS", {}))


def child(item):
    """``(class, parameters)`` or a class -> its context."""
    cls, params = item if isinstance(item, tuple) else (item, None)
    return context_of(cls, params)


def build_context(vehicle_class, params=None, *, hydrostatics, hydrodynamics, actuators=()):
    """The initialised context of a vehicle: ``hydrostatics`` and ``hydrodynamics`` are a class or a
    ``(class, parameters)``, ``actuators`` a list of them. Returns the vehicle plugin."""
    children = {"hydrostatics": [child(hydrostatics)], "hydrodynamics": [child(hydrodynamics)],
                "actuators": [child(a) for a in actuators]}
    context = context_of(vehicle_class, params, children)
    context.initialize()
    return context.get_instance()


def pick(cls, values, names=None):
    """The entries of ``values`` that ``cls`` declares as parameters (only ``names`` when given)."""
    declared = {d.name for d in cls.PARAMETERS}
    return {n: values[n] for n in (names or declared) if n in declared and n in values}


def remus_values(values=None):
    """The REMUS 100 set as the plugins' parameter names: the frozen file plus the site's latitude and the one
    kinematic viscosity."""
    v = dict(REMUS if values is None else values)
    v["diameter"] = v["beam"]
    v["latitude"] = REMUS_FILE["gravity"]["latitude_rad"]
    v["kinematic_viscosity"] = KINEMATIC_VISCOSITY
    return v


def remus(*, full=False, coriolis="kirchhoff_full", current="yaw_rate_terms", values=None, actuators=None,
          vehicle_params=None):
    """REMUS 100: the hull alone (a prescribed wrench, ``full=False``) or with the fin pairs and the propeller.
    ``values`` replaces the frozen parameter set (a perturbed set in a test)."""
    p = plugins()
    v = remus_values(values)
    vehicle_p = {**pick(p.SpheroidAuv, v), **({"coriolis_form": coriolis} if coriolis else {}),
                 "current_form": current, **(vehicle_params or {})}
    if actuators is None:
        actuators = ([(p.FinPairsDeflectionOnly, pick(p.FinPairsDeflectionOnly, v)), (p.Propeller, pick(p.Propeller, v))]
                     if full else [p.PrescribedWrench])
    return build_context(p.SpheroidAuv, vehicle_p,
                         hydrostatics=(p.SubmergedRestoring, pick(p.SubmergedRestoring, v, ("center_of_buoyancy",))),
                         hydrodynamics=(p.AuvHullLoads, pick(p.AuvHullLoads, v)), actuators=actuators)


def otter(*, current="none", hydrostatics=None, hydrodynamics=None, actuators=(), vehicle_params=None):
    """The Otter as a catamaran: the frozen MSS set, no force producer unless ``actuators`` gives some."""
    p = plugins()
    o = dict(OTTER)
    o["kinematic_viscosity"] = KINEMATIC_VISCOSITY
    vehicle_p = {**pick(p.Catamaran, o), "current_form": current, **(vehicle_params or {})}
    return build_context(
        p.Catamaran, vehicle_p,
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


def mss_to_ours(x):
    """MSS state [nu; eta] to ours [eta; nu]."""
    x = np.asarray(x)
    return np.concatenate([x[..., 6:12], x[..., 0:6]], axis=-1)


ours_to_mss = mss_to_ours
