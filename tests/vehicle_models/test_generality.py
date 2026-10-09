"""Rule 20 (AGENTS.md, E-107): every vehicle type builds and runs with a second parameter
set that shares none of the first's numbers -- proof that the plugin's own graph is driven
by its declared parameters, not by a hardcoded vehicle number buried in the code.

"Set A" is each vehicle type's own cited ``PARAMETERS`` default (``DEFAULTS.md``, a named
part's numbers). "Set B" is a synthetic, uniformly scaled variant: no second identified
vehicle exists for every one of these four forms (``DerivativeAuv`` has only the one
published nondimensional-derivative table on disk), so this file scales set A by a fixed
factor (and nudges an exact zero off zero) instead of guessing a second vehicle's numbers
(rule 21) -- labelled synthetic throughout, never presented as a vehicle. ``hull_count``
(a 1-or-2 choice that decides which inputs the hull form even has, not a scalable physical
quantity) is pinned at the vehicle type's own structural choice in both sets instead of
scaled. Each test builds the vehicle hull alone (``PrescribedWrench``, Luka's/Fossen's
stage-1 force producer,
no parameters of its own) and checks both sets build, step to a finite state, and disagree
-- selectors (site/current/coriolis/added-mass forms, the hydrodynamics forms) are left at
the plugin's own default in both sets: this test is about numbers, not forms.

Author:    Enio Krizman
Date:      2026-10-09
"""
import numpy as np

import vehicle_contract as vc

SCALE = 0.8  # arbitrary, < 1 (several declared fractions are bounded above by 1), applied
             # to every numeric leaf of set A to build set B
ZERO_BUMP = 0.01  # an exact 0.0 also needs to differ between set A and set B
REST = np.zeros(12)


def _scaled(value):
    """Set B from set A: every float/int leaf times ``SCALE`` (an exact zero nudged off
    zero too); strings (selectors), bools and nested structure otherwise untouched."""
    if isinstance(value, bool) or isinstance(value, str) or value is None:
        return value
    if isinstance(value, (int, float)):
        scaled = value * SCALE
        return scaled if scaled != 0 else ZERO_BUMP
    if isinstance(value, list):
        return [_scaled(v) for v in value]
    return value


def _defaults(cls):
    return {d.name: d.default_value for d in cls.PARAMETERS}


def _build(vehicle_cls, vehicle_params, hydrostatics_cls, hydrostatics_params,
          hydrodynamics_cls, hydrodynamics_params):
    p = vc.plugins()
    return vc.build_context(
        vehicle_cls, {**vehicle_params, "current_form": "none"},
        hydrostatics=(hydrostatics_cls, hydrostatics_params),
        hydrodynamics=(hydrodynamics_cls, hydrodynamics_params), actuators=[p.PrescribedWrench])


def _step(vessel):
    """A finite state after one step from rest with a small prescribed wrench (``desired_wrench``,
    ``PrescribedWrench``'s own input, Luka's ``ForceProducer`` contract)."""
    graph = vc.graph_of(vessel)
    u = np.full(vc.input_size(vessel), 10.0)
    x1 = np.asarray(graph.step(REST, u)).ravel()
    assert np.all(np.isfinite(x1)), "set must step to a finite state"
    return x1


def _assert_generic(vehicle_cls, hydrostatics_cls, hydrodynamics_cls, *, hydrostatics_fixed=None):
    """``hydrostatics_fixed`` names the hydrostatics plugin's own parameters that this
    *vehicle type* (not the shared hydrostatics plugin) pins structurally -- e.g.
    ``Monohull`` always builds a one-hull form, so ``hull_count`` stays 1 in both sets
    even though ``SurfaceRestoring`` (shared with ``Catamaran``) declares 2 as its own
    class default; scaling or flipping it would ask the hull form for an input
    (``hull_lateral_offset``) the vehicle's own hull form never produces."""
    a_vehicle, a_hydrostatics, a_hydrodynamics = (
        vc.pick(vehicle_cls, _defaults(vehicle_cls)), vc.pick(hydrostatics_cls, _defaults(hydrostatics_cls)),
        vc.pick(hydrodynamics_cls, _defaults(hydrodynamics_cls)))
    a_hydrostatics.update(hydrostatics_fixed or {})
    b_vehicle = {k: _scaled(v) for k, v in _defaults(vehicle_cls).items() if k in a_vehicle}
    b_hydrostatics = {k: _scaled(v) for k, v in _defaults(hydrostatics_cls).items() if k in a_hydrostatics}
    b_hydrodynamics = {k: _scaled(v) for k, v in _defaults(hydrodynamics_cls).items() if k in a_hydrodynamics}
    b_hydrostatics.update(hydrostatics_fixed or {})
    assert a_vehicle != b_vehicle and b_vehicle, "set B must differ from set A (vehicle)"

    vessel_a = _build(vehicle_cls, a_vehicle, hydrostatics_cls, a_hydrostatics, hydrodynamics_cls, a_hydrodynamics)
    vessel_b = _build(vehicle_cls, b_vehicle, hydrostatics_cls, b_hydrostatics, hydrodynamics_cls, b_hydrodynamics)
    assert type(vessel_a).__name__ == type(vessel_b).__name__ == vehicle_cls.__name__

    x_a, x_b = _step(vessel_a), _step(vessel_b)
    assert np.abs(x_a - x_b).max() > 1e-6, "set A and set B must give different trajectories"


def test_spheroid_auv_builds_and_runs_with_a_second_synthetic_parameter_set():
    p = vc.plugins()
    _assert_generic(p.SpheroidAuv, p.SubmergedRestoring, p.AuvHullLoads)


def test_catamaran_builds_and_runs_with_a_second_synthetic_parameter_set():
    p = vc.plugins()
    _assert_generic(p.Catamaran, p.SurfaceRestoring, p.SurfaceHullLoads, hydrostatics_fixed={"hull_count": 2})


def test_monohull_builds_and_runs_with_a_second_synthetic_parameter_set():
    p = vc.plugins()
    _assert_generic(p.Monohull, p.SurfaceRestoring, p.SurfaceHullLoads, hydrostatics_fixed={"hull_count": 1})


def test_derivative_auv_builds_and_runs_with_a_second_synthetic_parameter_set():
    p = vc.plugins()
    _assert_generic(p.DerivativeAuv, p.SubmergedRestoring, p.CoefficientHullLoads)
