"""Hull hydrodynamic loads of a body whose damping, lift and the rigid
body's own Coriolis-centripetal terms are given together as one set of
published nondimensional derivatives (``models/coefficient_loads``) plus a
cross-flow strip integral with a combined, speed-normalised drag term
(``models/cross_flow``), as Luka's ``HydrodynamicsModel``: the first output
is the signed generalized force the vehicle adds (``vehicle_model.capnp``
18-21).

Inputs by name, from the vehicle: ``relative_velocity``, ``mass``,
``inertia_diagonal``, ``inertia_products``, ``center_of_gravity``,
``length``, ``water_density``. Own parameters: the 52 hull-only
nondimensional derivatives and the cross-flow geometry (the 14 acceleration
derivatives are read by the vehicle's added-mass block, the remaining 27
published ones by no plugin yet: ``models/coefficient_loads`` module
docstring).

**No rigid-body Coriolis matrix of the vehicle's own** applies alongside
this block: a vehicle using this hull-load form folds the rigid-body
Coriolis-centripetal terms into the hull force directly
(``models/coefficient_loads`` module docstring; ``DerivativeAuv.coriolis_matrix``
is zero).

References
----------
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/npsauv.m 119-122, 218-236,
    242-275 @ 6e78e98.

Author:    Enio Krizman
Date:      2026-10-09
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import unset

from more_dynamics.models.coefficient_loads.coefficient_loads import (
    coefficient_hull_loads, coefficient_hull_loads_parameters)
from more_dynamics.models.cross_flow.cross_flow_strip import (
    COMBINED_DRAG_CROSS_FLOW_PARAMETERS, cross_flow_strip_combined_drag)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, read_value, refuses_unset


class CoefficientHullLoads(HydrodynamicsModel):
    PARAMETERS = [
        # 52 hull-only derivatives (npsauv.m 127-157): no generic physics default exists for an empirical
        # coefficient table, a composition gives every one of them, each with its kind (measured, predicted,
        # assumed, zero with its reason); declared here once, by the part that reads them.
        *[ParameterDescription(d.name, unset()) for d in coefficient_hull_loads_parameters()],
        ParameterDescription("cross_flow_drag_coefficients", unset()),  # npsauv.m:221 [Cdy, Cdz]
        ParameterDescription("cross_flow_section_dimensions", unset()),  # npsauv.m:221-222 [Hx, Bx]
        ParameterDescription("cross_flow_sections", unset()),  # npsauv.m:219, number of strips (vehicle data)
    ]

    def __init__(self) -> None:
        self._hull = self._cross_flow = None

    @refuses_unset
    def initialize(self, context: ComponentContext) -> None:
        sections = int(read_value(context, "cross_flow_sections", "1", "number of strips of the cross-flow integral",
                                 placeholder=1))
        if sections < 1:
            raise ValueError(f"CoefficientHullLoads: cross_flow_sections must be >= 1, got {sections}")
        hull_declared = coefficient_hull_loads_parameters()
        self._hull = frozen_block(context, coefficient_hull_loads(), hull_declared)
        cross_flow_declared = tuple(d for d in COMBINED_DRAG_CROSS_FLOW_PARAMETERS if d.name != "cross_flow_sections")
        self._cross_flow = frozen_block(context, cross_flow_strip_combined_drag(sections=sections),
                                        cross_flow_declared)

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._hull is None:
            raise RuntimeError("CoefficientHullLoads must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        rename = {"nu_r": "relative_velocity"}
        hull = io.call(self._hull)
        cross_flow = io.call(self._cross_flow, rename=rename)
        io.output("hydrodynamic_force", hull["hydrodynamic_force"] + cross_flow["tau"],
                  "signed generalized force added by the vehicle, BODY")
        io.output("coefficient_force", hull["hydrodynamic_force"])
        io.output("cross_flow_force", cross_flow["tau"])
        return io.payload()
