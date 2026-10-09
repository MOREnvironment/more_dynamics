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
nondimensional derivatives and the cross-flow geometry (the remaining 41 of
the published 93 derivatives are read by the added-mass block and the
actuators: ``models/coefficient_loads`` module docstring).

**No rigid-body Coriolis matrix of the vehicle's own** applies alongside
this block: the NPS AUV II folds it into the hull force directly
(``models/coefficient_loads`` module docstring; ``DerivativeAuv.coriolis_matrix``
is zero).

References
----------
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/npsauv.m 118-121, 217-235,
    241-274 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-09
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.coefficient_loads.coefficient_loads import (
    coefficient_hull_loads, coefficient_hull_loads_parameters)
from more_dynamics.models.cross_flow.cross_flow_strip import (
    COMBINED_DRAG_CROSS_FLOW_PARAMETERS, cross_flow_strip_combined_drag)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block


class CoefficientHullLoads(HydrodynamicsModel):
    PARAMETERS = [
        # 52 hull-only derivatives (npsauv.m 126-156), each 0.0 by default: no generic physics default exists for
        # an empirical coefficient table (contrast AuvHullLoads's cited REMUS defaults); a composition gives every
        # one of them (the named "NPS AUV II" part does, from nps_auv_ii_parameters.json).
        *[ParameterDescription(d.name, 0.0) for d in coefficient_hull_loads_parameters()],
        ParameterDescription("cross_flow_drag_coefficients", [0.5, 0.6]),  # npsauv.m:220 [Cdy, Cdz]
        ParameterDescription("cross_flow_section_dimensions", [0.53, 0.53]),  # npsauv.m:220-221 [Hx, Bx]
        ParameterDescription("cross_flow_sections", 10),  # npsauv.m:218, number of strips
    ]

    def __init__(self) -> None:
        self._hull = self._cross_flow = None

    def initialize(self, context: ComponentContext) -> None:
        sections = int(context.get_parameter("cross_flow_sections"))
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
