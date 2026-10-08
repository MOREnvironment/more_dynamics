"""Slender-body hull lift and induced drag (forceLiftDrag.m) as Luka's
``HydrodynamicsModel``; span from the hull form, water density from the site.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.hydrodynamic_load_parts import hull_lift_drag, hull_lift_drag_parameters
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block


class HullLiftDrag(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("planform_area", 0.2128),  # remus100.m:133 S = 0.7 * L_auv * D_auv
        ParameterDescription("parasitic_drag_coefficient", 0.05595961914206819),  # remus100.m:143-144 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry)
        ParameterDescription("oswald_efficiency", 0.3),  # coeffLiftDrag.m:54 e = 0.3
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, hull_lift_drag(), hull_lift_drag_parameters())

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("HullLiftDrag must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("hydrodynamic_force", out["tau"], "signed generalized force added by the vehicle, BODY")
        return io.payload()
