"""Uniform NED current; body acceleration with the yaw-rate terms only (remus100.m 118-122, an MSS form).

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import CurrentModel
from rpp_py.context import ComponentContext

from more_dynamics.models.vehicles.hull_parts.current.current import horizontal_current_yaw_rate_terms, HORIZONTAL_CURRENT_PARAMETERS
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name
from rpp_py.parameter_description import ParameterDescription


class HorizontalCurrentYawRateTerms(CurrentModel):
    PARAMETERS = [
        ParameterDescription("current_speed", 0.0),  # m/s, horizontal speed; still water unless the composition sets one
        ParameterDescription("current_direction", 0.0),  # rad, current set (NED, from north, clockwise)
        ParameterDescription("current_vertical_speed", 0.0),  # m/s, positive down (NED)
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, horizontal_current_yaw_rate_terms(), HORIZONTAL_CURRENT_PARAMETERS)

    def graph(self) -> CurrentModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("HorizontalCurrentYawRateTerms must be initialized before graph()")
        io = PayloadBuilder(CurrentModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
