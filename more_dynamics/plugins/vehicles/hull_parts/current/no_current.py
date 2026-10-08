"""No current: zero current velocity and acceleration.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import CurrentModel
from rpp_py.context import ComponentContext

from more_dynamics.models.vehicles.hull_parts.current.current import no_current, NO_CURRENT_PARAMETERS
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class NoCurrent(CurrentModel):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, no_current(), NO_CURRENT_PARAMETERS)

    def graph(self) -> CurrentModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("NoCurrent must be initialized before graph()")
        io = PayloadBuilder(CurrentModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
