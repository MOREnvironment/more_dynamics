"""Added-mass Coriolis matrix with every term kept (m2c.m 33-48), the physics form.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import AddedMassCoriolisModel
from rpp_py.context import ComponentContext

from more_dynamics.models.vehicles.hull_parts.added_mass_coriolis.added_mass_coriolis_parts import kirchhoff_full, KIRCHHOFF_FULL_PARAMETERS
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class KirchhoffFull(AddedMassCoriolisModel):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, kirchhoff_full(), KIRCHHOFF_FULL_PARAMETERS)

    def graph(self) -> AddedMassCoriolisModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("KirchhoffFull must be initialized before graph()")
        io = PayloadBuilder(AddedMassCoriolisModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
