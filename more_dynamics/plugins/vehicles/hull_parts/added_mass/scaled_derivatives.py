"""Added mass scaled from the hull's mass and inertia (otter.m 152-159).

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import AddedMassModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.added_mass.added_mass_parts import scaled_derivatives, SCALED_DERIVATIVES_PARAMETERS
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class ScaledDerivatives(AddedMassModel):
    PARAMETERS = [
        ParameterDescription("added_mass_coefficients", [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]),  # otter.m 152-157
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, scaled_derivatives(), SCALED_DERIVATIVES_PARAMETERS)

    def graph(self) -> AddedMassModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("ScaledDerivatives must be initialized before graph()")
        io = PayloadBuilder(AddedMassModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
