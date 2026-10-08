"""Added mass of the ideal prolate spheroid, Lamb's derivatives (imlay61.m 31-59).

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import AddedMassModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.added_mass.added_mass_parts import lamb_spheroid, lamb_spheroid_parameters
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class LambSpheroid(AddedMassModel):
    PARAMETERS = [
        ParameterDescription("roll_added_inertia_ratio", 0.3),  # remus100.m:136 r44
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, lamb_spheroid(), lamb_spheroid_parameters())

    def graph(self) -> AddedMassModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("LambSpheroid must be initialized before graph()")
        io = PayloadBuilder(AddedMassModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
