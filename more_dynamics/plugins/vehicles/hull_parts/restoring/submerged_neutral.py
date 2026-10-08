"""Restoring of a submerged body held neutral, W = m g = B (remus100.m 214),
as Luka's ``HydrostaticsModel``: the first output is the signed generalized
force the vehicle adds (``vehicle_model.capnp`` 11-13), so ``restoring_force
= -g`` (Fossen 2011, eq. 4.6, p. 60, through the submerged block); weight,
buoyancy and the centre of buoyancy follow for the parts that consume them.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HydrostaticsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.restoring.restoring_parts import SUBMERGED_NEUTRAL_PARAMETERS, submerged_neutral
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class SubmergedNeutral(HydrostaticsModel):
    PARAMETERS = [
        ParameterDescription("center_of_buoyancy", [0, 0, 0]),  # remus100.m:138 r_bB
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, submerged_neutral(), SUBMERGED_NEUTRAL_PARAMETERS)

    def graph(self) -> HydrostaticsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("SubmergedNeutral must be initialized before graph()")
        io = PayloadBuilder(HydrostaticsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("restoring_force", -out["g"], "signed generalized force added by the vehicle, BODY, N and N m")
        for name in self._model.name_out():
            if name != "g":
                io.output(payload_name(name), out[name])
        return io.payload()
