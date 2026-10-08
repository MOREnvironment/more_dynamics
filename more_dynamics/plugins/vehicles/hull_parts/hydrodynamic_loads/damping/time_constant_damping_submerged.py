"""Time-constant linear damping of a submerged body (Dmtrx.m; remus100.m 218)
as Luka's ``HydrodynamicsModel``: the first output is the signed force the
vehicle adds (``vehicle_model.capnp`` 18-21). It consumes the rigid-body and
added-mass matrices, the weight and both centres by their payload names.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.hydrodynamic_load_parts import (
    time_constant_damping_submerged, time_constant_damping_submerged_parameters)
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class TimeConstantDampingSubmerged(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("time_constants", [20, 20, 1]),  # remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217)
        ParameterDescription("damping_ratios", [0.3, 0.8]),  # remus100.m:194, 195 [zeta4 zeta5]
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, time_constant_damping_submerged(),
                                   time_constant_damping_submerged_parameters())

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("TimeConstantDampingSubmerged must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("hydrodynamic_force", out["tau"], "signed generalized force added by the vehicle, BODY")
        for name in ("D", "damping_coefficients"):
            io.output(payload_name(name), out[name])
        return io.payload()
