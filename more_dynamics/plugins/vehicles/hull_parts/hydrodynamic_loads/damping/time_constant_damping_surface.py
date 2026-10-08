"""Time-constant linear damping of a surface craft, an MSS calibration from
top speed and thrust (otter.m 195-240), as Luka's ``HydrodynamicsModel``. Its
block calls the restoring matrix at the centre of flotation
``restoring_matrix``; on the payload it is ``restoring_matrix_at_flotation``,
the name its producer gives it.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.hydrodynamic_load_parts import (
    time_constant_damping_surface, time_constant_damping_surface_parameters)
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class TimeConstantDampingSurface(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("max_forward_thrust", 239.364),  # otter.m, 24.4 kgf = 24.4 * 9.81 N
        ParameterDescription("max_speed", 3.0864),  # otter.m, 6 knots = 6 * 0.5144 m/s
        ParameterDescription("time_constants", [1.0, 1.0]),  # otter.m 99-101, 202-240
        ParameterDescription("damping_ratios", [0.3, 0.2, 0.4]),  # otter.m 99-101, 202-240
        ParameterDescription("yaw_damping_nonlinearity", 10.0),  # otter.m 99-101, 202-240
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, time_constant_damping_surface(),
                                   time_constant_damping_surface_parameters())

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("TimeConstantDampingSurface must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        out = io.call(self._model, rename={"restoring_matrix": "restoring_matrix_at_flotation"})
        io.output("hydrodynamic_force", out["tau"], "signed generalized force added by the vehicle, BODY")
        for name in ("D", "damping_derivatives"):
            io.output(payload_name(name), out[name])
        return io.payload()
