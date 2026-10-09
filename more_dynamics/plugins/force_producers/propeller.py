"""Propeller with a linearised open-water characteristic and given
coefficients (remus100.m 148-176): thrust ``(1 - t) X`` along the shaft axis
plus the shaft torque about it, the shaft speed (rpm) saturated at
``max_shaft_speed`` (remus100.m 115). Luka's ``ForceProducer``: the first
output is the wrench on the vehicle; relative velocity and water density come
from the vehicle by name, the shaft speed is the command.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 115,
    148-176, 249-252 (propeller); the equations and their lines are cited in
    ``more_dynamics.models.propeller.propeller``.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.propeller.propeller_linearized_open_water import (
    propeller_linearized_open_water, propeller_linearized_open_water_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block


class Propeller(ForceProducer):
    PARAMETERS = [
        ParameterDescription("propeller_diameter", 0.14),  # remus100.m:148 D_prop
        ParameterDescription("max_shaft_speed", 1525),  # remus100.m:110 n_max
        ParameterDescription("thrust_deduction", 0.1),  # remus100.m:149 t_prop
        ParameterDescription("wake_fraction", 0.05600000000000005),  # remus100.m:150 Va = 0.944 * U_r, w = 1 - 0.944
        ParameterDescription("pitch_diameter_ratio", 1),  # remus100.m:156 wageningen(0,1,0.718,3), second argument
        ParameterDescription("blade_area_ratio", 0.718),  # remus100.m:155-156 blade-area ratio 0.718
        ParameterDescription("blade_count", 3),  # remus100.m:155-156 3 blades
        ParameterDescription("max_advance_number", 0.6632),  # remus100.m:153 Ja_max
        ParameterDescription("roll_moment_scale", 0.1),  # remus100.m:252 tau(4) = K_prop / 10
        ParameterDescription("position", [0, 0, 0]),  # remus100.m:249-252 thrust on x_b through the CO (no moment arm)
        ParameterDescription("orientation", [0, 0, 0]),  # remus100.m:249, 252 shaft along x_b (thrust in tau(1), torque in tau(4))
        ParameterDescription("thrust_torque_coefficients", [0.4566, 0.07, 0.1798, 0.0312]),  # remus100.m:157, 158, 160, 161 [KT_0 KQ_0 KT_max KQ_max]
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, propeller_linearized_open_water(),
                                   propeller_linearized_open_water_parameters())

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("Propeller must be initialized before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        out = io.call(self._model, rename={"command": "shaft_speed_command"})
        io.output("generated_force", out["tau"], "propeller wrench about the CO, BODY, N and N m")
        for name in ("shaft_axis", "coefficients_in_use"):
            if name in self._model.name_out():
                io.output(name, out[name])
        return io.payload()
