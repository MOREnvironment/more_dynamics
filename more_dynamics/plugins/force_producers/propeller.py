"""Propeller with a linearised open-water characteristic and given
coefficients (remus100.m 149-177): thrust ``(1 - t) X`` along the shaft axis
plus the shaft torque about it, the shaft speed (rpm) saturated at
``max_shaft_speed`` (remus100.m 116). Luka's ``ForceProducer``: the first
output is the wrench on the vehicle; relative velocity and water density come
from the vehicle by name, the shaft speed is the command.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ ed2ce7b: CRAFT/AUV/models/remus100.m 116,
    149-177, 250-253 (propeller); the equations and their lines are cited in
    ``more_dynamics.models.propeller.propeller``.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import unset

from more_dynamics.models.propeller.propeller_linearized_open_water import (
    propeller_linearized_open_water, propeller_linearized_open_water_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, refuses_unset


class Propeller(ForceProducer):
    PARAMETERS = [
        ParameterDescription("propeller_diameter", unset()),  # remus100.m:149 D_prop
        ParameterDescription("max_shaft_speed", unset()),  # remus100.m:111 n_max
        ParameterDescription("thrust_deduction", unset()),  # remus100.m:150 t_prop
        ParameterDescription("wake_fraction", unset()),  # remus100.m:151 Va = 0.944 * U_r, w = 1 - 0.944
        ParameterDescription("pitch_diameter_ratio", unset()),  # remus100.m:157 wageningen(0,1,0.718,3), second argument
        ParameterDescription("blade_area_ratio", unset()),  # remus100.m:156-157 blade-area ratio 0.718
        ParameterDescription("blade_count", unset()),  # remus100.m:156-157 3 blades
        ParameterDescription("max_advance_number", unset()),  # remus100.m:154 Ja_max
        ParameterDescription("roll_moment_scale", unset()),  # remus100.m:253 tau(4) = K_prop / 10
        ParameterDescription("position", unset()),  # remus100.m:250-253 thrust on x_b through the CO (no moment arm)
        ParameterDescription("orientation", unset()),  # remus100.m:250, 253 shaft along x_b (thrust in tau(1), torque in tau(4))
        ParameterDescription("thrust_torque_coefficients", unset()),  # remus100.m:158, 159, 161, 162 [KT_0 KQ_0 KT_max KQ_max]
    ]

    def __init__(self) -> None:
        self._model = None

    @refuses_unset
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
