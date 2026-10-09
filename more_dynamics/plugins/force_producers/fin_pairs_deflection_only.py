"""Two deflection-only fin pairs (rudder and stern plane), the form of
remus100.m 228-245: each pair's lift ``0.5 rho U^2 A C_L delta`` at its
position, the deflections saturated at ``max_deflection`` (remus100.m
113-114). Luka's ``ForceProducer``: the first output is the wrench on the
vehicle; relative velocity and water density come from the vehicle by name,
the two deflections (rudder, stern plane) are the command.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m
    113-114, 179-189, 228-245 (fin pairs); the equations and their lines
    are cited in ``more_dynamics.models.fin.fins``.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.fin.fin_pairs_deflection_only import (fin_pairs_deflection_only,
                                                                    fin_pairs_deflection_only_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block


class FinPairsDeflectionOnly(ForceProducer):
    PARAMETERS = [
        ParameterDescription("rudder_area", 0.0133),  # remus100.m:183 A_r = 2 * S_fin (S_fin :179)
        ParameterDescription("stern_plane_area", 0.0133),  # remus100.m:188 A_s = 2 * S_fin (S_fin :179)
        ParameterDescription("rudder_lift_coefficient", 0.5),  # remus100.m:182 CL_delta_r
        ParameterDescription("stern_plane_lift_coefficient", 0.7),  # remus100.m:187 CL_delta_s
        ParameterDescription("rudder_position", -0.8),  # remus100.m:184 x_r = -a, a = L_auv / 2 (one geometry)
        ParameterDescription("stern_plane_position", -0.8),  # remus100.m:189 x_s = -a, a = L_auv / 2 (one geometry)
        ParameterDescription("max_deflection", 0.3490658503988659),  # remus100.m:109 delta_max = deg2rad(20)
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, fin_pairs_deflection_only(), fin_pairs_deflection_only_parameters())

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinPairsDeflectionOnly must be initialized before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        out = io.call(self._model, rename={"command": "fin_deflection_command"})
        io.output("generated_force", out["tau"], "fin wrench about the CO, BODY, N and N m")
        return io.payload()
