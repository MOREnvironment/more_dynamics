"""Two deflection-only fin pairs (rudder and stern plane), the form of
remus100.m 228-245: each pair's lift ``0.5 rho U^2 A C_L delta`` at its
position, the deflections saturated at ``max_deflection`` (remus100.m
113-114). Luka's ``ForceProducer``: the first output is the wrench on the
vehicle; relative velocity and water density come from the vehicle by name,
the two deflections (rudder, stern plane) are the command.

``fin_positions_method`` sets how the two fin positions are known:
``"computed"`` (default) puts both pairs at the tail of the body, ``x = -a``
with ``a`` the body's half length (``semi_major_axis``, read from the vehicle
by name; ``remus100.m`` 184, 189), so a change of the hull moves them;
``"given"`` takes ``rudder_position`` and ``stern_plane_position`` as values.
Raise: the measured positions of the fin axes.

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

from more_dynamics.models.fin.fin_pairs_deflection_only import (fin_pair_positions, fin_pairs_deflection_only,
                                                                    fin_pairs_deflection_only_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block


FIN_POSITIONS_METHODS = ("computed", "given")


class FinPairsDeflectionOnly(ForceProducer):
    PARAMETERS = [
        ParameterDescription("rudder_area", 0.0133),  # remus100.m:183 A_r = 2 * S_fin (S_fin :179)
        ParameterDescription("stern_plane_area", 0.0133),  # remus100.m:188 A_s = 2 * S_fin (S_fin :179)
        ParameterDescription("rudder_lift_coefficient", 0.5),  # remus100.m:182 CL_delta_r
        ParameterDescription("stern_plane_lift_coefficient", 0.7),  # remus100.m:187 CL_delta_s
        ParameterDescription("fin_positions_method", "computed"),  # x_r = x_s = -a (remus100.m:184, 189)
        ParameterDescription("rudder_position", -0.8),  # given value, read with "given": remus100.m:184 x_r = -a, a = L_auv / 2 (one geometry)
        ParameterDescription("stern_plane_position", -0.8),  # given value, read with "given": remus100.m:189 x_s = -a, a = L_auv / 2 (one geometry)
        ParameterDescription("max_deflection", 0.3490658503988659),  # remus100.m:109 delta_max = deg2rad(20)
    ]

    def __init__(self) -> None:
        self._model = self._positions = None

    def initialize(self, context: ComponentContext) -> None:
        method = context.get_parameter("fin_positions_method")
        if method not in FIN_POSITIONS_METHODS:
            raise ValueError(f"FinPairsDeflectionOnly: fin_positions_method must be one of {FIN_POSITIONS_METHODS}, "
                             f"got {method!r}")
        given = method == "given"
        self._model = frozen_block(context, fin_pairs_deflection_only(),
                                   fin_pairs_deflection_only_parameters(positions_given=given))
        self._positions = None if given else fin_pair_positions()

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinPairsDeflectionOnly must be initialized before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        positions = {} if self._positions is None else dict(io.call(self._positions).items())
        out = io.call(self._model, rename={"command": "fin_deflection_command"}, given=positions)
        io.output("generated_force", out["tau"], "fin wrench about the CO, BODY, N and N m")
        return io.payload()
