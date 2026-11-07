"""Two deflection-only fin pairs (rudder and stern plane), the form of
remus100.m 229-246: each pair's lift ``0.5 rho U^2 A C_L delta`` at its
position, the deflections saturated at ``max_deflection`` (remus100.m
114-115). Luka's ``ForceProducer``: the first output is the wrench on the
vehicle; relative velocity and water density come from the vehicle by name,
the two deflections (rudder, stern plane) are the command.

``convention`` declares the sign convention of the model, ``"starboard_down_positive"``: a positive rudder deflection
gives a force along -y (a turn to starboard), a positive stern-plane deflection a force along -z (bow down). It is the
only value the model has; another convention is refused (a ``Fin`` takes its sign from its ``lift_axis``).

``fin_positions_method`` sets how the two fin positions are known:
``"computed"`` (default) puts both pairs at the tail of the body, ``x = -a``
with ``a`` the body's half length (``semi_major_axis``, read from the vehicle
by name; ``remus100.m`` 185, 190), so a change of the hull moves them;
``"given"`` takes ``rudder_position`` and ``stern_plane_position`` as values.
Raise: the measured positions of the fin axes.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ 6e78e98: CRAFT/AUV/models/remus100.m
    114-115, 180-190, 229-246 (fin pairs); the equations and their lines
    are cited in ``more_dynamics.models.fin.fins``.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription
from rpp_py.lifecycle import CallbackResult

from more_transformations.more_casadi_transformations import unset

from more_dynamics.models.fin.fin_pairs_deflection_only import (FIN_PAIRS_CONVENTION, fin_pair_positions,
                                                                    fin_pairs_deflection_only,
                                                                    fin_pairs_deflection_only_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, refuses_unset


FIN_POSITIONS_METHODS = ("computed", "given")


class FinPairsDeflectionOnly(ForceProducer):
    PARAMETERS = [
        ParameterDescription("rudder_area", unset()),  # remus100.m:184 A_r = 2 * S_fin (S_fin :179)
        ParameterDescription("stern_plane_area", unset()),  # remus100.m:189 A_s = 2 * S_fin (S_fin :179)
        ParameterDescription("rudder_lift_coefficient", unset()),  # remus100.m:183 CL_delta_r
        ParameterDescription("stern_plane_lift_coefficient", unset()),  # remus100.m:188 CL_delta_s
        ParameterDescription("convention", FIN_PAIRS_CONVENTION),  # the model's own: refused when another
        ParameterDescription("fin_positions_method", "computed"),  # x_r = x_s = -a (remus100.m:185, 190)
        ParameterDescription("rudder_position", unset()),  # given value, read with "given": remus100.m:185 x_r = -a, a = L_auv / 2 (one geometry)
        ParameterDescription("stern_plane_position", unset()),  # given value, read with "given": remus100.m:190 x_s = -a, a = L_auv / 2 (one geometry)
        ParameterDescription("max_deflection", unset()),  # remus100.m:110 delta_max = deg2rad(20)
    ]

    def __init__(self) -> None:
        self._model = self._positions = None

    @refuses_unset
    def on_configure(self, context: ComponentContext) -> CallbackResult:
        convention = context.get_parameter("convention")
        if convention != FIN_PAIRS_CONVENTION:
            raise ValueError(f"FinPairsDeflectionOnly: its model has the convention {FIN_PAIRS_CONVENTION!r}, got {convention!r}")
        method = context.get_parameter("fin_positions_method")
        if method not in FIN_POSITIONS_METHODS:
            raise ValueError(f"FinPairsDeflectionOnly: fin_positions_method must be one of {FIN_POSITIONS_METHODS}, "
                             f"got {method!r}")
        given = method == "given"
        self._model = frozen_block(context, fin_pairs_deflection_only(),
                                   fin_pairs_deflection_only_parameters(positions_given=given))
        self._positions = None if given else fin_pair_positions()
        return CallbackResult()

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinPairsDeflectionOnly must be configured before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        positions = {} if self._positions is None else dict(io.call(self._positions).items())
        out = io.call(self._model, rename={"command": "fin_deflection_command"}, given=positions)
        io.output("generated_force", out["tau"], "fin wrench about the CO, BODY, N and N m")
        return io.payload()
