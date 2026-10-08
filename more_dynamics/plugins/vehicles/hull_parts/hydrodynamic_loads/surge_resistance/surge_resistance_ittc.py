"""ITTC-1957 surge resistance with a form factor (Fossen 2011, eqs.
6.82-6.85, p. 125; XuuITTC.m 32-39) and the linear-to-quadratic blend, as
Luka's ``HydrodynamicsModel``: the first output is the signed force the
vehicle adds. It takes the mass, wetted surface, length, water density and
kinematic viscosity from the vehicle by name; a restoring part that gives no
``wetted_surface`` is refused by the vehicle.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 6.82-6.85, p. 125.
[MSS] Fossen, T. I. MSS, MIT: LIBRARY/modeling/XuuITTC.m 32-39 @ cc07579;
    forceSurgeDamping.m 57 @ ac77394; addedMassSurge.m 33 @ cc07579;
    CRAFT/SHIP/models/osv.m 128 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.hydrodynamic_load_parts import (surge_resistance_ittc,
                                                                         surge_resistance_ittc_parameters)
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block


class SurgeResistanceIttc(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("time_constant", 100.0),  # osv.m 128, vessel.T1 = 100 s
        ParameterDescription("form_factor", 0.1),  # XuuITTC.m 33, k = 0.1; Fossen 2011, p. 125
        ParameterDescription("crossover_speed", 2.0),  # forceSurgeDamping.m 57 @ ac77394, u_cross = 2 m/s
        ParameterDescription("surge_added_mass_factor", 2.7),  # addedMassSurge.m 33, A11 = 2.7 rho nabla^(5/3) / L^2
        ParameterDescription("ittc_reynolds_floor", 100000.0),  # XuuITTC.m 35, Re_min = 1e5
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, surge_resistance_ittc(), surge_resistance_ittc_parameters())

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("SurgeResistanceIttc must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("hydrodynamic_force", out["tau"], "signed generalized force added by the vehicle, BODY")
        return io.payload()
