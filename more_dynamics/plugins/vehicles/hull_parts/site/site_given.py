"""Site with a given gravity, one water density and one viscosity for every part.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import SiteModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.site.site import site_given, SITE_GIVEN_PARAMETERS
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class SiteGiven(SiteModel):
    PARAMETERS = [
        ParameterDescription("gravity", 9.81),  # otter.m 90-91
        ParameterDescription("water_density", 1025.0),  # otter.m 90-91
        ParameterDescription("kinematic_viscosity", 1e-06),  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, site_given(), SITE_GIVEN_PARAMETERS)

    def graph(self) -> SiteModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("SiteGiven must be initialized before graph()")
        io = PayloadBuilder(SiteModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
