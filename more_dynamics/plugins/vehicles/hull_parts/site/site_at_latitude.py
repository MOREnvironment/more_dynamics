"""Site with WGS-84 normal gravity at a latitude, one water density and one viscosity for every part.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import SiteModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.site.site import site_at_latitude, SITE_AT_LATITUDE_PARAMETERS
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class SiteAtLatitude(SiteModel):
    PARAMETERS = [
        ParameterDescription("latitude", 1.1073560310932362),  # remus100.m 96-97, mu = deg2rad(63.446827)
        ParameterDescription("water_density", 1026),  # remus100.m:98 rho = 1026, the one value (MSS: 1026 at imlay61.m:31 and forceLiftDrag.m:26, 1025 at crossFlowDrag.m:36)
        ParameterDescription("kinematic_viscosity", 1e-06),  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, site_at_latitude(), SITE_AT_LATITUDE_PARAMETERS)

    def graph(self) -> SiteModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("SiteAtLatitude must be initialized before graph()")
        io = PayloadBuilder(SiteModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
