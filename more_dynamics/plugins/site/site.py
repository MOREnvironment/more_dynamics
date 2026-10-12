"""The site a vehicle floats in: the water density, the acceleration of
gravity and the kinematic viscosity, as a ``SiteModel`` (``site_model.capnp``).
A stateless graph with no input; the vehicle connects every input of its
other children to the output of the same name (``water_density``,
``gravity``, ``kinematic_viscosity``). One quantity, one value: a vehicle
holds no copy of these, it takes them from the component in its ``site``
slot, and several vehicles in the same water link the one component.

``site_form`` sets where the gravity comes from, each with the measurement
that raises it to the next level:

* ``"latitude"``: WGS-84 normal gravity at the geodetic ``latitude``
  (``gravity.m`` 11-12), through ``more_transformations``. Raise: the local
  gravity measured at the site.
* ``"given"``: ``gravity`` as a value, 9.81 m/s^2 to match MSS from release
  2.0.4 (``LIBRARY/mssConstants.m`` 14 @ 49e03e3; ``otter.m`` 90-91).

``latitude`` is read only with ``"latitude"`` and ``gravity`` only with
``"given"``; the water density and the kinematic viscosity are read in both.
Every vehicle number is unset in this plugin and given by the composition (a named part of the
workspace); each reference part's values and their provenance are in ``DEFAULTS.md``. A value the
selected forms read and the composition leaves unset refuses the build, naming it with its unit and
meaning.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. WGS-84 normal gravity, through
    more_transformations' ``ECEFNEDtransform.gravity`` (cited there).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: INS/functions/gravity.m 11-12;
    CRAFT/AUV/models/remus100.m 96-98; CRAFT/USV/models/otter.m 90-91; and,
    @ 49e03e3 (release 2.0.5), LIBRARY/mssConstants.m 12-14.
[HydrostaticsModel] Mandic, L. more_dynamics ``LinearSurfaceHydrostatics``
    (linear_surface_hydrostatics.py): the form of a stateless plugin whose
    ``graph()`` carries named outputs.

Author:    Enio Krizman
Date:      2026-10-10
"""
from rpp_plugin_types.more_dynamics import SiteModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import unset

from more_dynamics.models.site.site import (
    SITE_AT_LATITUDE_PARAMETERS, SITE_GIVEN_PARAMETERS, site_at_latitude, site_given)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, refuses_unset
from more_dynamics.plugins.shared.vehicle_graph import CompositionError

SITE_FORMS = {"latitude": (site_at_latitude, SITE_AT_LATITUDE_PARAMETERS),
              "given": (site_given, SITE_GIVEN_PARAMETERS)}
OUTPUT_DESCRIPTIONS = {"water_density": "water density, kg/m^3", "gravity": "acceleration of gravity, m/s^2",
                       "kinematic_viscosity": "kinematic viscosity of the water, m^2/s"}


class Site(SiteModel):
    PARAMETERS = [
        ParameterDescription("site_form", "given"),  # otter.m 90-91; gravity.m 11-12 for "latitude"
        ParameterDescription("latitude", unset()),  # rad, remus100.m 96-97 mu = deg2rad(63.446827); read with site_form "latitude"
        ParameterDescription("gravity", unset()),  # m/s^2, mssConstants.m:14 @ 49e03e3 (MSS 2.0.5), otter.m 90-91; read with site_form "given"
        ParameterDescription("water_density", unset()),  # kg/m^3, otter.m 90-91
        ParameterDescription("kinematic_viscosity", unset()),  # m^2/s, cylinderDrag.m 78-80 nu_water = 1e-6; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
    ]

    def __init__(self) -> None:
        self._model = None

    @refuses_unset
    def initialize(self, context: ComponentContext) -> None:
        form = context.get_parameter("site_form")
        if form not in SITE_FORMS:
            raise CompositionError(f"Site: site_form must be one of {sorted(SITE_FORMS)}, got {form!r}")
        function, declared = SITE_FORMS[form]
        self._model = frozen_block(context, function(), declared)

    def graph(self) -> SiteModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("Site must be initialized before graph()")
        io = PayloadBuilder(SiteModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(name, out[name], OUTPUT_DESCRIPTIONS[name])
        return io.payload()
