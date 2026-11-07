"""Restoring of a submerged body, as Luka's ``HydrostaticsModel``: weight and
buoyancy are separate, with a centre of gravity and a centre of buoyancy of
their own. The first output is the signed generalized force the vehicle adds
(``vehicle_model.capnp`` 11-13), so ``restoring_force = -g`` (Fossen 2011,
eq. 4.5, p. 60, through the submerged block); weight, buoyancy and the centre
of buoyancy follow for the models that consume them.

``buoyancy_method`` sets how the buoyancy is known, each with the
measurement that raises it to the next level:

* ``"neutral"``: ``B = W`` (Fossen 2011, eq. 4.7, p. 61; MSS ``remus100.m``
  215), a neutrally buoyant body. Raise: weigh the vehicle in water at the
  quay (reserve buoyancy).
* ``"from_volume"``: ``B = rho g nabla`` with ``displaced_volume`` given
  (Fossen 2011, eq. 4.1, p. 59). Raise: take the volume from the measured
  hull, then ``"given"`` from a weighing.
* ``"given"``: ``B`` as a value, e.g. a measured reserve of buoyancy. A
  positively buoyant vehicle (``B > W``) floats up and trims: the buoyancy
  and the weight act at different points, ``center_of_buoyancy`` against the
  vehicle's ``center_of_gravity``.

* ``"immersed_body_of_revolution"``: the hull is a body of revolution given by
  a station table (``hull_profile_x``, ``hull_profile_radius``,
  ``hull_axis_origin``, ``hull_profile_slices``); ``B = rho g V(eta)`` and its
  centre follow the part of the hull under the calm free surface, section by
  section with exact circle segments, so the body floats; fully submerged it is
  ``"from_volume"`` exactly (``models/restoring/immersed.py``; Fossen 2011,
  eqs. 4.1, 4.5, pp. 59-60). The radii are scaled so the hull displaces
  ``displaced_volume``; the centroid of the table must lie within
  ``hull_profile_cb_tolerance`` (m) of ``center_of_buoyancy``, or the build is
  refused. A neutral buoyancy cannot float, so this method has no neutral
  form. Raise: the station table from the hull lines; the volume from a
  weighing in water.

``displaced_volume`` and ``buoyancy`` are read only by their own method.
Every vehicle number is unset in this plugin and given by the composition (a named part of the
workspace); each reference part's values and their provenance are in ``DEFAULTS.md``. A value the
selected forms read and the composition leaves unset refuses the build, naming it with its unit and
meaning.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 4.1, 4.5, 4.7, pp. 59-61.
[MSS] Fossen, T. I. MSS, MIT, @ 6e78e98: CRAFT/AUV/models/remus100.m 3, 96-97,
    132-133, 139, 215.

Author:    Enio Krizman
Date:      2026-10-09
"""
from rpp_plugin_types.more_dynamics import HydrostaticsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription
from rpp_py.lifecycle import CallbackResult

import casadi as ca

from more_transformations.more_casadi_transformations import unset

from more_dynamics.models.restoring.immersed import immersed_restoring, immersed_restoring_parameters

from more_dynamics.models.restoring.restoring_parts import submerged_restoring, submerged_restoring_parameters
from more_dynamics.plugins.shared.payload_io import (
    PayloadBuilder, frozen_block, payload_name, read_value, refuses_unset)

IMMERSED = "immersed_body_of_revolution"


class SubmergedRestoring(HydrostaticsModel):
    PARAMETERS = [
        ParameterDescription("buoyancy_method", "neutral"),  # remus100.m:215 B = W
        ParameterDescription("center_of_buoyancy", unset()),  # remus100.m:139 r_bB
        ParameterDescription("displaced_volume", unset()),  # m^3, 4/3 pi (L/2)(D/2)^2, remus100.m:132-133
        ParameterDescription("buoyancy", unset()),  # N, m g with m = 31.9 kg (remus100.m:3), g at remus100.m:96-97
        ParameterDescription("hull_profile_x", unset()),  # m, stations along BODY x from the axis point, read with "immersed_body_of_revolution"
        ParameterDescription("hull_profile_radius", unset()),  # m, section radius at each station, read with "immersed_body_of_revolution"
        ParameterDescription("hull_axis_origin", unset()),  # m, axis point at station 0, BODY FRD from the CO, read with "immersed_body_of_revolution"
        ParameterDescription("hull_profile_slices", 40),  # sections per interval between two stations (numerical resolution)
        ParameterDescription("hull_profile_cb_tolerance", unset()),  # m, largest accepted distance of the table's centroid from center_of_buoyancy
    ]

    def __init__(self) -> None:
        self._model = None

    @refuses_unset
    def on_configure(self, context: ComponentContext) -> CallbackResult:
        method = context.get_parameter("buoyancy_method")
        if method != IMMERSED:
            self._model = frozen_block(context, submerged_restoring(method), submerged_restoring_parameters(method))
            return CallbackResult()
        stations = read_value(context, "hull_profile_x", "m", "stations of the hull profile along BODY x",
                              placeholder=[0.0, 1.0])
        slices = int(context.get_parameter("hull_profile_slices"))
        self._model = frozen_block(context, immersed_restoring(stations, slices),
                                   immersed_restoring_parameters(len(stations)))
        tolerance = float(read_value(context, "hull_profile_cb_tolerance", "m",
                                     "largest distance of the profile's centroid from center_of_buoyancy"))
        probe = {n: ca.DM.zeros(*self._model.size_in(n)) for n in self._model.name_in()}
        probe["mass"], probe["gravity"], probe["water_density"] = 1.0, 1.0, 1.0
        offset = float(self._model(**probe)["profile_centroid_offset"])  # independent of the pose and the mass
        if not offset <= tolerance:
            raise ValueError(f"SubmergedRestoring: the centroid of the hull profile is {offset:.3g} m from "
                             f"center_of_buoyancy (hull_profile_cb_tolerance {tolerance:g} m): one centre of buoyancy")
        return CallbackResult()

    def graph(self) -> HydrostaticsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("SubmergedRestoring must be configured before graph()")
        io = PayloadBuilder(HydrostaticsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("restoring_force", -out["g"], "signed generalized force added by the vehicle, BODY, N and N m")
        for name in self._model.name_out():
            if name not in ("g", "profile_centroid_offset"):
                io.output(payload_name(name), out[name])
        return io.payload()
