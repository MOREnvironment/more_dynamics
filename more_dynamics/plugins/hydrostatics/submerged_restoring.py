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

``displaced_volume`` and ``buoyancy`` are read only by their own method.
Every vehicle number is unset in this plugin and given by the composition (a named part of the
workspace); each reference part's values and their provenance are in ``DEFAULTS.md``. A value the
selected forms read and the composition leaves unset refuses the build, naming it with its unit and
meaning.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 4.1, 4.5, 4.7, pp. 59-61.
[MSS] Fossen, T. I. MSS, MIT, @ ed2ce7b: CRAFT/AUV/models/remus100.m 3, 96-97,
    132-133, 139, 215.

Author:    Enio Krizman
Date:      2026-10-09
"""
from rpp_plugin_types.more_dynamics import HydrostaticsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import unset

from more_dynamics.models.restoring.restoring_parts import submerged_restoring, submerged_restoring_parameters
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, payload_name, refuses_unset


class SubmergedRestoring(HydrostaticsModel):
    PARAMETERS = [
        ParameterDescription("buoyancy_method", "neutral"),  # remus100.m:215 B = W
        ParameterDescription("center_of_buoyancy", unset()),  # remus100.m:139 r_bB
        ParameterDescription("displaced_volume", unset()),  # m^3, 4/3 pi (L/2)(D/2)^2, remus100.m:132-133
        ParameterDescription("buoyancy", unset()),  # N, m g with m = 31.9 kg (remus100.m:3), g at remus100.m:96-97
    ]

    def __init__(self) -> None:
        self._model = None

    @refuses_unset
    def initialize(self, context: ComponentContext) -> None:
        method = context.get_parameter("buoyancy_method")
        self._model = frozen_block(context, submerged_restoring(method), submerged_restoring_parameters(method))

    def graph(self) -> HydrostaticsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("SubmergedRestoring must be initialized before graph()")
        io = PayloadBuilder(HydrostaticsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("restoring_force", -out["g"], "signed generalized force added by the vehicle, BODY, N and N m")
        for name in self._model.name_out():
            if name != "g":
                io.output(payload_name(name), out[name])
        return io.payload()
