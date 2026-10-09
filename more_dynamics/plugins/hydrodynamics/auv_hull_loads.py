"""The hydrodynamic loads on the hull of a slender submerged vehicle, as
Luka's ``HydrodynamicsModel``: time-constant linear damping (``Dmtrx.m``;
``remus100.m`` 217-218), slender-body lift and induced drag
(``forceLiftDrag.m``) and the cross-flow strip integral with the
circular-cylinder section (``crossFlowDrag.m`` 54-69, ``cylinderDrag.m``
78-110). The first output is the sum, the signed generalized force the
vehicle adds (``vehicle_model.capnp`` 18-21); the three loads follow as
``damping_force``, ``lift_drag_force`` and ``cross_flow_force``, and the
damping matrix as ``damping_matrix``.

Inputs by name, from the vehicle and the hydrostatics: ``relative_velocity``,
the two mass matrices, ``weight``, ``center_of_gravity``,
``center_of_buoyancy``, ``span``, ``water_density``, ``length``, ``draft``,
``section_beam``, ``kinematic_viscosity``. Own parameters: the damping time
constants and ratios, the Oswald factor of the lift, and the two
quantities of the lift that the hull's geometry defines.

Derived quantities, each with a ``*_method`` (computed by its equation, the
default, or given as a value):

* ``planform_area_method``: ``"computed"`` is ``S = f L D`` with the
  ``planform_fraction`` ``f`` of the rectangle length x diameter
  (``remus100.m`` 133; ``length`` and ``span`` are the vehicle's, so a
  change of the hull moves ``S``), ``"given"`` is ``planform_area``. Raise:
  the planform from the hull lines.
* ``parasitic_drag_coefficient_method``: ``"computed"`` is ``C_D0 = Cd pi
  (D/2)^2 / S`` with the frontal drag coefficient ``Cd`` (``remus100.m``
  142-144; Allen et al. 2000, as MSS cites it, not read), ``"given"`` is
  ``parasitic_drag_coefficient``. Raise: a measured drag of the hull in a tow.

``smooth_speed`` is a numerical guard, not a property of the vehicle: with
it ``True`` (default) the angle of attack of the lift and the speed of the
surge-damping fade are regularised by ``smooth_speed_epsilon`` (m/s) so the
Jacobian of the vehicle is finite at rest (``models/lift_drag``,
``models/damping``); the loads equal the MSS form to about
``eps^2 / speed^2`` of their value away from rest. ``False`` is the MSS form,
whose Jacobian is not finite at zero relative velocity (a stiff integrator
cannot start from rest).

Fidelity: the time-constant damping is an MSS calibration (raise: damping
identified from decay tests); the lift and the cross-flow are physics forms
with printed coefficients (raise: towing-tank or log identification).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 6.91-6.92, p. 127 (strip integral).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/Dmtrx.m;
    CRAFT/AUV/models/remus100.m 133, 143-144, 192-196, 217-221;
    LIBRARY/modeling/forceLiftDrag.m; coeffLiftDrag.m 54;
    LIBRARY/modeling/crossFlowDrag.m 54-69; cylinderDrag.m 78-110.

Author:    Enio Krizman
Date:      2026-10-09
"""
import casadi as ca

from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.cross_flow.cross_flow_strip import cross_flow_strip_circular_cylinder_reynolds
from more_dynamics.models.damping.hydrodynamic_load_parts import (
    HULL_PARASITIC_DRAG_PARAMETERS, HULL_PLANFORM_PARAMETERS, hull_lift_drag, hull_lift_drag_parameters,
    hull_parasitic_drag, hull_planform_area, time_constant_damping_submerged,
    time_constant_damping_submerged_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, payload_name

DERIVED_METHODS = ("computed", "given")


class AuvHullLoads(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("time_constants", [20, 20, 1]),  # remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217)
        ParameterDescription("damping_ratios", [0.3, 0.8]),  # remus100.m:194, 195 [zeta4 zeta5]
        ParameterDescription("oswald_efficiency", 0.3),  # coeffLiftDrag.m:54 e = 0.3
        ParameterDescription("planform_area_method", "computed"),  # S = f L D (remus100.m:133)
        ParameterDescription("planform_fraction", 0.7),  # remus100.m:133 S = 70 % of the rectangle L_auv * D_auv
        ParameterDescription("planform_area", 0.2128),  # given value, read with "given": remus100.m:133 S = 0.7 * L_auv * D_auv
        ParameterDescription("parasitic_drag_coefficient_method", "computed"),  # CD_0 = Cd pi (D/2)^2 / S (remus100.m:142-144)
        ParameterDescription("hull_frontal_drag_coefficient", 0.42),  # remus100.m:142 Cd = 0.42, Allen et al. (2000) as MSS cites it
        ParameterDescription("parasitic_drag_coefficient", 0.05595961914206819),  # given value, read with "given": remus100.m:143-144 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry)
        ParameterDescription("smooth_speed", True),  # numerical guard, finite Jacobian at rest (models/lift_drag, models/damping)
        ParameterDescription("smooth_speed_epsilon", 1e-09),  # m/s, numerical regularisation of the guard, not a vehicle quantity
    ]

    def __init__(self) -> None:
        self._damping = self._lift_drag = self._cross_flow = self._planform = self._parasitic_drag = None
        self._given_planform = 0.0

    def initialize(self, context: ComponentContext) -> None:
        planform_method = context.get_parameter("planform_area_method")
        drag_method = context.get_parameter("parasitic_drag_coefficient_method")
        for name, method in (("planform_area_method", planform_method),
                             ("parasitic_drag_coefficient_method", drag_method)):
            if method not in DERIVED_METHODS:
                raise ValueError(f"AuvHullLoads: {name} must be one of {DERIVED_METHODS}, got {method!r}")
        smooth_speed = context.get_parameter("smooth_speed")
        self._damping = frozen_block(context, time_constant_damping_submerged(smooth_speed=smooth_speed),
                                     time_constant_damping_submerged_parameters(smooth_speed=smooth_speed))
        computed = {"planform_area": planform_method == "computed",
                    "parasitic_drag_coefficient": drag_method == "computed"}
        declared = tuple(d for d in hull_lift_drag_parameters(smooth_speed=smooth_speed) if not computed.get(d.name))
        self._lift_drag = frozen_block(context, hull_lift_drag(smooth_speed=smooth_speed), declared)
        self._planform = (frozen_block(context, hull_planform_area(), HULL_PLANFORM_PARAMETERS)
                          if computed["planform_area"] else None)
        self._parasitic_drag = (frozen_block(context, hull_parasitic_drag(), HULL_PARASITIC_DRAG_PARAMETERS)
                                if computed["parasitic_drag_coefficient"] else None)
        self._given_planform = None if computed["planform_area"] else float(context.get_parameter("planform_area"))
        self._cross_flow = cross_flow_strip_circular_cylinder_reynolds()

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._damping is None:
            raise RuntimeError("AuvHullLoads must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        rename = {"beam": "section_beam"}
        damping = io.call(self._damping)
        derived = {}  # the lift's quantities the hull's geometry defines, computed from the vehicle's length and span
        planform = ca.SX(self._given_planform) if self._planform is None else io.call(self._planform)["planform_area"]
        if self._planform is not None:
            derived["planform_area"] = planform
        if self._parasitic_drag is not None:
            derived["parasitic_drag_coefficient"] = io.call(
                self._parasitic_drag, given={"planform_area": planform})["parasitic_drag_coefficient"]
        lift_drag = io.call(self._lift_drag, given=derived)
        cross_flow = io.call(self._cross_flow, rename=rename)
        io.output("hydrodynamic_force", damping["tau"] + lift_drag["tau"] + cross_flow["tau"],
                  "signed generalized force added by the vehicle, BODY")
        io.output("damping_force", damping["tau"])
        io.output("lift_drag_force", lift_drag["tau"])
        io.output("cross_flow_force", cross_flow["tau"])
        for name in ("D", "damping_coefficients"):
            io.output(payload_name(name), damping[name])
        return io.payload()
