"""The hydrodynamic loads on the hull of a slender submerged vehicle, as
Luka's ``HydrodynamicsModel``: time-constant linear damping (``Dmtrx.m``;
``remus100.m`` 218-219), slender-body lift and induced drag
(``forceLiftDrag.m``) and the cross-flow strip integral with the
circular-cylinder section (``crossFlowDrag.m`` 55-70, ``cylinderDrag.m``
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
  (``remus100.m`` 134; ``length`` and ``span`` are the vehicle's, so a
  change of the hull moves ``S``), ``"given"`` is ``planform_area``. Raise:
  the planform from the hull lines.
* ``parasitic_drag_coefficient_method``: ``"computed"`` is ``C_D0 = Cd pi
  (D/2)^2 / S`` with the frontal drag coefficient ``Cd`` (``remus100.m``
  143-145; Allen et al. 2000, as MSS cites it, not read), ``"given"`` is
  ``parasitic_drag_coefficient``. Raise: a measured drag of the hull in a tow.

Forms of the loads, each a selector with its source (the next level of
fidelity needs the measurement named):

* ``body_lift_form``: ``"vertical_plane_wing"`` (default: the lift half of the
  slender-body wing block, on the vertical-plane angle of attack only, no sway
  force and no moment) or ``"point_force_both_planes"`` (the body lift in both
  planes as one point force at the pressure centre, Prestero 2001, eqs.
  4.27-4.36, pp. 30-31; the moment has the sign of the cross product, not of the
  printed eq. 4.36, ``models/lift_drag/body_lift.py``). Raise: lift measured in
  a tow or on a rotating arm.
* ``body_lift_slope_method`` (``point_force_both_planes``):
  ``"slender_body_hoerner"`` reads ``body_lift_slope_per_degree`` and refuses a
  hull outside ``6.7 <= l/d <= 10`` (Prestero 2001, eqs. 4.30-4.32, p. 30), or
  ``"given"`` reads ``body_lift_slope`` (1/rad). ``pressure_centre_method``:
  ``"fraction_of_length_from_nose"`` reads ``pressure_centre_fraction`` (the
  nose at ``+l/2``, Prestero 2001, eq. 4.35, p. 31) or ``"given"`` reads
  ``pressure_centre`` (m, BODY).
* ``hull_drag_form``: ``"wing_drag"`` (default: the drag half of the wing block,
  parasitic plus induced on the vertical-plane lift) or ``"quadratic_axial"``
  (``X = -1/2 rho c_d A_f u|u|`` with ``hull_frontal_drag_coefficient``, Prestero
  2001, eqs. 4.5-4.6, p. 25).
* ``roll_damping_form``: ``"linear"`` (default: the time-constant damping alone)
  or ``"linear_plus_quadratic"`` (adds ``K = K_p|p| p|p|`` from
  ``roll_quadratic_damping``, Prestero 2001, eqs. 4.10-4.11, p. 27, into
  ``damping_force``; ``damping_matrix`` stays the linear matrix).

Outputs ``body_lift_force`` and ``hull_drag_force`` are the two halves of
``lift_drag_force`` (their sum).

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
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc
    thesis. Eqs. 4.5-4.6, p. 25; eqs. 4.10-4.11, p. 27; eqs. 4.27-4.36, pp. 30-31.
[MSS] Fossen, T. I. MSS, MIT, @ 6e78e98: LIBRARY/modeling/Dmtrx.m;
    CRAFT/AUV/models/remus100.m 134, 144-145, 193-197, 218-222;
    LIBRARY/modeling/forceLiftDrag.m; coeffLiftDrag.m 54;
    LIBRARY/modeling/crossFlowDrag.m 55-70; cylinderDrag.m 78-110.

Author:    Enio Krizman
Date:      2026-10-09
"""
import casadi as ca

from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import check_values, unset

from more_dynamics.models.cross_flow.cross_flow_strip import cross_flow_strip_circular_cylinder_reynolds
from more_dynamics.models.damping.hydrodynamic_load_parts import (
    HULL_PARASITIC_DRAG_PARAMETERS, HULL_PLANFORM_PARAMETERS, hull_parasitic_drag, hull_planform_area,
    hull_wing_loads, hull_wing_loads_parameters, time_constant_damping_submerged,
    time_constant_damping_submerged_parameters)
from more_dynamics.models.damping.quadratic_loads import (
    QUADRATIC_AXIAL_DRAG_PARAMETERS, QUADRATIC_ROLL_DAMPING_PARAMETERS, quadratic_axial_drag, quadratic_roll_damping)
from more_dynamics.models.lift_drag.body_lift import (
    FRACTION_PRESSURE_CENTRE_PARAMETERS, GIVEN_LIFT_SLOPE_PARAMETERS, GIVEN_PRESSURE_CENTRE_PARAMETERS,
    SLENDER_BODY_LIFT_SLOPE_PARAMETERS, check_slenderness, fraction_of_length_pressure_centre,
    point_force_body_lift, slender_body_lift_slope)
from more_dynamics.plugins.shared.payload_io import (
    PayloadBuilder, frozen_block, payload_name, read_declared, read_value, refuses_unset)
from more_dynamics.plugins.shared.vehicle_graph import CompositionError

DERIVED_METHODS = ("computed", "given")
BODY_LIFT_FORMS = ("vertical_plane_wing", "point_force_both_planes")
BODY_LIFT_SLOPE_METHODS = ("slender_body_hoerner", "given")
PRESSURE_CENTRE_METHODS = ("fraction_of_length_from_nose", "given")
HULL_DRAG_FORMS = ("wing_drag", "quadratic_axial")
ROLL_DAMPING_FORMS = ("linear", "linear_plus_quadratic")


def _constant(expression):
    """The number an expression of the vehicle's geometry evaluates to, or ``None`` when it depends on a symbol."""
    try:
        return float(ca.evalf(expression))
    except RuntimeError:
        return None


class AuvHullLoads(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("time_constants", unset()),  # remus100.m:193, 194, 197 [T1 T2 T6] (Dmtrx.m call at :217)
        ParameterDescription("damping_ratios", unset()),  # remus100.m:195, 196 [zeta4 zeta5]
        ParameterDescription("oswald_efficiency", unset()),  # coeffLiftDrag.m:54 e = 0.3
        ParameterDescription("planform_area_method", "computed"),  # S = f L D (remus100.m:134)
        ParameterDescription("planform_fraction", unset()),  # remus100.m:134 S = 70 % of the rectangle L_auv * D_auv
        ParameterDescription("planform_area", unset()),  # given value, read with "given": remus100.m:134 S = 0.7 * L_auv * D_auv
        ParameterDescription("parasitic_drag_coefficient_method", "computed"),  # CD_0 = Cd pi (D/2)^2 / S (remus100.m:143-145)
        ParameterDescription("hull_frontal_drag_coefficient", unset()),  # remus100.m:143 Cd = 0.42, Allen et al. (2000) as MSS cites it
        ParameterDescription("parasitic_drag_coefficient", unset()),  # given value, read with "given": remus100.m:144-145 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry)
        ParameterDescription("body_lift_form", "vertical_plane_wing"),  # Prestero 2001, eqs. 4.27-4.36, pp. 30-31 (point force); the wing block's lift half otherwise
        ParameterDescription("body_lift_slope_method", "slender_body_hoerner"),  # Prestero 2001, eqs. 4.30-4.32, p. 30
        ParameterDescription("body_lift_slope_per_degree", unset()),  # 1/deg, Prestero 2001, eq. 4.31, p. 30 (c_yb, Hoerner's value for 6.7 <= l/d <= 10)
        ParameterDescription("body_lift_slope", unset()),  # 1/rad, given value, read with "given"
        ParameterDescription("pressure_centre_method", "fraction_of_length_from_nose"),  # Prestero 2001, eq. 4.35, p. 31
        ParameterDescription("pressure_centre_fraction", unset()),  # fraction of the length from the nose, Prestero 2001, eq. 4.35, p. 31 (0.65; Hoerner's range 0.6-0.7)
        ParameterDescription("pressure_centre", unset()),  # m, x from the CO, given value, read with "given"
        ParameterDescription("hull_drag_form", "wing_drag"),  # the wing block's drag half, or Prestero 2001, eqs. 4.5-4.6, p. 25
        ParameterDescription("roll_damping_form", "linear"),  # the time-constant damping alone, or plus Prestero 2001, eqs. 4.10-4.11, p. 27
        ParameterDescription("roll_quadratic_damping", unset()),  # N m s^2/rad^2, K_p|p| of eq. 4.11, p. 27; not positive
        ParameterDescription("smooth_speed", True),  # numerical guard, finite Jacobian at rest (models/lift_drag, models/damping)
        ParameterDescription("smooth_speed_epsilon", 1e-09),  # m/s, numerical regularisation of the guard, not a vehicle quantity
    ]

    def __init__(self) -> None:
        self._damping = self._wing = self._cross_flow = self._planform = self._parasitic_drag = None
        self._axial_drag = self._roll_damping = self._body_lift = self._lift_slope = self._pressure_centre = None
        self._given_planform = 0.0
        self._given_slope = self._given_centre = None
        self._lift_form = self._drag_form = self._slope_method = None

    @refuses_unset
    def initialize(self, context: ComponentContext) -> None:
        selectors = {}
        for name, allowed in (("planform_area_method", DERIVED_METHODS), ("parasitic_drag_coefficient_method", DERIVED_METHODS),
                              ("body_lift_form", BODY_LIFT_FORMS), ("hull_drag_form", HULL_DRAG_FORMS),
                              ("roll_damping_form", ROLL_DAMPING_FORMS)):
            selectors[name] = context.get_parameter(name)
            if selectors[name] not in allowed:
                raise ValueError(f"AuvHullLoads: {name} must be one of {allowed}, got {selectors[name]!r}")
        planform_method = selectors["planform_area_method"]
        drag_method = selectors["parasitic_drag_coefficient_method"]
        lift_form, drag_form = selectors["body_lift_form"], selectors["hull_drag_form"]
        smooth_speed = context.get_parameter("smooth_speed")
        self._lift_form, self._drag_form = lift_form, drag_form
        self._damping = frozen_block(context, time_constant_damping_submerged(smooth_speed=smooth_speed),
                                     time_constant_damping_submerged_parameters(smooth_speed=smooth_speed))
        # the wing block is built for the halves that are its own: the lift half (vertical_plane_wing) and the drag half (wing_drag)
        loads = tuple(load for load, on in (("lift", lift_form == "vertical_plane_wing"), ("drag", drag_form == "wing_drag")) if on)
        self._wing = self._planform = self._parasitic_drag = self._given_planform = None
        if loads:
            computed = {"planform_area": planform_method == "computed",
                        "parasitic_drag_coefficient": "drag" in loads and drag_method == "computed"}
            declared = tuple(d for d in hull_wing_loads_parameters(smooth_speed=smooth_speed, loads=loads)
                             if not computed.get(d.name))
            self._wing = frozen_block(context, hull_wing_loads(smooth_speed=smooth_speed, loads=loads), declared)
            self._planform = (frozen_block(context, hull_planform_area(), HULL_PLANFORM_PARAMETERS)
                              if computed["planform_area"] else None)
            self._parasitic_drag = (frozen_block(context, hull_parasitic_drag(), HULL_PARASITIC_DRAG_PARAMETERS)
                                    if computed["parasitic_drag_coefficient"] else None)
            self._given_planform = (None if computed["planform_area"]
                                    else float(read_value(context, "planform_area", "m^2", "given planform area of the hull")))
        self._axial_drag = (frozen_block(context, quadratic_axial_drag(), QUADRATIC_AXIAL_DRAG_PARAMETERS)
                            if drag_form == "quadratic_axial" else None)
        self._roll_damping = (frozen_block(context, quadratic_roll_damping(), QUADRATIC_ROLL_DAMPING_PARAMETERS)
                              if selectors["roll_damping_form"] == "linear_plus_quadratic" else None)
        self._body_lift = self._lift_slope = self._pressure_centre = self._given_slope = self._given_centre = None
        self._slope_method = None
        if lift_form == "point_force_both_planes":
            self._point_force_body_lift(context)
        self._cross_flow = cross_flow_strip_circular_cylinder_reynolds()

    def _point_force_body_lift(self, context):
        """The point-force body lift: its slope and pressure centre each from its method (computed or given)."""
        slope_method = context.get_parameter("body_lift_slope_method")
        centre_method = context.get_parameter("pressure_centre_method")
        for name, method, allowed in (("body_lift_slope_method", slope_method, BODY_LIFT_SLOPE_METHODS),
                                      ("pressure_centre_method", centre_method, PRESSURE_CENTRE_METHODS)):
            if method not in allowed:
                raise ValueError(f"AuvHullLoads: {name} must be one of {allowed}, got {method!r}")
        self._slope_method = slope_method
        if slope_method == "slender_body_hoerner":
            self._lift_slope = frozen_block(context, slender_body_lift_slope(), SLENDER_BODY_LIFT_SLOPE_PARAMETERS)
        else:
            self._given_slope = float(check_values(GIVEN_LIFT_SLOPE_PARAMETERS, read_declared(
                context, GIVEN_LIFT_SLOPE_PARAMETERS))["body_lift_slope"])
        if centre_method == "fraction_of_length_from_nose":
            self._pressure_centre = frozen_block(context, fraction_of_length_pressure_centre(),
                                                 FRACTION_PRESSURE_CENTRE_PARAMETERS)
        else:
            self._given_centre = float(check_values(GIVEN_PRESSURE_CENTRE_PARAMETERS, read_declared(
                context, GIVEN_PRESSURE_CENTRE_PARAMETERS))["pressure_centre"])
        self._body_lift = point_force_body_lift()

    def check_vehicle(self, known) -> None:
        """Called by the vehicle with the quantities it gives: the slender-body slope is stated for
        ``6.7 <= l/d <= 10`` only (Prestero 2001, eq. 4.31, p. 30), so a vehicle outside that range is refused here,
        naming the ratio."""
        if self._slope_method != "slender_body_hoerner":
            return
        length, span = _constant(known["length"]), _constant(known["span"])
        if length is None or span is None:
            return
        try:
            check_slenderness(length, span)
        except ValueError as error:
            raise CompositionError(f"AuvHullLoads: {error}") from None

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._damping is None:
            raise RuntimeError("AuvHullLoads must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        rename = {"beam": "section_beam"}
        damping = io.call(self._damping)
        damping_force = damping["tau"]
        if self._roll_damping is not None:
            damping_force = damping_force + io.call(self._roll_damping)["tau"]  # K_p|p| p|p| (Prestero 2001, eq. 4.10, p. 27)
        derived = {}  # the wing's quantities the hull's geometry defines, computed from the vehicle's length and span
        wing = None
        if self._wing is not None:
            planform = ca.SX(self._given_planform) if self._planform is None else io.call(self._planform)["planform_area"]
            if self._planform is not None:
                derived["planform_area"] = planform
            if self._parasitic_drag is not None:
                derived["parasitic_drag_coefficient"] = io.call(
                    self._parasitic_drag, given={"planform_area": planform})["parasitic_drag_coefficient"]
            wing = io.call(self._wing, given=derived)
        if self._lift_form == "vertical_plane_wing":
            body_lift = wing["lift_tau"]
        else:
            slope = (ca.SX(self._given_slope) if self._lift_slope is None
                     else io.call(self._lift_slope)["body_lift_slope"])
            centre = (ca.SX(self._given_centre) if self._pressure_centre is None
                      else io.call(self._pressure_centre)["pressure_centre"])
            body_lift = io.call(self._body_lift, given={"body_lift_slope": slope, "pressure_centre": centre})["tau"]
        hull_drag = wing["drag_tau"] if self._drag_form == "wing_drag" else io.call(self._axial_drag)["tau"]
        # both halves from the wing block: its own sum, the block of lift_drag_casadi to the last bit
        lift_drag = wing["tau"] if (self._lift_form == "vertical_plane_wing" and self._drag_form == "wing_drag") \
            else body_lift + hull_drag
        cross_flow = io.call(self._cross_flow, rename=rename)
        io.output("hydrodynamic_force", damping_force + lift_drag + cross_flow["tau"],
                  "signed generalized force added by the vehicle, BODY")
        io.output("damping_force", damping_force)
        io.output("lift_drag_force", lift_drag)
        io.output("body_lift_force", body_lift, "lift half of lift_drag_force, BODY, N and N m about the CO")
        io.output("hull_drag_force", hull_drag, "drag half of lift_drag_force, BODY, N and N m about the CO")
        io.output("cross_flow_force", cross_flow["tau"])
        for name in ("D", "damping_coefficients"):
            io.output(payload_name(name), damping[name])
        return io.payload()
