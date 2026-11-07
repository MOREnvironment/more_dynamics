"""The hydrodynamic loads on a fully submerged hull whose potential damping comes from frequency tables, as Luka's
``HydrodynamicsModel``: the viscous damping that completes the reduced potential damping, its speed-dependent
factors, the hull's lift and drag and the cross-flow strip integral. The first output is the sum, the signed generalized
force the vehicle adds; the loads follow as ``damping_force``, ``lift_drag_force`` and ``cross_flow_force``, the
damping matrix as ``damping_matrix`` and the viscous part alone as ``viscous_damping``.

Inputs by name, from the vehicle (a ``SubmergedBemAuv``): ``rigid_body_mass_matrix``, ``added_mass_matrix``,
``potential_damping_matrix`` (about the CO), ``restoring_stiffness`` and ``table_reference_point`` (the point R the
tables are about), ``relative_velocity``, ``length``, ``span``, ``draft``, ``section_beam``, ``water_density``,
``kinematic_viscosity``. Frames: BODY FRD about the CO, SI units.

Equations (keys in References):

* the total damping about R, ``D_R = B_eq + Bv``: ``Bv_ii = M_ii / T_i - B_eq,ii`` for surge, sway, heave and yaw
  (``T_1236``) and ``Bv_ii = 2 delta_zeta_i sqrt(M_ii G_ii)`` for roll and pitch (``delta_zeta_45``), ``M = M_RB +
  A_eq`` about R (Fossen 2025; ``models/frequency_reduction/viscous_damping.py``); the matrices are taken from the CO
  to R by ``H(-r)^T X H(-r)`` and the damping back by ``H(r)^T D H(r)`` (Fossen 2011, eqs. 3.24, 3.27, pp. 49-50);
* its speed factors at R: ``D_11 exp(-k_u |u_r|)`` (``surge_damping_blend`` ``k_u``) and ``D_ii (1 + kappa_i
  |nu_r,i|)`` for roll, pitch and yaw (``kappa_456``); the damping force is ``-D(nu_r) nu_r``;
* lift and drag of the hull as a slender wing of span ``height``: ``S = p_factor L H``, ``C_D0 = c_d pi (H/2)^2 / S``,
  ``alpha = atan2(w_r, u_r)``, the Oswald factor ``e_oswald`` (``models/lift_drag/lift_drag.py``, the
  sources cited there);
* the cross-flow strips (Fossen 2011, eqs. 6.91-6.92, p. 127): the 2-D section drag of a circular cylinder at the
  Reynolds number of the section breadth (``cross_flow_model = "cylinder"``) or Hoerner's rectangle (``"hoerner"``).

A time constant that asks for a negative viscous damping (``T_i > M_ii / B_eq,ii``) and a hull with heave stiffness
(``G33`` not zero: a surface-piercing craft) are refused when the vehicle is assembled.

``smooth_speed`` is a numerical guard of the lift and drag, not a property of the hull (``auv_hull_loads.py``).

Fidelity: the time constants and the damping ratios are tuned values of the manoeuvring model (raise: decay tests or
an identification from logs); the lift, drag and cross-flow are physics forms with printed coefficients (raise:
towing-tank or log identification).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Eqs. 3.24, 3.27,
    pp. 49-50; eqs. 6.91-6.92, p. 127.
[Fossen 2025] Fossen, T. I. (2025). Maneuvering coefficient estimation from frequency-dependent added mass and damping:
    a power-based approach. Ocean Engineering 341:122494 (the viscous damping; not read, not on disk).

Author:    Enio Krizman
Date:      2026-11-05
"""
import casadi as ca
import numpy as np

from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription
from rpp_py.lifecycle import CallbackResult

from more_transformations.more_casadi_transformations import Parameter, freeze, unset

from more_dynamics.models.cross_flow.cross_flow_strip import (
    cross_flow_strip_circular_cylinder_reynolds, cross_flow_strip_rectangular_section_hoerner)
from more_dynamics.models.damping.hydrodynamic_load_parts import (
    HULL_PARASITIC_DRAG_PARAMETERS, HULL_PLANFORM_PARAMETERS, hull_lift_drag, hull_lift_drag_parameters,
    hull_parasitic_drag, hull_planform_area)
from more_dynamics.models.frequency_reduction.viscous_damping import (
    SUBMERGED_VISCOUS_DAMPING_PARAMETERS, move_to_table_reference, scaled_table_damping, submerged_viscous_damping)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, read_declared, refuses_unset
from more_dynamics.plugins.shared.vehicle_graph import CompositionError

CROSS_FLOW_MODELS = {"cylinder": cross_flow_strip_circular_cylinder_reynolds,
                     "hoerner": cross_flow_strip_rectangular_section_hoerner}
HEAVE_STIFFNESS_TOLERANCE = 1e-10      # G33 below this (relative to max(1, |G|_F)) is an unrestrained heave
SPEED_FACTOR_PARAMETERS = (
    Parameter("kappa_456", (3, 1), "s/rad", "quadratic factors of the roll, pitch and yaw damping", 0.0),
    Parameter("surge_damping_blend", (1, 1), "s/m", "exponent k_u of the surge damping D_11 exp(-k_u |u_r|)", 0.0),
)
LIFT_DRAG_OWN = (
    Parameter("p_factor", (1, 1), "1", "planform area factor, S = p_factor length height", 0.0, 1.0,
              minimum_exclusive=True),
    Parameter("c_d", (1, 1), "1", "drag coefficient of the hull on its frontal area", 0.0, minimum_exclusive=True),
    Parameter("e_oswald", (1, 1), "1", "Oswald efficiency factor of the induced drag", 0.0, 1.0,
              minimum_exclusive=True),
)


def _number(expression):
    try:
        return np.asarray(ca.evalf(expression), float)
    except RuntimeError:
        return None


class BemHullLoads(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("T_1236", unset()),  # s, time constants of the total damping [surge, sway, heave, yaw] (Fossen 2025)
        ParameterDescription("delta_zeta_45", unset()),  # damping-ratio increments of roll and pitch (Fossen 2025)
        ParameterDescription("kappa_456", unset()),  # s/rad, D_ii (1 + kappa_i |nu_i|) for roll, pitch, yaw
        ParameterDescription("surge_damping_blend", unset()),  # s/m, k_u of D_11 exp(-k_u |u_r|)
        ParameterDescription("p_factor", unset()),  # planform area S = p_factor length height
        ParameterDescription("c_d", unset()),  # frontal drag coefficient, C_D0 = c_d pi (height/2)^2 / S
        ParameterDescription("e_oswald", unset()),  # Oswald efficiency of the induced drag
        ParameterDescription("cross_flow_model", "cylinder"),  # 2-D section drag of the strips: "cylinder" or "hoerner"
        ParameterDescription("smooth_speed", True),  # numerical guard, finite Jacobian at rest (models/lift_drag)
        ParameterDescription("smooth_speed_epsilon", 1e-09),  # m/s, numerical regularisation of the guard
    ]

    def __init__(self) -> None:
        self._viscous = self._scaled = self._lift_drag = self._planform = self._drag = self._cross_flow = None

    @refuses_unset
    def on_configure(self, context: ComponentContext) -> CallbackResult:
        model = context.get_parameter("cross_flow_model")
        if model not in CROSS_FLOW_MODELS:
            raise ValueError(f"BemHullLoads: cross_flow_model must be one of {sorted(CROSS_FLOW_MODELS)}, got {model!r}")
        smooth = context.get_parameter("smooth_speed")
        self._viscous = freeze(submerged_viscous_damping(), SUBMERGED_VISCOUS_DAMPING_PARAMETERS,
                               read_declared(context, SUBMERGED_VISCOUS_DAMPING_PARAMETERS))
        self._scaled = freeze(scaled_table_damping(), SPEED_FACTOR_PARAMETERS,
                              read_declared(context, SPEED_FACTOR_PARAMETERS))
        own = read_declared(context, LIFT_DRAG_OWN)
        self._planform = freeze(hull_planform_area(), HULL_PLANFORM_PARAMETERS, {"planform_fraction": own["p_factor"]})
        self._drag = freeze(hull_parasitic_drag(), HULL_PARASITIC_DRAG_PARAMETERS,
                            {"hull_frontal_drag_coefficient": own["c_d"]})
        lift_drag_declared = tuple(d for d in hull_lift_drag_parameters(smooth_speed=smooth)
                                   if d.name in ("oswald_efficiency", "smooth_speed_epsilon"))
        values = {"oswald_efficiency": own["e_oswald"]}
        if smooth:
            values["smooth_speed_epsilon"] = context.get_parameter("smooth_speed_epsilon")
        self._lift_drag = freeze(hull_lift_drag(smooth_speed=smooth), lift_drag_declared, values)
        self._cross_flow = CROSS_FLOW_MODELS[model]()
        return CallbackResult()

    def check_vehicle(self, known) -> None:
        """Refuse a hull with heave stiffness, and time constants that need a negative viscous damping (the vehicle's
        matrices are constants here)."""
        g = _number(known["restoring_stiffness"])
        if g is not None and abs(g[2, 2]) > HEAVE_STIFFNESS_TOLERANCE * max(1.0, float(np.linalg.norm(g))):
            raise CompositionError("BemHullLoads: G33 is not zero, the hull is surface-piercing (a submerged hull has "
                                   "no heave stiffness)")
        r = known["table_reference_point"]
        to_r = move_to_table_reference()
        mass = to_r(known["rigid_body_mass_matrix"] + known["added_mass_matrix"], r)
        damping = to_r(known["potential_damping_matrix"], r)
        margin = _number(self._viscous(mass, damping, known["restoring_stiffness"])[1])
        if margin is not None and margin.min() < 0.0:
            m, b = _number(mass), _number(damping)
            bad = [(i + 1, m[i, i] / b[i, i] if b[i, i] > 0 else float("inf"))
                   for k, i in enumerate((0, 1, 2, 5)) if margin.ravel()[k] < 0.0]
            raise CompositionError("BemHullLoads: T_1236 needs a negative viscous damping for degree(s) of freedom "
                                   + ", ".join(f"{dof} (choose T <= M/B_eq = {limit:.4g} s)" for dof, limit in bad))

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._viscous is None:
            raise RuntimeError("BemHullLoads must be configured before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        r = io.input("table_reference_point", 3)
        mass = io.input("rigid_body_mass_matrix", (6, 6)) + io.input("added_mass_matrix", (6, 6))
        to_r = move_to_table_reference()
        mass_r = to_r(mass, r)  # M about R (Fossen 2011, eq. 3.24, p. 49)
        potential_r = to_r(io.input("potential_damping_matrix", (6, 6)), r)
        viscous = self._viscous(mass_r, potential_r, io.input("restoring_stiffness", (6, 6)))[0]  # Bv (Fossen 2025)
        nu_r = io.input("relative_velocity", 6)
        damping = self._scaled(potential_r + viscous, nu_r, r)  # D(nu_r) about the CO
        damping_force = -damping @ nu_r
        planform = io.call(self._planform)["planform_area"]  # S = p_factor L H
        drag = io.call(self._drag, given={"planform_area": planform})["parasitic_drag_coefficient"]
        lift_drag = io.call(self._lift_drag, given={"planform_area": planform,
                                                    "parasitic_drag_coefficient": drag})["tau"]
        cross_flow = io.call(self._cross_flow, rename={"beam": "section_beam"})["tau"]
        io.output("hydrodynamic_force", damping_force + lift_drag + cross_flow,
                  "signed generalized force added by the vehicle, BODY")
        io.output("damping_force", damping_force)
        io.output("lift_drag_force", lift_drag)
        io.output("cross_flow_force", cross_flow)
        io.output("damping_matrix", damping, "D(nu_r) about the CO, BODY")
        io.output("viscous_damping", viscous, "viscous part Bv about the table reference point")
        return io.payload()
