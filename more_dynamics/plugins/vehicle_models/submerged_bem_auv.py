"""A fully submerged hull whose added mass and potential damping are given as frequency tables (boundary-element data
of the hull at one speed and one depth), as Luka's ``VehicleModel3D`` in the form of ``HullVessel``
(``hull_vessel.py``): the tables are reduced to constant matrices once, at build, by the power-based method
(``models/frequency_reduction``), the rigid body is the full inertia tensor about the CO, and the slots ``site``,
``hydrostatics``, ``hydrodynamics`` and ``actuators`` come from the composition. The graph is the equation of motion of
a craft in a current over the children's payloads (``plugins/shared/vehicle_graph.py``)::

    nu_r    = nu - nu_c
    nu_dot  = nu_c_dot + M^-1 ( tau - C nu_r )      M = M_RB + M_A,  C = C_RB(nu) + C_A(nu_r)
    eta_dot = J(eta) nu

``M_A = sym(H(r)^T A_eq H(r))``: the reduced added mass ``A_eq`` about the table reference point R
(``database_reference``, from the CO) moved to the CO (Fossen 2011, eqs. 3.24, 3.27, pp. 49-50). The reduced
potential damping is moved the same way and handed on, with the restoring stiffness and R, to the hull loads (a
``BemHullLoads`` completes the damping); the vehicle adds no damping itself. ``C_RB(nu)`` is parametrised by the
angular velocity alone (Fossen 2011, eq. 3.57, p. 55; Property 8.1, p. 222), so the rigid body sees an irrotational
current only through ``nu_r``.

State ``[eta; nu]`` after the actuators' states: NED position, ZYX Euler angles, BODY FRD velocity about the CO, SI
units (13 values ``[p; q; nu]`` with ``kinematics_form = "unit_quaternion"``). Inputs: the current (when
``current_form`` is not ``"none"``), then the commands of the actuators in list order.

Primitives: the tables (``frequencies``, ``added_mass``, ``potential_damping``, ``hydrostatic_stiffness_345``,
``database_reference``), the body (``body_mass``, ``inertia_diagonal`` and ``inertia_products`` about the CO,
``body_center_of_gravity``), the section geometry (``length``, ``beam``, ``height``) and the spectrum that weights
the tables (``omega_p``, ``spectrum_gamma``: a modelling choice of the reduction, not a property of the hull).

Outputs by name for the children: ``added_mass_matrix``, ``rigid_body_mass_matrix``, ``potential_damping_matrix``
(all about the CO), ``restoring_stiffness`` (about R), ``table_reference_point`` (R), ``mass``,
``center_of_gravity``, ``length``, ``draft`` (= ``height``), ``span`` (= ``height``), ``section_beam`` (= ``beam``).

Options, each with its fidelity and what raises it:

* ``coriolis_form``: ``"full"`` (default: every term of ``C_A``, the physics form of a body in an ideal
  fluid, Fossen 2011, eq. 6.43, p. 120, with the Munk moment of eq. 6.52, p. 121) or ``"couplings_removed"``
  (the pitch-heave, pitch-surge, yaw-surge and yaw-sway entries and their transposes removed: a bare hull without
  the fin and body-lift derivatives that balance the Munk moment in reality). Raise: the straight-line stability of
  the hull with its fins and body lift, checked against a logged run.
* ``current_form``: ``"full_attitude"`` (default; Fossen 2011, eqs. 8.138-8.141, 8.157, pp. 221-225), ``"none"``,
  or the shortcut forms of ``plugins/shared/vehicle_options.py``.
* ``kinematics_form``: ``"euler"`` (default) or ``"unit_quaternion"`` (Fossen 2011, eq. 2.69, p. 30).

The reduction refuses a mass matrix ``M_RB + M_A`` that is not positive definite. Every vehicle number is unset in
this plugin and given by the composition; a value left unset refuses the build, naming it.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Eqs. 3.24,
    3.27, pp. 49-50; eq. 3.44, p. 52; eq. 3.57, p. 55; eqs. 6.43, 6.52, pp. 120-121; eqs. 8.138-8.141, 8.157 and
    Property 8.1, pp. 221-225.
[Fossen 2025] Fossen, T. I. (2025). Maneuvering coefficient estimation from frequency-dependent added mass and damping:
    a power-based approach. Ocean Engineering 341:122494 (the reduction; not read, not on disk).
[HullVessel] Mandic, L. more_dynamics ``HullVessel`` (hull_vessel.py): the form of the vehicle plugin this follows.

Author:    Enio Krizman
Date:      2026-11-05
"""
from typing import List

import casadi as ca
import numpy as np

from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import unset

from more_dynamics.models.coriolis.added_mass_coriolis_parts import full, couplings_removed
from more_dynamics.models.frequency_reduction.power_based import (
    power_based_reduction, power_based_reduction_parameters)
from more_dynamics.models.frequency_reduction.viscous_damping import (
    RESTORING_STIFFNESS_PARAMETERS, TABLE_REFERENCE_PARAMETERS, move_to_body_origin, restoring_stiffness_345,
    table_reference_point)
from more_dynamics.models.hull_form.hull_form import (
    SUBMERGED_SECTION_MAIN_DIMENSIONS_PARAMETERS, submerged_section_main_dimensions)
from more_dynamics.models.rigid_body.rigid_body_parts import (
    FULL_TENSOR_RIGID_BODY_PARAMETERS, full_tensor_rigid_body, full_tensor_rigid_body_coriolis)
from more_dynamics.plugins.shared.payload_io import read_value, refuses_unset
from more_dynamics.plugins.shared.vehicle_graph import (
    INTEGRATION_STEP, CompositionError, Craft, Parts, VehicleAssembly)
from more_dynamics.plugins.shared.vehicle_options import CURRENT_FORMS, current_inputs, kinematics_form

CORIOLIS_FORMS = {"full": full, "couplings_removed": couplings_removed}


def _number(expression):
    try:
        return np.asarray(ca.evalf(expression), float)
    except RuntimeError:
        return None


class SubmergedBemAuv(VehicleModel3D):

    COMPONENTS = {
        "actuators": "List[more_dynamics::ForceProducer]",
        "site": "more_dynamics::SiteModel",
        "hydrostatics": "more_dynamics::HydrostaticsModel",
        "hydrodynamics": "more_dynamics::HydrodynamicsModel",
    }

    PARAMETERS = [
        INTEGRATION_STEP,
        ParameterDescription("frequencies", unset()),  # rad/s, increasing; an entry >= 10 rad/s stands for infinite frequency
        ParameterDescription("added_mass", unset()),  # A(omega), 36 values (6 x 6 row-major) per frequency, about R
        ParameterDescription("potential_damping", unset()),  # B(omega), 36 values per frequency, about R
        ParameterDescription("hydrostatic_stiffness_345", unset()),  # heave-roll-pitch block of G about R, 3 x 3 row-major
        ParameterDescription("database_reference", unset()),  # m, R from the CO, BODY FRD
        ParameterDescription("body_mass", unset()),  # kg
        ParameterDescription("inertia_diagonal", unset()),  # kg*m^2, [Ix, Iy, Iz] about the CO (Fossen 2011, eq. 3.44, p. 52)
        ParameterDescription("inertia_products", unset()),  # kg*m^2, [Ixy, Iyz, Ixz] about the CO
        ParameterDescription("body_center_of_gravity", unset()),  # m, CO -> CG, BODY FRD
        ParameterDescription("length", unset()),  # m, hull length (the cross-flow strips span it)
        ParameterDescription("beam", unset()),  # m, breadth of the hull section
        ParameterDescription("height", unset()),  # m, height of the hull section
        ParameterDescription("omega_p", unset()),  # rad/s, peak frequency of the spectrum that weights the tables (Fossen 2025)
        ParameterDescription("spectrum_gamma", unset()),  # JONSWAP peak factor (DNVGL-RP-C205 2017, 3.5.5.3, p. 65)
        ParameterDescription("coriolis_form", "full"),  # Fossen 2011, eq. 6.43, p. 120 (every term kept)
        ParameterDescription("kinematics_form", "euler"),  # Fossen 2011, eq. 2.40, p. 26, or eq. 2.69, p. 30
        ParameterDescription("current_form", "full_attitude"),  # Fossen 2011, eqs. 8.138-8.141, 8.157, pp. 221-225
    ]

    def __init__(self):
        self._assembly = None

    @refuses_unset
    def initialize(self, context: ComponentContext):
        form = context.get_parameter("current_form")
        frequencies = [float(v) for v in read_value(context, "frequencies", "rad/s", "frequencies of the tables",
                                                    placeholder=[1.0, 2.0])]
        reduction = power_based_reduction(frequencies)
        reduction_parameters = power_based_reduction_parameters(len(frequencies))
        own_inputs = current_inputs(context, "SubmergedBemAuv", CompositionError)

        def compute(pose, velocity, current, site):
            known = {"pose": pose, "velocity": velocity, **current, **site}
            parts = Parts(context, "SubmergedBemAuv", known)
            parts.run("hull form", submerged_section_main_dimensions(), SUBMERGED_SECTION_MAIN_DIMENSIONS_PARAMETERS)
            parts.run("table reference point", table_reference_point(), TABLE_REFERENCE_PARAMETERS)
            parts.run("restoring stiffness", restoring_stiffness_345(), RESTORING_STIFFNESS_PARAMETERS)
            parts.run("power-based reduction", reduction, reduction_parameters)  # (Fossen 2025)
            parts.run("added mass at the CO", move_to_body_origin("added_mass_at_table_reference", "M_A",
                                                                  symmetric=True))  # (Fossen 2011, eq. 3.27, p. 50)
            parts.run("potential damping at the CO", move_to_body_origin("potential_damping_at_table_reference",
                                                                         "potential_damping_matrix"))
            parts.run("rigid body", full_tensor_rigid_body(), FULL_TENSOR_RIGID_BODY_PARAMETERS)  # (eq. 3.44, p. 52)
            parts.run("rigid-body Coriolis", full_tensor_rigid_body_coriolis("book"))  # (Fossen 2011, eq. 3.57, p. 55)
            parts.run("current", CURRENT_FORMS[form]())
            current_velocity, current_acceleration = known["current_velocity"], known["current_acceleration"]
            known["relative_velocity"] = velocity - current_velocity  # nu_r = nu - nu_c (Fossen 2011, eq. 8.157)
            parts.run("added-mass Coriolis", parts.select("coriolis_form", CORIOLIS_FORMS)())
            mass_matrix = known["rigid_body_mass_matrix"] + known["added_mass_matrix"]
            numbers = _number(mass_matrix)
            if numbers is not None and np.linalg.eigvalsh(0.5 * (numbers + numbers.T)).min() <= 0.0:
                raise CompositionError("SubmergedBemAuv: the reduced mass matrix M_RB + M_A is not positive definite")
            return Craft(mass_matrix=mass_matrix,
                         coriolis_matrix=known["rigid_body_coriolis_matrix"] + known["added_mass_coriolis_matrix"],
                         current_velocity=current_velocity, current_acceleration=current_acceleration, known=known)

        self._assembly = VehicleAssembly(context, "SubmergedBemAuv", compute, own_inputs,
                                         kinematics_form(context, "SubmergedBemAuv", CompositionError))

    # -- the plugin type's methods ---------------------------------------------------------------
    def graph(self):
        if self._assembly is None:
            raise RuntimeError("SubmergedBemAuv must be initialized before graph()")
        return self._assembly.payload()

    def getInputDescriptions(self):
        return list(self._assembly.input_descriptions)

    def signals(self, x, u=()):
        """Every named quantity of the vehicle and its children's outputs at a state and input (diagnostics)."""
        return self._assembly.signals(x, u)

    def step_vector(self, x, u, dt):
        return self._assembly.step_vector(x, u, dt)

    def step(self, state: VehicleModel3D.Odometry3D, command: List[VehicleModel3D.Command], t: float, dt: float,
             **kwargs) -> VehicleModel3D.Odometry3D:
        return self._assembly.step(state, command, t, dt)
