"""A torpedo-shaped autonomous underwater vehicle with a prolate-spheroid hull,
as Luka's ``VehicleModel3D`` in the form of ``HullVessel`` (``hull_vessel.py``):
its own parameters (geometry, mass, centres, forms), the rigid body, the
added mass, the Coriolis terms and the current computed here from the model
functions of ``models/``, and the slots ``hydrostatics``, ``hydrodynamics``
and ``actuators`` (a list of force producers) taken from the composition. The
graph is the equation of motion of a craft in a current over the children's
payloads (``plugins/shared/vehicle_graph.py``).

    nu_r    = nu - nu_c
    nu_dot  = nu_c_dot + M^-1 ( tau - C nu_r )      M = M_RB + M_A,  C = C_RB(nu) + C_A(nu_r)
    eta_dot = J(eta) nu

State ``[eta; nu]`` after the actuators' states: NED position, ZYX Euler
angles, BODY FRD velocity, SI units. Inputs: the current (``current_speed``,
``current_direction``, ``current_vertical_speed``) when ``current_form`` is
not ``"none"``, then the commands of the actuators in list order.

Primitives (rule: derived quantities are computed): ``latitude`` (gravity by
WGS-84 normal gravity), ``water_density``, ``kinematic_viscosity``, ``length``
and ``diameter`` (the hull, one geometry: beam, draft, span and the semi-axes
follow), ``body_density`` and ``body_center_of_gravity`` (mass and inertia of
the homogeneous spheroid), ``roll_added_inertia_ratio`` (Lamb added mass).

Options, each with its fidelity and what raises it:

* ``added_mass_form``: ``"lamb_spheroid"`` (the ideal fluid around a prolate
  spheroid, Fossen 2011 ch. 6, imlay61.m). Raise: a Capytaine or identified
  matrix.
* ``coriolis_form``: ``"kirchhoff_full"`` (every term of ``C_A``, the physics
  form, default) or ``"munk_couplings_removed"`` (MSS shortcut of
  ``remus100.m`` 207-210, a comparison form).
* ``current_form``: ``"none"`` (still water), ``"yaw_rate_terms"`` (MSS
  ``remus100.m`` 118-122) or ``"full_rotation_rate"`` (MSS ``otter.m``
  113-118). Raise: a measured current profile.

Defaults are the REMUS 100 (Prestero 2001; MSS ``remus100.m``), one density
and one geometry throughout; each default cites its line.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 2.40, p. 26; eq. 6.48, p. 120.
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis.
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 3, 96-98,
    118-125, 131-138, 207-210, 257-259; LIBRARY/modeling/spheroid.m 35-42;
    imlay61.m 31-59; INS/functions/gravity.m 11-12.
[HullVessel] Mandic, L. more_dynamics ``HullVessel`` (hull_vessel.py): the
    form of the vehicle plugin (parameters, slots, ``graph()``) this follows.

Author:    Enio Krizman
Date:      2026-10-09
"""
from typing import List

import casadi as ca

from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.added_mass.added_mass_parts import lamb_spheroid, lamb_spheroid_parameters
from more_dynamics.models.coriolis.added_mass_coriolis_parts import kirchhoff_full, munk_couplings_removed
from more_dynamics.models.current.current import (
    horizontal_current_full_rotation_rate, horizontal_current_yaw_rate_terms, no_current)
from more_dynamics.models.hull_form.hull_form import (
    PROLATE_SPHEROID_MAIN_DIMENSIONS_PARAMETERS, prolate_spheroid_main_dimensions)
from more_dynamics.models.rigid_body.rigid_body_parts import homogeneous_spheroid, homogeneous_spheroid_parameters
from more_dynamics.models.site.site import SITE_AT_LATITUDE_PARAMETERS, site_at_latitude
from more_dynamics.plugins.shared.vehicle_graph import (
    INTEGRATION_STEP, CompositionError, Craft, Parts, VehicleAssembly)

ADDED_MASS_FORMS = {"lamb_spheroid": (lamb_spheroid, lamb_spheroid_parameters)}
CORIOLIS_FORMS = {"kirchhoff_full": kirchhoff_full, "munk_couplings_removed": munk_couplings_removed}
CURRENT_FORMS = {"none": no_current, "yaw_rate_terms": horizontal_current_yaw_rate_terms,
                 "full_rotation_rate": horizontal_current_full_rotation_rate}
CURRENT_INPUTS = (("current_speed", 1, "horizontal current speed, m/s"),
                  ("current_direction", 1, "current set, NED, from north, clockwise, rad"),
                  ("current_vertical_speed", 1, "vertical current speed, NED, down positive, m/s"))


class SpheroidAuv(VehicleModel3D):

    COMPONENTS = {
        "actuators": "List[more_dynamics::ForceProducer]",
        "hydrostatics": "more_dynamics::HydrostaticsModel",
        "hydrodynamics": "more_dynamics::HydrodynamicsModel",
    }

    PARAMETERS = [
        INTEGRATION_STEP,
        ParameterDescription("latitude", 1.1073560310932362),  # remus100.m 96-97, mu = deg2rad(63.446827)
        ParameterDescription("water_density", 1026),  # remus100.m:98 rho = 1026, the one value (MSS: 1026 at imlay61.m:31 and forceLiftDrag.m:26, 1025 at crossFlowDrag.m:36)
        ParameterDescription("kinematic_viscosity", 1e-06),  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
        ParameterDescription("length", 1.6),  # remus100.m:131 L_auv, passed as L at :221
        ParameterDescription("diameter", 0.19),  # remus100.m:132 D_auv, passed as B at :221
        ParameterDescription("body_density", 1054.7872613500267),  # mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: 1025 at spheroid.m:35
        ParameterDescription("body_center_of_gravity", [0, 0, 0.02]),  # remus100.m:137 r_bG
        ParameterDescription("roll_added_inertia_ratio", 0.3),  # remus100.m:136 r44
        ParameterDescription("added_mass_form", "lamb_spheroid"),  # imlay61.m 31-59
        ParameterDescription("coriolis_form", "kirchhoff_full"),  # m2c.m 33-48 (every term kept)
        ParameterDescription("current_form", "none"),  # still water
    ]

    def __init__(self):
        self._assembly = None

    def initialize(self, context: ComponentContext):
        form = context.get_parameter("current_form")
        if form not in CURRENT_FORMS:
            raise CompositionError(f"SpheroidAuv: current_form must be one of {sorted(CURRENT_FORMS)}, got {form!r}")
        current_inputs = () if form == "none" else CURRENT_INPUTS

        def compute(pose, velocity, current):
            known = {"pose": pose, "velocity": velocity, **current}
            parts = Parts(context, "SpheroidAuv", known)
            parts.run("site", site_at_latitude(), SITE_AT_LATITUDE_PARAMETERS)  # (gravity.m 11-12; remus100.m 96-97)
            parts.run("hull form", prolate_spheroid_main_dimensions(), PROLATE_SPHEROID_MAIN_DIMENSIONS_PARAMETERS)
            parts.run("rigid body", homogeneous_spheroid(), homogeneous_spheroid_parameters())  # (spheroid.m 35-42)
            added_mass, added_mass_parameters = parts.select("added_mass_form", ADDED_MASS_FORMS)
            parts.run("added mass", added_mass(), added_mass_parameters())
            parts.run("current", CURRENT_FORMS[form]())
            current_velocity, current_acceleration = known["current_velocity"], known["current_acceleration"]
            known["relative_velocity"] = velocity - current_velocity  # nu_r = nu - nu_c (remus100.m 125)
            parts.run("added-mass Coriolis", parts.select("coriolis_form", CORIOLIS_FORMS)())
            return Craft(mass_matrix=known["rigid_body_mass_matrix"] + known["added_mass_matrix"],
                         coriolis_matrix=known["rigid_body_coriolis_matrix"] + known["added_mass_coriolis_matrix"],
                         current_velocity=current_velocity, current_acceleration=current_acceleration, known=known)

        self._assembly = VehicleAssembly(context, "SpheroidAuv", compute, current_inputs)

    # -- the plugin type's methods ---------------------------------------------------------------
    def graph(self):
        if self._assembly is None:
            raise RuntimeError("SpheroidAuv must be initialized before graph()")
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
