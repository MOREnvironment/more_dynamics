"""A single-hull surface craft (a monohull), as Luka's ``VehicleModel3D`` in the
form of ``HullVessel`` (``hull_vessel.py``): its own parameters (the hull, the
mass with a centre of gravity and radii of gyration, the added-mass
derivatives), the rigid body, the added mass, the Coriolis terms and the
current computed here from the model functions of ``models/``, and the slots
``hydrostatics``, ``hydrodynamics`` and ``actuators`` (a list of force
producers) taken from the composition. The graph is the equation of motion of
a craft in a current over the children's payloads
(``plugins/shared/vehicle_graph.py``).

    nu_r    = nu - nu_c
    nu_dot  = nu_c_dot + M^-1 ( tau - C nu_r )      M = M_RB + M_A,  C = C_RB(nu) + C_A(nu_r)
    eta_dot = J(eta) nu

State ``[eta; nu]`` after the actuators' states: NED position, ZYX Euler
angles, BODY FRD velocity, SI units. Inputs: the current (``current_speed``,
``current_direction``, ``current_vertical_speed``) when ``current_form`` is
not ``"none"``, then the commands of the actuators in list order.

The hull of one body is the ``single_hull`` form (length, beam, block and
waterplane coefficients; the restoring reads it as ``hull_beam``,
``hull_block_coefficient``, ``hull_waterplane_coefficient``, the cross-flow
strip as ``section_beam``). The rigid body is the hull mass plus an optional
point payload (``payload_mass`` 0 is a bare hull), the added mass the scaled
derivatives of the same form as the catamaran (MSS ``otter.m`` 152-159), and
the options are those of the catamaran: ``added_mass_form``,
``coriolis_form``, ``current_form``, each with its fidelity.

Defaults are those of the NTNU Mariner 5 USV Grethe (Maritime Robotics build
40402, Pioner 17 ft hull), from the Grethe parameter file
``scripts/vehicle_models/params/grethe_mariner5.yaml`` (the key is in each
comment). **Document** values are read off a cited page; **estimate** values
are starting points with no document behind them (earlier Grethe parameter
files, unverified): to be identified from logs, and none is a measurement.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 2.40, p. 26; eq. 6.48, p. 120.
[F003] Maritime Robotics, NTNU Mariner 5 USV Combined FAT/SAT, rev 1.5, signed
    8 Apr 2021. p. 2 (length, beam, draft, weight of the boat).
[M001] Maritime Robotics, Mariner 5 USV User Manual rev 1.0, Apr 2021. p. 7
    (hull dimensions, mass).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/USV/models/otter.m 113-118,
    121-128, 152-159, 261-263 (the forms the equations follow).
[HullVessel] Mandic, L. more_dynamics ``HullVessel`` (hull_vessel.py): the
    form of the vehicle plugin (parameters, slots, ``graph()``) this follows.

Author:    Enio Krizman
Date:      2026-10-09
"""
from typing import List

from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.added_mass.added_mass_parts import SCALED_DERIVATIVES_PARAMETERS, scaled_derivatives
from more_dynamics.models.coriolis.added_mass_coriolis_parts import kirchhoff_full, munk_couplings_removed
from more_dynamics.models.current.current import (
    horizontal_current_full_rotation_rate, horizontal_current_yaw_rate_terms, no_current)
from more_dynamics.models.hull_form.hull_form import SINGLE_HULL_PARAMETERS, single_hull
from more_dynamics.models.rigid_body.rigid_body_parts import hull_with_point_payload, hull_with_point_payload_parameters
from more_dynamics.models.site.site import SITE_GIVEN_PARAMETERS, site_given
from more_dynamics.plugins.shared.vehicle_graph import (
    INTEGRATION_STEP, CompositionError, Craft, Parts, VehicleAssembly)

ADDED_MASS_FORMS = {"scaled_derivatives": (scaled_derivatives, lambda: SCALED_DERIVATIVES_PARAMETERS)}
CORIOLIS_FORMS = {"kirchhoff_full": kirchhoff_full, "munk_couplings_removed": munk_couplings_removed}
CURRENT_FORMS = {"none": no_current, "yaw_rate_terms": horizontal_current_yaw_rate_terms,
                 "full_rotation_rate": horizontal_current_full_rotation_rate}
CURRENT_INPUTS = (("current_speed", 1, "horizontal current speed, m/s"),
                  ("current_direction", 1, "current set, NED, from north, clockwise, rad"),
                  ("current_vertical_speed", 1, "vertical current speed, NED, down positive, m/s"))


class Monohull(VehicleModel3D):

    COMPONENTS = {
        "actuators": "List[more_dynamics::ForceProducer]",
        "hydrostatics": "more_dynamics::HydrostaticsModel",
        "hydrodynamics": "more_dynamics::HydrodynamicsModel",
    }

    PARAMETERS = [
        INTEGRATION_STEP,
        ParameterDescription("gravity", 9.81),  # estimate: grethe_mariner5.yaml hydrostatics.gravity ("standard")
        ParameterDescription("water_density", 1025.0),  # estimate: grethe_mariner5.yaml hydrostatics.water_density ("standard sea water")
        ParameterDescription("kinematic_viscosity", 1e-06),  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
        ParameterDescription("length", 5.2),  # document: grethe_mariner5.yaml hull.length (F003 p2; M001 p7)
        ParameterDescription("beam", 2.15),  # document: grethe_mariner5.yaml hull.beam (F003 p2; M001 p7; D002 p1)
        ParameterDescription("hull_block_coefficient", 0.233),  # estimate: grethe_mariner5.yaml hydrostatics.block_coefficient (a mass knob: 806 / (1025 * 5.2 * 2.15 * 0.3) = 0.234)
        ParameterDescription("hull_waterplane_coefficient", 0.8),  # estimate: grethe_mariner5.yaml hydrostatics.waterplane_coefficient
        ParameterDescription("hull_mass", 806.0),  # document: grethe_mariner5.yaml mass.mass_delivered (F003 p2); the field mass (USBL mount, payload) is a gap
        ParameterDescription("payload_mass", 0.0),  # a bare hull; the field payload is a gap (grethe_mariner5.yaml mass.mass_field)
        ParameterDescription("hull_center_of_gravity", [0.0, 0.0, 0.025]),  # estimate: grethe_mariner5.yaml mass.center_of_gravity
        ParameterDescription("payload_position", [0.0, 0.0, 0.0]),  # unused while payload_mass is 0
        ParameterDescription("radii_of_gyration", [0.35, 0.25, 0.25]),  # estimate: grethe_mariner5.yaml mass.radii_of_gyration_scale (fractions of [B, L, L])
        ParameterDescription("added_mass_coefficients", [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]),  # estimate: grethe_mariner5.yaml mass.added_mass_scales workspace_rppws_value (open: the library value is -1.2 in yaw)
        ParameterDescription("added_mass_form", "scaled_derivatives"),  # otter.m 152-159
        ParameterDescription("coriolis_form", "kirchhoff_full"),  # m2c.m 33-48 (every term kept)
        ParameterDescription("current_form", "none"),  # still water
    ]

    def __init__(self):
        self._assembly = None

    def initialize(self, context: ComponentContext):
        form = context.get_parameter("current_form")
        if form not in CURRENT_FORMS:
            raise CompositionError(f"Monohull: current_form must be one of {sorted(CURRENT_FORMS)}, got {form!r}")
        current_inputs = () if form == "none" else CURRENT_INPUTS

        def compute(pose, velocity, current):
            known = {"pose": pose, "velocity": velocity, **current}
            parts = Parts(context, "Monohull", known)
            parts.run("site", site_given(), SITE_GIVEN_PARAMETERS)
            parts.run("hull form", single_hull(), SINGLE_HULL_PARAMETERS)
            parts.run("rigid body", hull_with_point_payload(), hull_with_point_payload_parameters())  # (otter.m 94-98, 122-128)
            added_mass, added_mass_parameters = parts.select("added_mass_form", ADDED_MASS_FORMS)
            parts.run("added mass", added_mass(), added_mass_parameters())  # (otter.m 152-159)
            parts.run("current", CURRENT_FORMS[form]())
            current_velocity, current_acceleration = known["current_velocity"], known["current_acceleration"]
            known["relative_velocity"] = velocity - current_velocity  # nu_r = nu - nu_c (otter.m 116)
            parts.run("added-mass Coriolis", parts.select("coriolis_form", CORIOLIS_FORMS)())
            return Craft(mass_matrix=known["rigid_body_mass_matrix"] + known["added_mass_matrix"],
                         coriolis_matrix=known["rigid_body_coriolis_matrix"] + known["added_mass_coriolis_matrix"],
                         current_velocity=current_velocity, current_acceleration=current_acceleration, known=known)

        self._assembly = VehicleAssembly(context, "Monohull", compute, current_inputs)

    # -- the plugin type's methods ---------------------------------------------------------------
    def graph(self):
        if self._assembly is None:
            raise RuntimeError("Monohull must be initialized before graph()")
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
