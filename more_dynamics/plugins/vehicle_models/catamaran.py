"""A twin-hull surface craft (a catamaran), as Luka's ``VehicleModel3D`` in the
form of ``HullVessel`` (``hull_vessel.py``): its own parameters (the two
pontoons, the hull mass and a point payload, the added-mass derivatives), the
rigid body, the added mass, the Coriolis terms and the current computed here
from the model functions of ``models/``, and the slots ``hydrostatics``,
``hydrodynamics`` and ``actuators`` (a list of force producers) taken from the
composition. The graph is the equation of motion of a craft in a current over
the children's payloads (``plugins/shared/vehicle_graph.py``).

    nu_r    = nu - nu_c
    nu_dot  = nu_c_dot + M^-1 ( tau - C nu_r )      M = M_RB + M_A,  C = C_RB(nu) + C_A(nu_r)
    eta_dot = J(eta) nu

State ``[eta; nu]`` after the actuators' states: NED position, ZYX Euler
angles, BODY FRD velocity, SI units. Inputs: the current (``current_speed``,
``current_direction``, ``current_vertical_speed``) when ``current_form`` is
not ``"none"``, then the commands of the actuators in list order.

Primitives: ``gravity``, ``water_density``, ``kinematic_viscosity``; the hull
(``length``, ``beam``, ``pontoon_beam``, ``pontoon_lateral_offset``,
``pontoon_block_coefficient``, ``pontoon_waterplane_coefficient``); the mass
(``hull_mass``, ``hull_center_of_gravity``, ``radii_of_gyration`` as
fractions of the beam and the length, a point payload ``payload_mass`` at
``payload_position``) and ``added_mass_coefficients`` (the scaled
derivatives). The pontoons' geometry is read by the restoring and the hull
loads by name (``hull_beam``, ``hull_block_coefficient``,
``hull_waterplane_coefficient``, ``hull_lateral_offset``, ``section_beam``).

Options, each with its fidelity and what raises it:

* ``added_mass_form``: ``"scaled_derivatives"`` (``M_A = -diag(c [A11, m, m,
  I11, I22, I33])``, MSS ``otter.m`` 152-159; raise: a Capytaine or
  identified matrix).
* ``coriolis_form``: ``"kirchhoff_full"`` (every term of ``C_A``, default) or
  ``"munk_couplings_removed"`` (MSS shortcut, comparison form).
* ``current_form``: ``"none"``, ``"yaw_rate_terms"`` or
  ``"full_rotation_rate"`` (MSS ``otter.m`` 113-118). Raise: a measured
  current profile.

Defaults are the Otter (MSS ``otter.m``), each citing its line.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 2.40, p. 26; eq. 6.48, p. 120.
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/USV/models/otter.m 90-98,
    104-107, 113-118, 121-128, 152-159, 261-263.
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
from more_dynamics.models.hull_form.hull_form import TWIN_PONTOONS_PARAMETERS, twin_pontoons
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
HULL_ALIASES = {"hull_beam": "pontoon_beam", "hull_block_coefficient": "pontoon_block_coefficient",
                "hull_waterplane_coefficient": "pontoon_waterplane_coefficient",
                "hull_lateral_offset": "pontoon_lateral_offset"}


class Catamaran(VehicleModel3D):

    COMPONENTS = {
        "actuators": "List[more_dynamics::ForceProducer]",
        "hydrostatics": "more_dynamics::HydrostaticsModel",
        "hydrodynamics": "more_dynamics::HydrodynamicsModel",
    }

    PARAMETERS = [
        INTEGRATION_STEP,
        ParameterDescription("gravity", 9.81),  # otter.m 90-91
        ParameterDescription("water_density", 1025.0),  # otter.m 90-91
        ParameterDescription("kinematic_viscosity", 1e-06),  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
        ParameterDescription("length", 2.0),  # otter.m 92-93, 104-107
        ParameterDescription("beam", 1.08),  # otter.m 92-93, 104-107
        ParameterDescription("pontoon_beam", 0.25),  # otter.m 104-107
        ParameterDescription("pontoon_lateral_offset", 0.395),  # otter.m 104-107
        ParameterDescription("pontoon_block_coefficient", 0.4),  # otter.m 104-107
        ParameterDescription("pontoon_waterplane_coefficient", 0.75),  # otter.m 104-107
        ParameterDescription("hull_mass", 55.0),  # otter.m 94-98
        ParameterDescription("payload_mass", 25.0),  # generate_catamaran_mss.m 61-62 (payload mass)
        ParameterDescription("hull_center_of_gravity", [0.2, 0.0, -0.2]),  # otter.m 94-98
        ParameterDescription("payload_position", [0.05, 0.0, -0.35]),  # generate_catamaran_mss.m 61-62 (payload position)
        ParameterDescription("radii_of_gyration", [0.4, 0.25, 0.25]),  # otter.m 94-98
        ParameterDescription("added_mass_coefficients", [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]),  # otter.m 152-157
        ParameterDescription("added_mass_form", "scaled_derivatives"),  # otter.m 152-159
        ParameterDescription("coriolis_form", "kirchhoff_full"),  # m2c.m 33-48 (every term kept)
        ParameterDescription("current_form", "none"),  # still water
    ]

    def __init__(self):
        self._assembly = None

    def initialize(self, context: ComponentContext):
        form = context.get_parameter("current_form")
        if form not in CURRENT_FORMS:
            raise CompositionError(f"Catamaran: current_form must be one of {sorted(CURRENT_FORMS)}, got {form!r}")
        current_inputs = () if form == "none" else CURRENT_INPUTS

        def compute(pose, velocity, current):
            known = {"pose": pose, "velocity": velocity, **current}
            parts = Parts(context, "Catamaran", known)
            parts.run("site", site_given(), SITE_GIVEN_PARAMETERS)  # (otter.m 90-91)
            parts.run("hull form", twin_pontoons(), TWIN_PONTOONS_PARAMETERS)  # (otter.m 92-93, 104-107)
            for alias, name in HULL_ALIASES.items():
                known[alias] = known[name]
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

        self._assembly = VehicleAssembly(context, "Catamaran", compute, current_inputs)

    # -- the plugin type's methods ---------------------------------------------------------------
    def graph(self):
        if self._assembly is None:
            raise RuntimeError("Catamaran must be initialized before graph()")
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
