"""A torpedo-shaped autonomous underwater vehicle whose mass is a full
(non-diagonal) inertia tensor and whose added mass and hydrodynamic forces
are given together as published nondimensional derivative tables (Healey and
Lienard 1993, via MSS ``npsauv.m``), rather than ``SpheroidAuv``'s formula
parts (Lamb spheroid added mass, Prestero lift/drag) — the given-derivatives
route an identified vehicle's own model later takes (E-69, E-103). As Luka's
``VehicleModel3D`` in the form of ``HullVessel``/``SpheroidAuv``: its own
parameters, the rigid body and the added mass computed here from
``models/``, and the slots ``site`` (the water and the gravity),
``hydrostatics``, ``hydrodynamics`` and ``actuators`` taken from the
composition.

**No rigid-body Coriolis matrix of its own.** Unlike ``remus100.m``/
``otter.m``, ``npsauv.m`` applies no separate ``-C(nu) nu_r`` term in its
equation of motion: the rigid-body Coriolis-centripetal force, evaluated
with the relative velocity, is written directly into the hull's
nondimensional-derivative force (``npsauv.m`` 241-274;
``models/coefficient_loads`` module docstring, which also flags a candidate
MSS sign inconsistency, reproduced bug-for-bug, not fixed). This vehicle's
``coriolis_matrix`` is therefore the zero matrix; the hydrodynamics slot
(``CoefficientHullLoads``) carries that physics instead, so it is not
counted twice.

    nu_r    = nu - nu_c
    nu_dot  = nu_c_dot + M^-1 ( tau )                M = M_RB + M_A
    eta_dot = J(eta) nu

State ``[eta; nu]`` after the actuators' states: NED position, ZYX Euler
angles, BODY FRD velocity, SI units. Inputs: the current
(``current_speed``, ``current_direction``, ``current_vertical_speed``) when
``current_form`` is not ``"none"``, then the commands of the actuators in
list order.

Primitives (rule 16: derived quantities are computed, except where
``npsauv.m`` gives a quantity directly and no defining equation of this
vehicle's own recomputes it): ``length`` (``npsauv.m`` 105),
``body_mass`` and the full inertia tensor (``inertia_diagonal``,
``inertia_products``, ``npsauv.m`` 106, 110-111), ``body_center_of_gravity``
(``npsauv.m`` 106). ``weight`` and ``buoyancy`` are declared too (``npsauv.m``
108, the test contract): both equal ``body_mass * gravity`` for this
reference (neutral buoyancy, checked: 5443.425... * 9.81 = 53400.0 to the
file's precision, the gravity of the ``site`` slot being 9.81 as ``npsauv.m`` 105 gives it) and are read by the hydrostatics child through its own
``mass``/``gravity`` couplings (``SubmergedRestoring``'s default
``buoyancy_method = "neutral"``), not recomputed here; they are declared on
the vehicle because ``npsauv.m`` gives them as its own primitives (``W``,
``B``), not because this vehicle's graph consumes them a second time.

Of the 93 published nondimensional derivatives (``nondim_<Name>``), this
vehicle declares only the 14 ``*dot`` (acceleration) derivatives: they are its
own (the added-mass step below). The 52 hull-only derivatives are declared on
``CoefficientHullLoads``, which alone reads them, and the remaining 27 couple
to a fin deflection or the propeller and are declared by no plugin until an
actuator design gives them a reader (flag, below); a derivative is declared
once, by the part that reads it. All 93 stay in the vehicle's parameter data
file and its ``SOURCE.md``.

**Flag (rule 21, not fixed, not guessed): the 27 actuator-coupled
derivatives** (e.g. ``Xqds``, ``Ydr``, ``Zqn``\\ ``*epsilon``) are
``npsauv.m``'s control-surface forces (``tau_control``, ``npsauv.m``
176-212): empirical polynomial terms in the deflection and the vehicle's own
relative velocity, several also scaled by a through-water correction
``epsilon`` shared across the stern plane, bow planes and roll moment, and by
the propeller's own sign. ``Fin``'s lift/drag model is a different physics
form (a lifting surface evaluated at the local flow, Prestero 2001) built for
the vehicles that use the ``Fin`` force producer; it does not reproduce this empirical polynomial, and the shared
``epsilon`` is a coupling between an actuator (a fin) and another actuator
(the propeller) that the existing architecture's actuator slot does not
carry (each force producer is independent, ``vehicle_graph.py``). Four
``Fin`` and one ``Propeller`` instance are composed here as the brief names
(their servo lag correctly reproduces the actuator time constant and the
lag-gate test), but their net force does not reproduce ``npsauv.m``'s
``tau_control`` term by term — the full-model G1/G3/G4 gates are not green.
This is a gap in the standard's actuator slot for this vehicle's physics
(rule 19): the smallest extension in the standard's own shape is not
designed here; it is an owner question for the next job.

Every vehicle number is unset in the plugin and comes from a named part
(``DEFAULTS.md`` gives each reference part's provenance); the underlying
paper itself is not read (needs-access, ``SOURCE.md``) — every value is read
from MSS's instrumented workspace, cited to its line.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/npsauv.m 104-172.
[HullVessel] Mandic, L. more_dynamics ``HullVessel`` (hull_vessel.py): the
    form of the vehicle plugin this follows.
[SpheroidAuv] Krizman, E. more_dynamics ``SpheroidAuv`` (spheroid_auv.py):
    the first instance of this form in this library.

Author:    Enio Krizman
Date:      2026-10-09
"""
from typing import List

import casadi as ca

from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import check_symmetric_positive_definite, unset

from more_dynamics.models.added_mass.added_mass_parts import (
    GIVEN_DERIVATIVE_ADDED_MASS_NAMES, given_derivative_added_mass, given_derivative_added_mass_parameters)
from more_dynamics.models.rigid_body.rigid_body_parts import FULL_TENSOR_RIGID_BODY_PARAMETERS, full_tensor_rigid_body
from more_dynamics.plugins.shared.payload_io import read_value, refuses_unset
from more_dynamics.plugins.shared.vehicle_graph import (
    INTEGRATION_STEP, CompositionError, Craft, Parts, VehicleAssembly)
from more_dynamics.plugins.shared.vehicle_options import CURRENT_FORMS, current_inputs

def _check_inertia(context: ComponentContext) -> None:
    """Refuse a composition whose full inertia tensor about the CO
    (``inertia_diagonal``, ``inertia_products``, npsauv.m 159-161) is not
    symmetric positive definite — a wrong sign convention on
    ``inertia_products`` (rule 11 of the A-62 trap list: a composition error
    silently accepted) gives a physically impossible body."""
    diagonal = [float(v) for v in read_value(context, "inertia_diagonal", "kg*m^2",
                                             "principal moments of inertia [Ix, Iy, Iz] about the CO",
                                             placeholder=[1.0, 1.0, 1.0])]
    ixy, iyz, ixz = (float(v) for v in read_value(context, "inertia_products", "kg*m^2",
                                                  "products of inertia [Ixy, Iyz, Ixz] about the CO",
                                                  placeholder=[0.0, 0.0, 0.0]))
    inertia = [[diagonal[0], -ixy, -ixz], [-ixy, diagonal[1], -iyz], [-ixz, -iyz, diagonal[2]]]
    check_symmetric_positive_definite({"inertia": inertia}, "inertia")


class DerivativeAuv(VehicleModel3D):

    COMPONENTS = {
        "actuators": "List[more_dynamics::ForceProducer]",
        "site": "more_dynamics::SiteModel",
        "hydrostatics": "more_dynamics::HydrostaticsModel",
        "hydrodynamics": "more_dynamics::HydrodynamicsModel",
    }

    PARAMETERS = [
        INTEGRATION_STEP,
        ParameterDescription("length", unset()),  # m, npsauv.m:105 L
        ParameterDescription("body_mass", unset()),  # kg, W / g (npsauv.m:109), W at npsauv.m:108
        ParameterDescription("inertia_diagonal", unset()),  # kg*m^2, npsauv.m:110 [Ix, Iy, Iz]
        ParameterDescription("inertia_products", unset()),  # kg*m^2, npsauv.m:111 [Ixy, Iyz, Ixz]
        ParameterDescription("body_center_of_gravity", unset()),  # m, npsauv.m:106 r_bg
        ParameterDescription("weight", unset()),  # N, npsauv.m:108 W (read by hydrostatics via mass * gravity)
        ParameterDescription("buoyancy", unset()),  # N, npsauv.m:108 B (neutral: the hydrostatics child's own method)
        ParameterDescription("center_of_buoyancy", unset()),  # m, npsauv.m:107 r_bb
        ParameterDescription("derivative_scaling", "nondimensional"),  # npsauv.m:116-124 (Fossen prime-scaling)
        ParameterDescription("current_form", "yaw_rate_terms"),  # npsauv.m:83-90 (the gate's own MSS form)
        *[ParameterDescription(f"nondim_{name}", unset()) for name in GIVEN_DERIVATIVE_ADDED_MASS_NAMES],  # npsauv.m 166-171
    ]

    def __init__(self):
        self._assembly = None

    @refuses_unset
    def initialize(self, context: ComponentContext):
        form = context.get_parameter("current_form")
        own_inputs = current_inputs(context, "DerivativeAuv", CompositionError)
        _check_inertia(context)

        def compute(pose, velocity, current, site):
            known = {"pose": pose, "velocity": velocity, **current, **site}
            parts = Parts(context, "DerivativeAuv", known)
            known["length"] = ca.SX(read_value(context, "length", "m", "overall length L"))  # npsauv.m:105 L (vehicle coupling, not a site/hull-form quantity here)
            parts.run("rigid body", full_tensor_rigid_body(), FULL_TENSOR_RIGID_BODY_PARAMETERS)
            parts.run("added mass", given_derivative_added_mass(), given_derivative_added_mass_parameters())
            parts.run("current", CURRENT_FORMS[form]())
            current_velocity, current_acceleration = known["current_velocity"], known["current_acceleration"]
            # No rigid-body Coriolis matrix: npsauv.m folds it into the hydrodynamics child's force (module docstring).
            return Craft(mass_matrix=known["rigid_body_mass_matrix"] + known["added_mass_matrix"],
                         coriolis_matrix=ca.SX.zeros(6, 6),
                         current_velocity=current_velocity, current_acceleration=current_acceleration, known=known)

        self._assembly = VehicleAssembly(context, "DerivativeAuv", compute, own_inputs)

    # -- the plugin type's methods ---------------------------------------------------------------
    def graph(self):
        if self._assembly is None:
            raise RuntimeError("DerivativeAuv must be initialized before graph()")
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
