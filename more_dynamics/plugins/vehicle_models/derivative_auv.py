"""A torpedo-shaped autonomous underwater vehicle whose mass is a full
(non-diagonal) inertia tensor and whose added mass and hydrodynamic forces
are given together as published nondimensional derivative tables (Healey and
Lienard 1993, via MSS ``npsauv.m``), rather than ``SpheroidAuv``'s formula
parts (Lamb spheroid added mass, Prestero lift/drag) — the given-derivatives
route an identified vehicle's own model later takes (E-69, E-103). As Luka's
``VehicleModel3D`` in the form of ``HullVessel``/``SpheroidAuv``: its own
parameters, the rigid body and the added mass computed here from
``models/``, and the slots ``hydrostatics``, ``hydrodynamics`` and
``actuators`` taken from the composition.

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
vehicle's own recomputes it): ``gravity`` (given, as ``npsauv.m`` 105, not
derived from ``latitude``: unlike ``remus100.m``, this reference does not
compute it from a site), ``water_density``, ``length`` (``npsauv.m`` 105),
``body_mass`` and the full inertia tensor (``inertia_diagonal``,
``inertia_products``, ``npsauv.m`` 106, 110-111), ``body_center_of_gravity``
(``npsauv.m`` 106). ``weight`` and ``buoyancy`` are declared too (``npsauv.m``
108, the test contract): both equal ``body_mass * gravity`` for this
reference (neutral buoyancy, checked: 5443.425... * 9.81 = 53400.0 to the
file's precision) and are read by the hydrostatics child through its own
``mass``/``gravity`` couplings (``SubmergedRestoring``'s default
``buoyancy_method = "neutral"``), not recomputed here; they are declared on
the vehicle because ``npsauv.m`` gives them as its own primitives (``W``,
``B``), not because this vehicle's graph consumes them a second time.

Every published nondimensional derivative (93, ``nondim_<Name>``) is
declared here (rule 16: no omitted primitive; test contract). The 14
``*dot`` (acceleration) derivatives are this vehicle's own (the added-mass
step below); the 52 hull-only derivatives are declared a second time, with
the same values, on ``CoefficientHullLoads`` (which alone reads them: the
contract requires every published quantity named on the vehicle, while the
quantity it is used by stays the hydrodynamics child, as ``AuvHullLoads``'s
own derivatives come from a different named part); the remaining 27 couple to a fin deflection or
the propeller and are not yet read by any plugin of this vehicle (flag,
below).

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

Default values and their provenance are in ``DEFAULTS.md``; the underlying
paper itself is not read (needs-access, ``SOURCE.md``) — every default is
read from MSS's instrumented workspace, cited to its line.

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

from more_transformations.more_casadi_transformations import check_symmetric_positive_definite

from more_dynamics.models.added_mass.added_mass_parts import (
    GIVEN_DERIVATIVE_ADDED_MASS_NAMES, given_derivative_added_mass, given_derivative_added_mass_parameters)
from more_dynamics.models.coefficient_loads.coefficient_loads import HULL_ONLY_DERIVATIVE_NAMES
from more_dynamics.models.rigid_body.rigid_body_parts import FULL_TENSOR_RIGID_BODY_PARAMETERS, full_tensor_rigid_body
from more_dynamics.plugins.shared.vehicle_graph import (
    INTEGRATION_STEP, CompositionError, Craft, Parts, VehicleAssembly)
from more_dynamics.plugins.shared.vehicle_options import CURRENT_FORMS, current_inputs, select_site

# npsauv.m 126-156: every published nondimensional derivative, by line.
_ACTUATOR_COUPLED_DERIVATIVES = {
    "Xqds": -1, "Xqdb2": -1, "Xrdr": -1, "Xvdr": -1, "Xwds": -1, "Xwdb2": -1, "Xdsds": -1, "Xdrdr": -1,
    "Xqdsn": -1, "Xwdsn": -1, "Xdsdsn": -1, "Ydr": -1, "Zds": -1, "Zdb2": -1, "Zqn": -1, "Zwn": -1, "Zdsn": -1,
    "Kdb2": -1, "Kpn": -1, "Kprop": -1, "Mds": -1, "Mdb2": -1, "Mqn": -1, "Mwn": -1, "Mdsn": -1, "Ndr": -1,
    "Nprop": -1,
}
_DERIVATIVE_VALUES = {
    "Xpp": 7.0e-3, "Xqq": -1.5e-2, "Xrr": 4.0e-3, "Xpr": 7.5e-4, "Xudot": -7.6e-3, "Xwq": -2.0e-1,
    "Xvp": -3.0e-3, "Xvr": 2.0e-2, "Xqds": 2.5e-2, "Xqdb2": -1.3e-3, "Xrdr": -1.0e-3, "Xvv": 5.3e-2,
    "Xww": 1.7e-1, "Xvdr": 1.7e-3, "Xwds": 4.6e-2, "Xwdb2": 0.5e-2, "Xdsds": -1.0e-2, "Xdrdr": -1.0e-2,
    "Xqdsn": 2.0e-3, "Xwdsn": 3.5e-3, "Xdsdsn": -1.6e-3,
    "Ypdot": 1.2e-4, "Yrdot": 1.2e-3, "Ypq": 4.0e-3, "Yqr": -6.5e-3, "Yvdot": -5.5e-2, "Yp": 3.0e-3,
    "Yr": 3.0e-2, "Yvq": 2.4e-2, "Ywp": 2.3e-1, "Ywr": -1.9e-2, "Yv": -1.0e-1, "Yvw": 6.8e-2, "Ydr": 2.7e-2,
    "Zqdot": -6.8e-3, "Zpp": 1.3e-4, "Zpr": 6.7e-3, "Zrr": -7.4e-3, "Zwdot": -2.4e-1, "Zq": -1.4e-1,
    "Zvp": -4.8e-2, "Zvr": 4.5e-2, "Zw": -3.0e-1, "Zvv": -6.8e-2, "Zds": -7.3e-2, "Zdb2": -1.3e-2,
    "Zqn": -2.9e-3, "Zwn": -5.1e-3, "Zdsn": -1.0e-2,
    "Kpdot": -1.0e-3, "Krdot": -3.4e-5, "Kpq": -6.9e-5, "Kqr": 1.7e-2, "Kvdot": 1.2e-4, "Kp": -1.1e-2,
    "Kr": -8.4e-4, "Kvq": -5.1e-3, "Kwp": -1.3e-4, "Kwr": 1.4e-2, "Kv": 3.1e-3, "Kvw": -1.9e-1, "Kdb2": 0.0,
    "Kpn": -5.7e-4, "Kprop": 0.0,
    "Mqdot": -1.7e-2, "Mpp": 5.3e-5, "Mpr": 5.0e-3, "Mrr": 2.9e-3, "Mwdot": -6.8e-3, "Muq": -6.8e-2,
    "Mvp": 1.2e-3, "Mvr": 1.7e-2, "Muw": 1.0e-1, "Mvv": -2.6e-2, "Mds": -4.1e-2, "Mdb2": 3.5e-3,
    "Mqn": -1.6e-3, "Mwn": -2.9e-3, "Mdsn": -5.2e-3,
    "Npdot": -3.4e-5, "Nrdot": -3.4e-3, "Npq": -2.1e-2, "Nqr": 2.7e-3, "Nvdot": 1.2e-3, "Np": -8.4e-4,
    "Nr": -1.6e-2, "Nvq": -1.0e-2, "Nwp": -1.7e-2, "Nwr": 7.4e-3, "Nv": -7.4e-3, "Nvw": -2.7e-2,
    "Ndr": -1.3e-2, "Nprop": 0.0,
}
_ALL_DERIVATIVE_NAMES = (GIVEN_DERIVATIVE_ADDED_MASS_NAMES + tuple(n for g in HULL_ONLY_DERIVATIVE_NAMES.values()
                         for n in g) + tuple(_ACTUATOR_COUPLED_DERIVATIVES))
assert len(_ALL_DERIVATIVE_NAMES) == 93 and set(_ALL_DERIVATIVE_NAMES) == set(_DERIVATIVE_VALUES)


def _check_inertia(context: ComponentContext) -> None:
    """Refuse a composition whose full inertia tensor about the CO
    (``inertia_diagonal``, ``inertia_products``, npsauv.m 159-161) is not
    symmetric positive definite — a wrong sign convention on
    ``inertia_products`` (rule 11 of the A-62 trap list: a composition error
    silently accepted) gives a physically impossible body."""
    diagonal = [float(v) for v in context.get_parameter("inertia_diagonal")]
    ixy, iyz, ixz = (float(v) for v in context.get_parameter("inertia_products"))
    inertia = [[diagonal[0], -ixy, -ixz], [-ixy, diagonal[1], -iyz], [-ixz, -iyz, diagonal[2]]]
    check_symmetric_positive_definite({"inertia": inertia}, "inertia")


class DerivativeAuv(VehicleModel3D):

    COMPONENTS = {
        "actuators": "List[more_dynamics::ForceProducer]",
        "hydrostatics": "more_dynamics::HydrostaticsModel",
        "hydrodynamics": "more_dynamics::HydrodynamicsModel",
    }

    PARAMETERS = [
        INTEGRATION_STEP,
        ParameterDescription("gravity", 9.81),  # m/s^2, npsauv.m:105 (given directly, not from a site)
        ParameterDescription("water_density", 1025),  # kg/m^3, npsauv.m:105 rho
        ParameterDescription("length", 5.3),  # m, npsauv.m:105 L
        ParameterDescription("body_mass", 5443.425076452599),  # kg, W / g (npsauv.m:109), W at npsauv.m:108
        ParameterDescription("inertia_diagonal", [2038, 13587, 13587]),  # kg*m^2, npsauv.m:110 [Ix, Iy, Iz]
        ParameterDescription("inertia_products", [-13.58, -13.58, -13.58]),  # kg*m^2, npsauv.m:111 [Ixy, Iyz, Ixz]
        ParameterDescription("body_center_of_gravity", [0, 0, 0.061]),  # m, npsauv.m:106 r_bg
        ParameterDescription("weight", 53400),  # N, npsauv.m:108 W (read by hydrostatics via mass * gravity)
        ParameterDescription("buoyancy", 53400),  # N, npsauv.m:108 B (neutral: SubmergedRestoring's own default)
        ParameterDescription("center_of_buoyancy", [0, 0, 0]),  # m, npsauv.m:107 r_bb
        ParameterDescription("derivative_scaling", "nondimensional"),  # npsauv.m:116-124 (Fossen prime-scaling)
        ParameterDescription("current_form", "yaw_rate_terms"),  # npsauv.m:83-90 (the gate's own MSS form)
        ParameterDescription("site_form", "given"),  # npsauv.m:105 g = 9.81 (no latitude in this reference)
        ParameterDescription("latitude", 0.0),  # rad; unread with site_form "given" (declared for the shared form)
        ParameterDescription("kinematic_viscosity", 1e-06),  # m^2/s; unread by this vehicle's own graph (site's own input)
        *[ParameterDescription(f"nondim_{name}", value) for name, value in _DERIVATIVE_VALUES.items()],
    ]

    def __init__(self):
        self._assembly = None

    def initialize(self, context: ComponentContext):
        form = context.get_parameter("current_form")
        own_inputs = current_inputs(context, "DerivativeAuv", CompositionError)
        _check_inertia(context)

        def compute(pose, velocity, current):
            known = {"pose": pose, "velocity": velocity, **current}
            parts = Parts(context, "DerivativeAuv", known)
            site, site_parameters = select_site(parts)  # (npsauv.m 105; otter.m 90-91)
            parts.run("site", site(), site_parameters)
            known["length"] = ca.SX(context.get_parameter("length"))  # npsauv.m:105 L (vehicle coupling, not a site/hull-form quantity here)
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
