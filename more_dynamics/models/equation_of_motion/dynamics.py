"""Equation of motion of a 6-DOF marine craft in a current: the vehicle's own
graph, built from the inertia, Coriolis and force quantities its parts give.
The model layer of ``MarineCraft6DOF`` (the way ``models/hull_vessel/
dynamics.py`` is the model layer of ``HullVessel``); it imports no plugin.

    nu_r    = nu - nu_c
    nu_dot  = nu_c_dot + M^-1 ( tau - C nu_r )      M = M_RB + M_A,  C = C_RB(nu) + C_A(nu_r)
    eta_dot = J(eta) nu

or, with unit quaternions (``marine_craft_6dof_quaternion_casadi``), the kinematics of the state
``[p; q; nu]`` (NED position, unit quaternion ``q = [eta, eps1, eps2, eps3]``, BODY FRD velocity)

    p_dot = R(q) v,   q_dot = T(q) omega        (13 states, the same nu_dot)

``tau`` is the sum of every signed generalized force the parts add (restoring
included, as in ``vehicle_model.capnp`` 11-21). State ``[eta; nu]``: NED
position, ZYX Euler angles, BODY FRD velocity, SI units.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 2.40, p. 26 (eta_dot = J(eta) nu, through
    more_transformations' ``EulerNEDBodyTransforms.body_to_ned_kinematics``);
    eq. 2.69, p. 30 (p_dot = R(q) v, q_dot = T(q) omega, through
    ``QuatNEDBodyTransforms.body_to_ned_kinematics_quat``);
    eq. 6.48, p. 120 (M = M_RB + M_A).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 257-259,
    CRAFT/USV/models/otter.m 261-263 (relative velocity and the current's
    body acceleration in the equation of motion).
[MSS 13] Fossen, T. I. MSS, MIT, @ ed2ce7b: CRAFT/AUV/models/remus100.m 225-227, 258-259 (the 13-state form of the
    same equation: ``quatern`` for J, ``xdot = [nu_dot; J nu]``).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import EulerNEDBodyTransforms, QuatNEDBodyTransforms


def marine_craft_6dof_casadi() -> ca.Function:
    """``(state, current_velocity, current_acceleration, mass_matrix,
    coriolis_matrix, force) -> state_dot``."""
    state = ca.SX.sym("state", 12)
    current_velocity = ca.SX.sym("current_velocity", 6)
    current_acceleration = ca.SX.sym("current_acceleration", 6)
    mass_matrix = ca.SX.sym("mass_matrix", 6, 6)
    coriolis_matrix = ca.SX.sym("coriolis_matrix", 6, 6)
    force = ca.SX.sym("force", 6)
    eta, nu = state[0:6], state[6:12]
    nu_r = nu - current_velocity  # (remus100.m 125; otter.m 116)
    nu_dot = current_acceleration + ca.solve(mass_matrix, force - coriolis_matrix @ nu_r)  # (remus100.m 257-258; otter.m 261-262)
    eta_dot = EulerNEDBodyTransforms.body_to_ned_kinematics(eta, nu)  # (Fossen 2011, eq. 2.40, p. 26; remus100.m 259)
    return ca.Function("marine_craft_6dof",
                       [state, current_velocity, current_acceleration, mass_matrix, coriolis_matrix, force],
                       [ca.vertcat(eta_dot, nu_dot)],
                       ["state", "current_velocity", "current_acceleration", "mass_matrix", "coriolis_matrix",
                        "force"],
                       ["state_dot"])


def marine_craft_6dof_quaternion_casadi() -> ca.Function:
    """The same as ``marine_craft_6dof_casadi`` for the state ``[p; q; nu]`` of 13 values: ``(state,
    current_velocity, current_acceleration, mass_matrix, coriolis_matrix, force) -> state_dot`` (13 values)."""
    state = ca.SX.sym("state", 13)
    current_velocity = ca.SX.sym("current_velocity", 6)
    current_acceleration = ca.SX.sym("current_acceleration", 6)
    mass_matrix = ca.SX.sym("mass_matrix", 6, 6)
    coriolis_matrix = ca.SX.sym("coriolis_matrix", 6, 6)
    force = ca.SX.sym("force", 6)
    q, nu = state[3:7], state[7:13]
    nu_r = nu - current_velocity  # (remus100.m 125)
    nu_dot = current_acceleration + ca.solve(mass_matrix, force - coriolis_matrix @ nu_r)  # (remus100.m 257-258)
    p_dot_q_dot = QuatNEDBodyTransforms.body_to_ned_kinematics_quat(q, nu)  # (Fossen 2011, eq. 2.69, p. 30; remus100.m 225-227, 259)
    return ca.Function("marine_craft_6dof_quaternion",
                       [state, current_velocity, current_acceleration, mass_matrix, coriolis_matrix, force],
                       [ca.vertcat(p_dot_q_dot, nu_dot)],
                       ["state", "current_velocity", "current_acceleration", "mass_matrix", "coriolis_matrix",
                        "force"],
                       ["state_dot"])
