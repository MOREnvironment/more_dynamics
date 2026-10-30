"""Payload trim: the shift of the equilibrium of the linear restoring that a point payload off the centre of
flotation causes (CasADi).

A craft whose restoring is the linear one, ``g = G eta`` (Fossen 2011, eq. 4.25, p. 65), has its equilibrium at
``eta = 0`` for the displacement the hull was drawn for. A payload of mass ``m_P`` at ``r_P`` adds a weight that is
not on the vertical through the centre of flotation: it pushes the craft down and tilts it, and the equilibrium moves
to ``eta_0`` with ``G eta_0 = g_0`` over heave, roll and pitch (Fossen 2011, §4.4, eqs. 4.94-4.100, pp. 74-76: the
equilibrium of a craft carrying a weight ``W`` at ``r``, ``Z = W``, ``K = y W``, ``M = -x W``), and the restoring force is
``g = G (eta - eta_0)``.

``payload_trim_offset(form)`` builds the block ``(eta, restoring_matrix, payload_mass, payload_position, gravity)
-> (trim_offset, kinematic_attitude_offset)``::

    f_P  = R(eta)^T [0 0 m_P g]^T                      the payload weight in BODY      (otter.m 265)
    m_P' = r_P x f_P                                    its moment about the CO         (otter.m 266)
    g_0  = [f_P(3); m_P'(1); m_P'(2)]                   heave, roll, pitch              (otter.m 267-268)
    eta_0 = [0; 0; G(3:5, 3:5)^-1 g_0; 0]                                               (otter.m 271)

with the rotation at the actual attitude (otter.m 265 reads ``eta`` before it is shifted). The two forms differ in
what ``g_0`` holds in heave and in where the kinematics are taken:

* ``"lever_arm_moment"`` (physics form, default): the displacement of the hydrostatic block is the loaded one,
  ``nabla = (m + m_P) / rho``, so buoyancy already balances the weight of hull and payload at ``eta = 0``; the payload
  enters the equilibrium only through the moment of its lever arm. ``g_0 = [0; m_P'(1); m_P'(2)]``; the kinematics are
  evaluated at the actual attitude (``kinematic_attitude_offset`` zero). Raise: the payload weighed and its position
  measured, and the trim read against a measured static roll and pitch.
* ``"lever_arm_and_heave_shift"`` (MSS ``otter.m`` 265-275, a comparison form): ``g_0`` also holds the payload weight
  in heave, so the same weight enters the heave equation twice (once through ``nabla``, once through ``eta_0``: at
  ``eta = eta_0`` a net vertical force ``m_P g`` is left), and the kinematics are evaluated at the shifted attitude
  ``eta - eta_0`` (``otter.m`` 272, 275; ``kinematic_attitude_offset = eta_0``). Both are inconsistencies of that
  file; the form keeps them so that a craft built on it equals ``otter.m`` with a payload.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Eq. 4.25, p. 65
    (the linear restoring); §4.4, eqs. 4.94-4.100, pp. 74-76 (the equilibrium under an added weight).
[MSS] Fossen, T. I. MSS, MIT, @ 6e78e98: CRAFT/USV/models/otter.m 127-128 (loaded displacement), 259-269 (payload force
    and moment, trim condition, shifted attitude); LIBRARY/kinematics/Rzyx.m, Smtrx.m.

Author:    Enio Krizman
Date:      2026-10-26
"""
import casadi as ca

from more_transformations.more_casadi_transformations import MatrixTransforms

PAYLOAD_TRIM_FORMS = ("lever_arm_moment", "lever_arm_and_heave_shift")


def check_payload_trim_form(form):
    if form not in PAYLOAD_TRIM_FORMS:
        raise ValueError(f"payload_trim_form must be one of {PAYLOAD_TRIM_FORMS}, got {form!r}")


def payload_trim_offset(form) -> ca.Function:
    """``(eta, restoring_matrix, payload_mass, payload_position, gravity) -> (trim_offset, kinematic_attitude_offset)``
    for ``form`` (one of ``PAYLOAD_TRIM_FORMS``; the module docstring): the shift ``eta_0`` (6, zero in surge, sway and
    yaw) of the equilibrium of ``g = G (eta - eta_0)``, and the shift of the attitude the kinematics are evaluated
    at (``eta_0`` for the MSS form, zero otherwise)."""
    check_payload_trim_form(form)
    eta = ca.SX.sym("eta", 6)
    restoring_matrix = ca.SX.sym("restoring_matrix", 6, 6)
    payload_mass = ca.SX.sym("payload_mass")
    payload_position = ca.SX.sym("payload_position", 3)
    gravity = ca.SX.sym("gravity")
    weight_body = MatrixTransforms.Rzyx(eta[3:6]).T @ ca.vertcat(0.0, 0.0, payload_mass * gravity)  # (otter.m 265)
    moment = ca.cross(payload_position, weight_body)  # (otter.m 266; Smtrx(r_P) f_P)
    heave = weight_body[2] if form == "lever_arm_and_heave_shift" else ca.SX(0.0)  # g_0(3) (otter.m 267-268)
    g0 = ca.vertcat(heave, moment[0], moment[1])  # g_0(3:5): heave, roll, pitch
    trim = ca.vertcat(0.0, 0.0, ca.solve(restoring_matrix[2:5, 2:5], g0), 0.0)  # G(3:5,3:5)^-1 g_0 (otter.m 271)
    kinematic = trim if form == "lever_arm_and_heave_shift" else ca.SX.zeros(6, 1)  # (otter.m 272, 275)
    return ca.Function("payload_trim_offset", [eta, restoring_matrix, payload_mass, payload_position, gravity],
                       [trim, kinematic],
                       ["eta", "restoring_matrix", "payload_mass", "payload_position", "gravity"],
                       ["trim_offset", "kinematic_attitude_offset"])
