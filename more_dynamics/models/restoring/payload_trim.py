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

    f_P  = R(eta)^T [0 0 m_P g]^T                      the payload weight in BODY      (otter_legacy.m 265)
    m_P' = r_P x f_P                                    its moment about the CO         (otter_legacy.m 266)
    g_0  = [f_P(3); m_P'(1); m_P'(2)]                   heave, roll, pitch              (otter_legacy.m 267-268)
    eta_0 = [0; 0; G(3:5, 3:5)^-1 g_0; 0]                                               (otter_legacy.m 271)

with the rotation at the actual attitude. The reference of ``eta`` is the **unloaded** waterline: ``z = 0`` is the
equilibrium of the bare hull, the displaced volume ``nabla = (m + m_P) / rho`` of the hydrostatic block only sizes the
draft and the metacentric heights (``G33 = rho g A_wp`` does not read it), and no buoyancy force is in the model, so
the payload weight is a wrench ``g_0`` on that reference, which the linear restoring has to balance. The wrench of a
weight at ``r_P`` about the CO equals the force at the centre of flotation plus the couple about it, so ``G eta_0 =
g_0`` is the textbook equilibrium of a floating body with a small added weight (Fossen 2011, §4.2 and §4.4): a
parallel sinkage ``m_P / (rho A_wp)`` of the centre of flotation and a rotation about it from the couple of the weight
about the centre of flotation. A payload above the centre of flotation sinks the craft level; the sinkage of the
centre of flotation, the roll and the pitch do not depend on where the body origin is (``tests/restoring/
test_payload_trim_origin.py``). Dropping the heave term, i.e. reading the payload weight as already carried by the
displacement, moves the reference to the loaded waterline and takes the moment about the origin instead of the centre
of flotation: the trim then depends on where the origin sits (PT1, E-211, GC2 section E; it supersedes the E-109
reading).

Forms (they differ in the kinematics, and the third in ``g_0``):

* ``"weight_wrench"`` (physics form, default): ``g_0`` as above; the kinematics are evaluated at the actual attitude
  (``kinematic_attitude_offset`` zero). Raise: the payload weighed, its position measured, and the trim read against a
  measured static draft, roll and pitch.
* ``"lever_arm_and_heave_shift"`` (MSS ``otter_legacy.m`` 264-275 at 6e78e98, a comparison form): the same ``g_0`` and
  ``eta_0``, with the kinematics evaluated at the shifted attitude ``eta - eta_0`` (``kinematic_attitude_offset =
  eta_0``; F-MSS-24, repaired upstream in 8667294), so a craft built on it equals that file with a payload.
* ``"lever_arm_moment"`` (legacy, **origin-dependent, not physical**): ``g_0 = [0; m_P'(1); m_P'(2)]``, the heave term
  dropped (E-109, withdrawn by E-211). Kept only to reproduce the frozen references of the retired model file whose
  payload trim was written that way (``grethe.m`` 412-419 at 679a6e8) and of the withdrawn ``otter-payload-trim``
  branch. Do not select it for a vehicle.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control, 1st ed. Wiley. Eq. 4.25,
    p. 65 (the linear restoring); §4.2 and §4.4, eqs. 4.94-4.100, pp. 74-76 (the equilibrium under an added weight).
[MSS] Fossen, T. I. MSS, MIT, @ 6e78e98: CRAFT/USV/models/otter_legacy.m 127-128 (loaded displacement), 264-275 (payload
    force and moment, trim condition, shifted attitude); LIBRARY/kinematics/Rzyx.m, Smtrx.m. Header of otter_legacy.m:
    z = 0 is the unloaded waterline (draft 13.4 cm for m_P = 0).

Author:    Enio Krizman
Date:      2026-10-26 (PT1 correction 2026-11-02)
"""
import casadi as ca

from more_transformations.more_casadi_transformations import MatrixTransforms

PAYLOAD_TRIM_FORMS = ("weight_wrench", "lever_arm_and_heave_shift", "lever_arm_moment")


def check_payload_trim_form(form):
    if form not in PAYLOAD_TRIM_FORMS:
        raise ValueError(f"payload_trim_form must be one of {PAYLOAD_TRIM_FORMS}, got {form!r}")


def payload_trim_offset(form) -> ca.Function:
    """``(eta, restoring_matrix, payload_mass, payload_position, gravity) -> (trim_offset, kinematic_attitude_offset)``
    for ``form`` (one of ``PAYLOAD_TRIM_FORMS``; the module docstring): the shift ``eta_0`` (6, zero in surge, sway and
    yaw) of the equilibrium of ``g = G (eta - eta_0)``, and the shift of the attitude the kinematics are evaluated
    at (``eta_0`` for ``"lever_arm_and_heave_shift"``, zero otherwise)."""
    check_payload_trim_form(form)
    eta = ca.SX.sym("eta", 6)
    restoring_matrix = ca.SX.sym("restoring_matrix", 6, 6)
    payload_mass = ca.SX.sym("payload_mass")
    payload_position = ca.SX.sym("payload_position", 3)
    gravity = ca.SX.sym("gravity")
    weight_body = MatrixTransforms.Rzyx(eta[3:6]).T @ ca.vertcat(0.0, 0.0, payload_mass * gravity)  # (otter_legacy.m 265)
    moment = ca.cross(payload_position, weight_body)  # (otter_legacy.m 266; Smtrx(r_P) f_P)
    heave = ca.SX(0.0) if form == "lever_arm_moment" else weight_body[2]  # g_0(3): the weight; the legacy form drops it
    g0 = ca.vertcat(heave, moment[0], moment[1])  # g_0(3:5): heave, roll, pitch
    trim = ca.vertcat(0.0, 0.0, ca.solve(restoring_matrix[2:5, 2:5], g0), 0.0)  # G(3:5,3:5)^-1 g_0 (otter_legacy.m 271)
    kinematic = trim if form == "lever_arm_and_heave_shift" else ca.SX.zeros(6, 1)  # (otter_legacy.m 272, 275)
    return ca.Function("payload_trim_offset", [eta, restoring_matrix, payload_mass, payload_position, gravity],
                       [trim, kinematic],
                       ["eta", "restoring_matrix", "payload_mass", "payload_position", "gravity"],
                       ["trim_offset", "kinematic_attitude_offset"])
