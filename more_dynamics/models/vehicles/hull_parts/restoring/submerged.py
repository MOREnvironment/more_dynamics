"""Submerged hydrostatics: restoring vector ``g(eta)`` (CasADi).

``submerged_hydrostatics_parameters()`` declares the four parameters (name,
shape, SI unit, meaning, admissible range); ``submerged_hydrostatics_casadi()``
builds the block ``(eta, weight, buoyancy, center_of_gravity,
center_of_buoyancy) -> g``. No numbers live in the module: a plugin freezes
them into the block once, identification calls it with symbols.

Equations (keys in References):

* ``W = m g``, ``B = rho g nabla`` (Fossen 2011, eq. 4.1, p. 59); a
  neutrally buoyant body has ``W = B`` (Fossen 2011, eq. 4.7, p. 61).
* ``g(eta) = -[f_g + f_b; r_bg x f_g + r_bb x f_b]`` with ``f_g = R^T [0, 0,
  W]``, ``f_b = -R^T [0, 0, B]`` (Fossen 2011, eq. 4.5, p. 60). With
  ``R(3,:) = [-s(theta), c(theta) s(phi), c(theta) c(phi)]`` (Fossen 2011,
  eq. 2.18, p. 22; MSS ``Rzyx.m`` 19) this is MSS ``gRvect.m`` 26-32:
  ``g[0:3] = -(W - B) R(3,:)^T``, ``g[3:6] = -m x R(3,:)^T`` with the first
  moment ``m = r_bg W - r_bb B``; expanded in Euler angles it is Fossen 2011,
  eq. 4.6, p. 60 (MSS ``gvect.m`` 22-28).

``g`` is the left-hand-side restoring vector: the equation of motion is
``M nu_dot + ... + g = tau`` (Fossen 2011, eq. 4.5, p. 60; MSS
``remus100.m`` 258, ``... - g``), so the force applied to the vehicle is
``-g``. For ``B > W`` at level attitude ``g_z = B - W > 0`` and the applied
``-g_z`` points up (negative z in NED).

``R(3,:)`` is ``MatrixTransforms.Rzyx_row3`` of
``more_transformations.more_casadi_transformations``. Conventions: NED, z
down; ``eta = [x y z phi theta psi]`` with zyx Euler angles;
``center_of_gravity`` (``r_bg``) and ``center_of_buoyancy`` (``r_bb``) are
measured from the CO in BODY. Yaw and position do not enter ``g``. One
vector serves the general, the ``gRvect.m`` and the neutrally buoyant
(``buoyancy = weight``) cases.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 2,
    eq. 2.18, p. 22; Ch. 4, eqs. 4.1-4.7, pp. 59-61.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579`` (these lines are the same at ``72656d1``):
    ``LIBRARY/modeling/gRvect.m`` 26-32, ``gvect.m`` 22-28;
    ``LIBRARY/kinematics/Rzyx.m`` 16-19; ``CRAFT/AUV/models/remus100.m``
    231, 258.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter, symbols

_DECLARED = (
    Parameter("weight", (1, 1), "N", "weight W = m g", 0.0, minimum_exclusive=True),
    Parameter("buoyancy", (1, 1), "N", "buoyancy B = rho g nabla", 0.0, minimum_exclusive=True),
    Parameter("center_of_gravity", (3, 1), "m", "CO -> CG r_bg, body axes (FRD)"),
    Parameter("center_of_buoyancy", (3, 1), "m", "CO -> CB r_bb, body axes (FRD)"),
)


def submerged_hydrostatics_parameters():
    """The declared parameter set: a tuple of ``Parameter`` (name, shape, SI
    unit, meaning, admissible range), in the order of the block's inputs."""
    return _DECLARED


def submerged_hydrostatics_casadi():
    """The block: ``(eta, weight, buoyancy, center_of_gravity,
    center_of_buoyancy) -> g`` (6x1, N and N m); Fossen 2011, eq. 4.5, p. 60;
    MSS ``gRvect.m`` 26-32.

    Contract
    --------
    Inputs: ``eta`` (6x1, ``[x y z phi theta psi]``, NED, m and rad) and the
    parameters of ``submerged_hydrostatics_parameters()``, each under its own
    name and declared shape. Output: ``g``. The force on the vehicle is
    ``-g``. Called with numbers it returns numbers; with a symbol for a
    parameter, expressions in it. No value is checked inside the graph
    (``check_values`` of the helper does that at the plugin boundary).
    """
    p = symbols(_DECLARED)
    eta = ca.SX.sym("eta", 6)
    phi, theta = eta[3], eta[4]
    weight, buoyancy = p["weight"], p["buoyancy"]

    # Third row of R = Rzyx(phi, theta, psi), body to NED (Fossen 2011, eq. 2.18, p. 22; Rzyx.m 19)
    row3 = MatrixTransforms.Rzyx_row3(phi, theta)
    r31, r32, r33 = row3[0], row3[1], row3[2]

    net = weight - buoyancy  # W - B (Fossen 2011, eq. 4.5, p. 60; gRvect.m 27-29)
    moment = p["center_of_gravity"] * weight - p["center_of_buoyancy"] * buoyancy  # r_bg W - r_bb B (gRvect.m 30-32)
    mx, my, mz = moment[0], moment[1], moment[2]

    # (Fossen 2011, eqs. 4.5-4.6, p. 60; gRvect.m 26-32, one line per row)
    g = ca.vertcat(
        -net * r31,  # (gRvect.m 27)
        -net * r32,  # (gRvect.m 28)
        -net * r33,  # (gRvect.m 29)
        -my * r33 + mz * r32,  # (gRvect.m 30)
        -mz * r31 + mx * r33,  # (gRvect.m 31)
        -mx * r32 + my * r31,  # (gRvect.m 32)
    )
    names = [d.name for d in _DECLARED]
    return ca.Function(
        "submerged_hydrostatics",
        [eta, *[p[name] for name in names]],
        [g],
        ["eta", *names],
        ["g"],
    )
