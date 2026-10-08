"""Rigid-body + added-mass kinetics (CasADi): ``M_RB``, ``M_A``, ``C_RB(nu)``,
``C_A(nu_r)`` and the block function ``rigid_body_casadi``.

The block is one CasADi function. Its inputs are the state (``nu``,
``nu_r``, 6x1 each) and one named input per declared parameter of the
chosen form (``mass_properties.rigid_body_parameters``); its outputs are
named: ``M``, ``C_RB``, ``C_A``, ``M_RB``, ``M_A`` (6x6), ``mass``,
``center_of_gravity`` (CO -> CG), ``inertia`` (about the CG), and per form ``displaced_volume`` and
``wetted_surface`` (``"displacement_hull"``) or ``lamb_k_factors``
(``"spheroid"``). The form selectors (``mass_properties``, ``coriolis``,
``stabilize_added_mass_coriolis``) are keywords of the builder: they choose
which graph is built and are not inputs. Numbers are checked once where
they enter (``check_rigid_body_values``) and may be frozen into the function
with ``more_transformations.more_casadi_transformations.freeze``.

Body frame z down, CO is the origin, ``r_g`` points CO -> CG, ``I`` is the
inertia about the CG. ``S`` and ``H`` come from
``more_transformations.more_casadi_transformations``.

Equations (keys in References):

* ``H(r) = [[I3, S(r)^T], [0, I3]]`` (Fossen 2011, eq. 3.24, p. 49;
  MSS ``Hmtrx.m`` 16-18).
* ``M_RB = H(r_g)^T diag(m I3, I) H(r_g)`` (Fossen 2011, eqs. 3.26, 3.28,
  p. 50, and 3.44, p. 52; MSS ``rbody.m`` 36-44).
* Three ``C_RB`` parametrisations, one force ``C_RB(nu) nu``:

  - ``"co"``: ``C_RB = H^T diag(m S(w), -S(I w)) H`` (Fossen 2011, eqs. 3.21,
    p. 49, and 3.27, p. 50; MSS ``rbody.m`` 39-45, ``spheroid.m`` 46-52),
    depends on ``w = nu[3:6]`` only.
  - ``"book"``: ``C_RB = [[m S(w), -m S(w) S(r)], [m S(r) S(w), -S(I_O w)]]``,
    ``I_O = I - m S(r)^2`` (Fossen 2011, eqs. 3.29, p. 50, and 3.57, p. 55;
    parallel axes eq. 3.34, p. 50). ``co - book = m |r|^2 S(w)`` in the moment
    block only; the two give the same ``C_RB(nu) nu`` (Fossen 2011, eq. 3.30,
    p. 50).
  - ``"lagrangian"``: ``C = m2c(M_RB, nu) = [[0, -S(M11 v1 + M12 v2)],
    [-S(M11 v1 + M12 v2), -S(M21 v1 + M22 v2)]]`` (Fossen 2011, Theorem 3.2,
    eq. 3.46, p. 53, and eq. 3.55, p. 55; MSS ``m2c.m`` 35-48), depends on
    the linear velocity too.

* ``M = M_RB + M_A`` (Fossen 2011, eq. 6.48, p. 120), symmetric and positive
  definite (Fossen 2011, Property 3.1, eq. 3.43, p. 52; Property 6.1,
  p. 118; Sec. 7.5.2, p. 171).
* ``C_A(nu_r) = m2c(M_A, nu_r)`` with ``M_A`` symmetrised,
  ``M = (M + M^T)/2`` (Fossen 2011, eq. 6.43, p. 120, and eq. 6.40, p. 118;
  MSS ``m2c.m`` 33). For a diagonal ``M_A`` it is Fossen 2011 eq. 6.54,
  p. 121.
* Stabilised ``C_A``: the pitch-heave, pitch-surge, yaw-surge and yaw-sway
  entries set to zero in both triangles (MSS ``remus100.m`` 207-210). This is
  a modelling choice of that file, not physics: it removes the Munk moment
  of a slender body (Fossen 2011, eq. 6.52, p. 121), which MSS states in
  ``remus100.m`` 202-206.

The Lagrangian ``C_A`` symmetrises ``M_A`` as MSS ``m2c.m`` 33 does (no
difference for the symmetric ``M_A`` every form builds).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    eqs. 3.21-3.57, pp. 49-55; Ch. 6, eqs. 6.38-6.54, pp. 118-121; Ch. 7,
    Sec. 7.5.2, p. 171.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/kinematics/Hmtrx.m`` 16-18; ``LIBRARY/modeling/rbody.m`` 35-45,
    ``m2c.m`` 33-48, ``spheroid.m`` 45-52, ``imlay61.m`` 42;
    ``CRAFT/AUV/models/remus100.m`` 202-210.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import (
    MatrixTransforms,
    check_symmetric_positive_definite,
    check_values,
    symbols,
)

from more_dynamics.models.vehicles.hull_parts.added_mass.added_mass import STABILIZED_ADDED_MASS_CORIOLIS_PAIRS, added_mass_matrix
from .mass_properties import (
    COMMON_QUANTITIES,
    FORM_QUANTITIES,
    mass_properties_casadi,
    rigid_body_parameters,
)

CORIOLIS_FORMS = ("co", "book", "lagrangian")

# Named outputs of the block before the mass properties of the form.
MATRIX_OUTPUTS = ("M", "C_RB", "C_A", "M_RB", "M_A")


def rigid_body_mass_matrix(mass, inertia, r_g):
    """``M_RB = H(r_g)^T diag(m I3, I) H(r_g)`` (Fossen 2011, eqs. 3.26, 3.28,
    p. 50; MSS ``rbody.m`` 36-44).

    Equals ``[[m I3, -m S(r)], [m S(r), I - m S(r)^2]]`` (Fossen 2011,
    eq. 3.44, p. 52) to rounding; one function serves every form.
    """
    zeros = ca.SX.zeros(3, 3)
    # M_RB^CG = diag(m I3, I) (Fossen 2011, eq. 3.21, p. 49; rbody.m 36-37)
    mass_at_cg = ca.vertcat(
        ca.horzcat(mass * ca.SX.eye(3), zeros),
        ca.horzcat(zeros, inertia),
    )
    transform = MatrixTransforms.H_matrix(r_g)  # (Fossen 2011, eq. 3.24, p. 49; Hmtrx.m 16-18)
    return transform.T @ mass_at_cg @ transform  # (Fossen 2011, eq. 3.26, p. 50; rbody.m 43-44)


def lagrangian_coriolis_casadi(mass_matrix, nu):
    """``m2c(M, nu)`` 6-DOF: ``[[0, -S(M11 v1 + M12 v2)], [-S(M11 v1 + M12 v2),
    -S(M21 v1 + M22 v2)]]`` with ``M`` symmetrised and ``M21 = M12^T``
    (Fossen 2011, Theorem 3.2, eq. 3.46, p. 53; MSS ``m2c.m`` 33-48)."""
    matrix = MatrixTransforms.symmetrize(mass_matrix)  # M = (M + M^T)/2 (Fossen 2011, eq. 6.40, p. 118; m2c.m 33)
    v1, v2 = nu[0:3], nu[3:6]
    linear = matrix[0:3, 0:3] @ v1 + matrix[0:3, 3:6] @ v2  # (Fossen 2011, eq. 3.49, p. 54; m2c.m 44)
    angular = matrix[0:3, 3:6].T @ v1 + matrix[3:6, 3:6] @ v2  # (Fossen 2011, eq. 3.50, p. 54; m2c.m 39, 45)
    s_linear = MatrixTransforms.skew(linear)
    # (Fossen 2011, eq. 3.46, p. 53; m2c.m 47-48)
    return ca.vertcat(
        ca.horzcat(ca.SX.zeros(3, 3), -s_linear),
        ca.horzcat(-s_linear, -MatrixTransforms.skew(angular)),
    )


def rigid_body_coriolis_casadi(mass, inertia, r_g, nu, coriolis):
    """``C_RB(nu)`` in the chosen parametrisation (module docstring)."""
    omega = nu[3:6]
    s_omega = MatrixTransforms.skew(omega)
    if coriolis == "co":
        zeros = ca.SX.zeros(3, 3)
        # C_RB^CG = diag(m S(w), -S(I w)) (Fossen 2011, eq. 3.21, p. 49; rbody.m 39-40)
        coriolis_at_cg = ca.vertcat(
            ca.horzcat(mass * s_omega, zeros),
            ca.horzcat(zeros, -MatrixTransforms.skew(inertia @ omega)),
        )
        transform = MatrixTransforms.H_matrix(r_g)  # (Fossen 2011, eq. 3.24, p. 49; Hmtrx.m 16-18)
        return transform.T @ coriolis_at_cg @ transform  # (Fossen 2011, eq. 3.27, p. 50; rbody.m 43-45)
    if coriolis == "book":
        s_r = MatrixTransforms.skew(r_g)
        inertia_at_co = inertia - mass * s_r @ s_r  # I_O = I - m S(r)^2 (Fossen 2011, Theorem 3.1, eq. 3.34, p. 50)
        # (Fossen 2011, eqs. 3.29, p. 50, and 3.57, p. 55)
        return ca.vertcat(
            ca.horzcat(mass * s_omega, -mass * s_omega @ s_r),
            ca.horzcat(mass * s_r @ s_omega, -MatrixTransforms.skew(inertia_at_co @ omega)),
        )
    if coriolis == "lagrangian":
        # m2c(M_RB, nu) (Fossen 2011, eq. 3.55, p. 55; m2c.m 35-48)
        return lagrangian_coriolis_casadi(rigid_body_mass_matrix(mass, inertia, r_g), nu)
    raise ValueError(f"coriolis must be one of {CORIOLIS_FORMS}, got {coriolis!r}")


def added_mass_coriolis_casadi(added_mass, nu_r, stabilize):
    """``C_A(nu_r) = m2c(M_A, nu_r)`` (Fossen 2011, eq. 6.43, p. 120; MSS
    ``m2c.m`` 33-48); ``stabilize`` zeroes the pitch-heave, pitch-surge,
    yaw-surge and yaw-sway entries (MSS ``remus100.m`` 207-210).

    The stabilised form is a modelling choice, not physics: it removes the
    Munk moment of a slender body (Fossen 2011, eq. 6.52, p. 121), as MSS
    states in ``remus100.m`` 202-206.
    """
    coriolis = lagrangian_coriolis_casadi(added_mass, nu_r)  # (Fossen 2011, eq. 6.43, p. 120)
    if stabilize:  # (remus100.m 207-210)
        for row, column in STABILIZED_ADDED_MASS_CORIOLIS_PAIRS:
            coriolis[row, column] = 0.0  # (remus100.m 207-210)
            coriolis[column, row] = 0.0  # (remus100.m 207-210)
    return coriolis


def _check_selectors(coriolis, stabilize_added_mass_coriolis):
    if coriolis not in CORIOLIS_FORMS:
        raise ValueError(f"coriolis must be one of {CORIOLIS_FORMS}, got {coriolis!r}")
    if not isinstance(stabilize_added_mass_coriolis, bool):
        raise ValueError("stabilize_added_mass_coriolis must be True or False")


def rigid_body_outputs(mass_properties="hull_with_payload"):
    """The output names of ``rigid_body_casadi`` for one form, in order."""
    rigid_body_parameters(mass_properties)  # refuses an unknown form
    return MATRIX_OUTPUTS + COMMON_QUANTITIES + FORM_QUANTITIES[mass_properties]


def rigid_body_casadi(
    *,
    mass_properties="hull_with_payload",
    coriolis="co",
    stabilize_added_mass_coriolis=False,
):
    """The block: ``(nu, nu_r, <each declared parameter by name>) -> (M,
    C_RB(nu), C_A(nu_r), M_RB, M_A, mass, center_of_gravity, inertia,
    <quantities of the form>)``.

    Contract
    --------
    Keywords choose the graph: ``mass_properties`` one of
    ``MASS_PROPERTIES_FORMS``, ``coriolis`` one of ``CORIOLIS_FORMS``,
    ``stabilize_added_mass_coriolis`` a bool; an unknown value raises
    ``ValueError``. Inputs: ``nu`` and ``nu_r`` (6x1, body-fixed velocity and
    velocity relative to the water, SI) and the parameters of
    ``rigid_body_parameters(mass_properties)``, each under its own name and
    declared shape. Outputs: the names of ``rigid_body_outputs
    (mass_properties)``. Called with numbers it returns numbers; called with
    a symbol for a parameter it returns expressions in that symbol
    (identification). No value is checked inside the graph.
    """
    declared = rigid_body_parameters(mass_properties)
    _check_selectors(coriolis, stabilize_added_mass_coriolis)
    p = symbols(declared)
    nu = ca.SX.sym("nu", 6)
    nu_r = ca.SX.sym("nu_r", 6)

    quantities = mass_properties_casadi(mass_properties, p)
    mass, r_g, inertia = quantities["mass"], quantities["center_of_gravity"], quantities["inertia"]
    outputs = {
        "M_RB": rigid_body_mass_matrix(mass, inertia, r_g),  # (Fossen 2011, eq. 3.26, p. 50; rbody.m 43-44)
        "M_A": added_mass_matrix(quantities["added_mass_derivatives"]),  # (Fossen 2011, eq. 6.53, p. 121)
    }
    outputs["M"] = outputs["M_RB"] + outputs["M_A"]  # (Fossen 2011, eq. 6.48, p. 120)
    outputs["C_RB"] = rigid_body_coriolis_casadi(mass, inertia, r_g, nu, coriolis)
    outputs["C_A"] = added_mass_coriolis_casadi(outputs["M_A"], nu_r, stabilize_added_mass_coriolis)
    outputs.update(quantities)

    names = [d.name for d in declared]
    output_names = rigid_body_outputs(mass_properties)
    return ca.Function(
        "rigid_body",
        [nu, nu_r, *[p[name] for name in names]],
        [outputs[name] for name in output_names],
        ["nu", "nu_r", *names],
        list(output_names),
    )


def check_rigid_body_values(values, *, mass_properties="hull_with_payload"):
    """Numbers for one form, checked: ``{name: ca.DM}`` of the declared shape.

    Contract
    --------
    ``values`` maps every declared name of ``mass_properties`` to a number,
    a (nested) list or an array. Raises ``ValueError`` naming the parameter
    and its unit when a name is missing or not declared, a shape differs, an
    entry is not finite or lies outside its declared range
    (``check_values``); when a ``"spheroid"`` is not prolate
    (``semi_major_axis > semi_minor_axis`` is required, MSS ``imlay61.m``
    42); and when the total mass matrix ``M`` these numbers give is not
    symmetric positive definite (Fossen 2011, Property 3.1, eq. 3.43, p. 52;
    Sec. 7.5.2, p. 171). This is the check a plugin runs once before it
    freezes the numbers into the block.
    """
    declared = rigid_body_parameters(mass_properties)
    numbers = check_values(declared, values)
    if mass_properties == "spheroid":
        a, b = float(numbers["semi_major_axis"]), float(numbers["semi_minor_axis"])
        if not a > b:  # prolate spheroid, 0 < e < 1 (imlay61.m 42)
            raise ValueError(
                f"parameter 'semi_major_axis' [m] must be larger than 'semi_minor_axis' [m] "
                f"(prolate spheroid), got {a:g} and {b:g}"
            )
    block = rigid_body_casadi(mass_properties=mass_properties)
    out = block(nu=ca.DM.zeros(6), nu_r=ca.DM.zeros(6), **numbers)
    check_symmetric_positive_definite(out, "M")  # M = M^T > 0 (Fossen 2011, eq. 3.43, p. 52; Sec. 7.5.2, p. 171)
    return numbers
