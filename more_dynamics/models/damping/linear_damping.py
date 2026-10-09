"""Linear damping: ``D`` and ``tau = -D nu_r`` (CasADi).

Three forms of one block, chosen by the selector ``form`` (keys in
References):

* ``"submerged"``: MSS ``Dmtrx.m``'s submerged branch with the surge fade of
  ``remus100.m``;
* ``"floating"``: ``Dmtrx.m``'s surface-craft branch, all six terms from time
  constants and damping ratios;
* ``"surface"``: the form of MSS ``otter.m``, surge from a thrust and a top
  speed, yaw with a quadratic term.

``linear_damping_parameters(form, smooth_speed=...)`` declares the numbers of
a form (name, shape, SI unit, meaning, admissible range);
``linear_damping_couplings(form)``
declares the matrices it takes from other blocks (couplings: the rigid
body's ``M_RB``, ``M_A`` or ``M``, the surface hydrostatics' ``G``), wired by
name inside a vehicle class, never typed as numbers;
``linear_damping_casadi(form=...)`` builds ``(nu_r, <couplings>,
<parameters>) -> (D, tau, <coefficients>)``. No numbers live in the module.

Selectors and parameters. The keywords ``form``, ``sway_damping_fade`` and
``smooth_speed`` decide the structure of the graph (which branch, which
terms); every number in an equation is a declared parameter and a named
input of the block, the regularisation ``smooth_speed_epsilon`` included.

Equations:

* ``D = -diag(Xu, Yv, Zw, Kp, Mq, Nr)`` positive, ``tau = -D nu_r``
  (Fossen 2011, eq. 6.62, p. 123; dissipative, Property 6.3, eqs.
  6.58-6.60, p. 123).
* Surge, sway, yaw from time constants ``D_ii = M_ii / T_i``; heave, roll,
  pitch from damping ratios ``D_ii = 2 zeta_i w_i M_ii`` (Fossen 2011,
  eqs. 6.71, p. 124, and 6.76-6.81, p. 125; MSS ``Dmtrx.m`` 30, 48-49,
  62-63).
* Natural frequencies: submerged ``w4 = sqrt(W (z_g - z_b) / M44)``,
  ``w5 = sqrt(W (z_g - z_b) / M55)``, ``T3 = T2`` (MSS ``Dmtrx.m`` 44-46);
  floating ``w_i = sqrt(G_ii / M_ii)``, i = 3, 4, 5 (Fossen 2011,
  eqs. 4.51-4.53, p. 68; MSS ``Dmtrx.m`` 58-60, ``otter.m`` 197-199).
* Submerged surge fade ``D11 = D11 exp(-3 U_r)``, ``U_r = sqrt(u_r^2 + v_r^2
  + w_r^2)`` (MSS ``remus100.m`` 127, 218).
* Surface: ``Xu = -thrust_max / U_max``, ``Yv = -M22 / T_sway``,
  ``Zw = -2 zeta3 w3 M33``, ``Kp = -2 zeta4 w4 M44``, ``Mq = -2 zeta5 w5 M55``,
  ``Nr = -M66 / T_yaw`` (MSS ``otter.m`` 202-207);
  ``tau = [Xu u_r, Yv v_r, Zw w_r, Kp p_r, Mq q_r, Nr (1 + k |r_r|) r_r]``
  (MSS ``otter.m`` 235-240).

Two sign conventions meet here. ``Dmtrx.m`` gives a positive ``D`` that is
subtracted. ``otter.m`` gives negative derivatives ``[Xu Yv Zw Kp Mq Nr]``
whose product with ``nu_r`` is added. Every form returns a positive ``D``
and the force ``tau`` on the vehicle, so the caller adds ``tau`` to the
right-hand side in each case.

Deviations from MSS:

* Inputs outside their domain are refused by ``check_linear_damping_values``
  naming the input, where MSS returns a complex, infinite or NaN number (a
  centre of gravity not below the centre of buoyancy, a negative restoring
  stiffness, a non-positive mass diagonal, time constant, weight or top
  speed). Negative ``damping_ratios`` are refused too (zero is allowed: an
  undamped mode); MSS does not reject them, and a negative ratio there gives
  a negative damping diagonal, energy fed into the mode against a positive
  restoring stiffness (Fossen 2011, Property 6.3, p. 123).
* MSS fixes ``zeta3 = 0.2`` inside the surface-craft branch (``Dmtrx.m`` 57);
  here it is a parameter.
* ``smooth_speed`` (selector of the submerged form, default ``False`` = MSS
  exactly). The speed ``U_r`` of ``remus100.m`` 127 has the derivative 0/0 at
  ``u_r = v_r = w_r = 0``, so the CasADi Jacobian of ``tau`` has non-finite
  entries there. ``smooth_speed=True`` replaces it by
  ``sqrt(u_r^2 + v_r^2 + w_r^2 + eps^2)``, which is smooth everywhere and
  departs from MSS by at most ``eps`` in the speed (a construction of this
  module, no published source). ``eps`` is then the declared parameter
  ``smooth_speed_epsilon`` (m/s, ``> 0``): a numerical regularisation chosen
  by the composition author, not a property of the vehicle.
* ``sway_damping_fade=True`` (submerged form) fades sway as surge, the
  setting of the earlier template reference generator; the default fades
  surge only, as ``remus100.m`` 218.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 4,
    eqs. 4.51-4.53, p. 68; Ch. 6, §6.4.1, eqs. 6.58-6.81, pp. 123-125.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``LIBRARY/modeling/Dmtrx.m`` 30-63;
    ``CRAFT/AUV/models/remus100.m`` 127, 217-218;
    ``CRAFT/USV/models/otter.m`` 197-207, 235-240.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, check_values, symbols

LINEAR_DAMPING_FORMS = ("submerged", "floating", "surface")

# Named outputs per form, in order.
LINEAR_DAMPING_OUTPUTS = {
    "submerged": ("D", "tau", "damping_coefficients"),
    "floating": ("D", "tau", "damping_coefficients"),
    "surface": ("D", "tau", "damping_derivatives"),
}


def _matrix(name, meaning):
    return Parameter(name, (6, 6), "kg, kg*m, kg*m^2 (6x6 blocks)", meaning)


_RIGID_BODY_MASS = _matrix("rigid_body_mass_matrix", "M_RB of the rigid-body block, about the CO")
_ADDED_MASS = _matrix("added_mass_matrix", "M_A of the rigid-body block")
_MASS = _matrix("mass_matrix", "M = M_RB + M_A of the rigid-body block")
_RESTORING = Parameter("restoring_matrix", (6, 6), "N/m, N, N*m/rad (6x6 blocks)",
                       "G of the hydrostatics block (only G33, G44, G55 used)")

_COUPLINGS = {
    "submerged": (_RIGID_BODY_MASS, _ADDED_MASS),
    "floating": (_RIGID_BODY_MASS, _ADDED_MASS, _RESTORING),
    "surface": (_MASS, _RESTORING),
}

_SMOOTH_SPEED_EPSILON = Parameter(
    "smooth_speed_epsilon", (1, 1), "m/s",
    "regularisation eps of sqrt(u_r^2 + v_r^2 + w_r^2 + eps^2) (smooth_speed=True); "
    "numerical, not a vehicle quantity",
    0.0, minimum_exclusive=True,
)

_DECLARATIONS = {
    "submerged": (
        Parameter("weight", (1, 1), "N", "weight W = m g", 0.0, minimum_exclusive=True),
        Parameter("center_of_gravity", (3, 1), "m", "CO -> CG r_bg, body axes (FRD)"),
        Parameter("center_of_buoyancy", (3, 1), "m", "CO -> CB r_bb, body axes (FRD); above the CG"),
        Parameter("time_constants", (3, 1), "s", "[T1, T2, T6] of surge, sway, yaw (T3 = T2)",
                  0.0, minimum_exclusive=True),
        Parameter("damping_ratios", (2, 1), "1", "[zeta4, zeta5] of roll and pitch", 0.0),
    ),
    "floating": (
        Parameter("time_constants", (3, 1), "s", "[T1, T2, T6] of surge, sway, yaw",
                  0.0, minimum_exclusive=True),
        Parameter("damping_ratios", (3, 1), "1", "[zeta3, zeta4, zeta5] of heave, roll, pitch", 0.0),
    ),
    "surface": (
        Parameter("max_forward_thrust", (1, 1), "N", "thrust at the top speed (Xu = -thrust / speed)", 0.0),
        Parameter("max_speed", (1, 1), "m/s", "top forward speed U_max", 0.0, minimum_exclusive=True),
        Parameter("time_constants", (2, 1), "s", "[T_sway, T_yaw]", 0.0, minimum_exclusive=True),
        Parameter("damping_ratios", (3, 1), "1", "[zeta3, zeta4, zeta5] of heave, roll, pitch", 0.0),
        Parameter("yaw_damping_nonlinearity", (1, 1), "s/rad", "k of Nr (1 + k |r|) r", 0.0),
    ),
}


def _check_form(form):
    if form not in LINEAR_DAMPING_FORMS:
        raise ValueError(f"form must be one of {LINEAR_DAMPING_FORMS}, got {form!r}")


def linear_damping_parameters(form, *, smooth_speed=False):
    """The declared parameters of one form: a tuple of ``Parameter`` (name,
    shape, SI unit, meaning, admissible range), in the order of the block's
    inputs after the couplings; with ``smooth_speed=True`` (submerged form)
    ``smooth_speed_epsilon`` is the last one."""
    _check_selectors(form, False, smooth_speed)
    return _DECLARATIONS[form] + ((_SMOOTH_SPEED_EPSILON,) if smooth_speed else ())


def linear_damping_couplings(form):
    """The inputs of one form that are outputs of other blocks (``M_RB``,
    ``M_A`` or ``M`` of the rigid body; ``G`` of the hydrostatics), declared
    as ``Parameter`` (name, shape, unit, meaning) in the order of the block's
    inputs after ``nu_r``."""
    _check_form(form)
    return _COUPLINGS[form]


def _relative_speed(nu_r, smooth_speed_epsilon):
    """``U_r = sqrt(u_r^2 + v_r^2 + w_r^2)`` (MSS ``remus100.m`` 127).

    With a ``smooth_speed_epsilon`` symbol (``smooth_speed=True``) its square
    is added under the root; ``None`` keeps the MSS line.
    """
    speed_squared = nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2  # (remus100.m 127)
    if smooth_speed_epsilon is not None:  # construction of this module, departure from remus100.m 127 (module docstring)
        speed_squared = speed_squared + smooth_speed_epsilon**2
    return ca.sqrt(speed_squared)


def _submerged(nu_r, c, p, sway_damping_fade, smooth_speed_epsilon):
    """MSS ``Dmtrx.m`` 30-49, submerged branch (Fossen 2011, eqs. 6.76-6.80,
    p. 125), and the surge fade of ``remus100.m`` 218."""
    mass = c["rigid_body_mass_matrix"] + c["added_mass_matrix"]  # M = MRB + MA (Dmtrx.m 30)
    t1, t2, t6 = p["time_constants"][0], p["time_constants"][1], p["time_constants"][2]
    zeta4, zeta5 = p["damping_ratios"][0], p["damping_ratios"][1]
    metacentric_height = p["center_of_gravity"][2] - p["center_of_buoyancy"][2]  # z_g - z_b (Dmtrx.m 45-46)
    t3 = t2  # (Dmtrx.m 44)
    w4 = ca.sqrt(p["weight"] * metacentric_height / mass[3, 3])  # (Dmtrx.m 45)
    w5 = ca.sqrt(p["weight"] * metacentric_height / mass[4, 4])  # (Dmtrx.m 46)
    # M_ii / T_i and 2 zeta_i w_i M_ii (Fossen 2011, eqs. 6.76-6.81, p. 125; Dmtrx.m 48-49)
    coefficients = ca.vertcat(
        mass[0, 0] / t1,
        mass[1, 1] / t2,
        mass[2, 2] / t3,
        mass[3, 3] * 2 * zeta4 * w4,
        mass[4, 4] * 2 * zeta5 * w5,
        mass[5, 5] / t6,
    )
    fade = ca.exp(-3 * _relative_speed(nu_r, smooth_speed_epsilon))  # (remus100.m 218)
    diagonal = [coefficients[i] for i in range(6)]
    diagonal[0] = diagonal[0] * fade  # (remus100.m 218)
    if sway_damping_fade:  # the template setting: sway faded as surge (module docstring)
        diagonal[1] = diagonal[1] * fade
    damping = ca.diag(ca.vertcat(*diagonal))
    return {"D": damping, "tau": -damping @ nu_r, "damping_coefficients": coefficients}  # tau = -D nu_r (Fossen 2011, eq. 6.62, p. 123; remus100.m 258)


def _floating(nu_r, c, p):
    """MSS ``Dmtrx.m`` 51-63, the surface-craft branch (Fossen 2011,
    eqs. 6.76-6.81, p. 125). MSS fixes ``zeta3`` inside the function
    (``Dmtrx.m`` 57); here it is a parameter."""
    mass = c["rigid_body_mass_matrix"] + c["added_mass_matrix"]  # M = MRB + MA (Dmtrx.m 30)
    restoring = c["restoring_matrix"]
    t1, t2, t6 = p["time_constants"][0], p["time_constants"][1], p["time_constants"][2]
    zeta3, zeta4, zeta5 = p["damping_ratios"][0], p["damping_ratios"][1], p["damping_ratios"][2]
    # w_i = sqrt(G_ii / M_ii) (Fossen 2011, eqs. 4.51-4.53, p. 68; Dmtrx.m 58-60)
    w3 = ca.sqrt(restoring[2, 2] / mass[2, 2])  # (Dmtrx.m 58)
    w4 = ca.sqrt(restoring[3, 3] / mass[3, 3])  # (Dmtrx.m 59)
    w5 = ca.sqrt(restoring[4, 4] / mass[4, 4])  # (Dmtrx.m 60)
    # (Fossen 2011, eqs. 6.76-6.81, p. 125; Dmtrx.m 62-63)
    coefficients = ca.vertcat(
        mass[0, 0] / t1,
        mass[1, 1] / t2,
        mass[2, 2] * 2 * zeta3 * w3,
        mass[3, 3] * 2 * zeta4 * w4,
        mass[4, 4] * 2 * zeta5 * w5,
        mass[5, 5] / t6,
    )
    damping = ca.diag(coefficients)
    return {"D": damping, "tau": -damping @ nu_r, "damping_coefficients": coefficients}  # (Fossen 2011, eq. 6.62, p. 123)


def _surface(nu_r, c, p):
    """MSS ``otter.m`` 197-207 and 235-240: ``D = -diag([Xu Yv Zw Kp Mq
    Nr])`` (positive) and ``tau`` = the derivatives times ``nu_r``, yaw times
    ``1 + k |r|``; ``restoring_matrix`` about the centre of flotation."""
    mass, restoring = c["mass_matrix"], c["restoring_matrix"]
    t_sway, t_yaw = p["time_constants"][0], p["time_constants"][1]
    zeta3, zeta4, zeta5 = p["damping_ratios"][0], p["damping_ratios"][1], p["damping_ratios"][2]
    # (Fossen 2011, eqs. 4.51-4.53, p. 68; otter.m 197-199)
    w3 = ca.sqrt(restoring[2, 2] / mass[2, 2])  # (otter.m 197)
    w4 = ca.sqrt(restoring[3, 3] / mass[3, 3])  # (otter.m 198)
    w5 = ca.sqrt(restoring[4, 4] / mass[4, 4])  # (otter.m 199)
    derivatives = ca.vertcat(
        -p["max_forward_thrust"] / p["max_speed"],  # Xu (otter.m 202)
        -mass[1, 1] / t_sway,  # Yv (Fossen 2011, eq. 6.77, p. 125; otter.m 203)
        -2 * zeta3 * w3 * mass[2, 2],  # Zw (Fossen 2011, eq. 6.78, p. 125; otter.m 204)
        -2 * zeta4 * w4 * mass[3, 3],  # Kp (Fossen 2011, eq. 6.79, p. 125; otter.m 205)
        -2 * zeta5 * w5 * mass[4, 4],  # Mq (Fossen 2011, eq. 6.80, p. 125; otter.m 206)
        -mass[5, 5] / t_yaw,  # Nr (Fossen 2011, eq. 6.81, p. 125; otter.m 207)
    )
    forces = [derivatives[i] * nu_r[i] for i in range(5)]  # (otter.m 235-239)
    forces.append(  # (otter.m 240)
        derivatives[5] * (1 + p["yaw_damping_nonlinearity"] * ca.fabs(nu_r[5])) * nu_r[5]
    )
    damping = ca.diag(-derivatives)  # D = -diag(derivatives) (Fossen 2011, eq. 6.62, p. 123)
    return {"D": damping, "tau": ca.vertcat(*forces), "damping_derivatives": derivatives}


def _check_selectors(form, sway_damping_fade, smooth_speed):
    _check_form(form)
    if not isinstance(sway_damping_fade, bool):
        raise ValueError("sway_damping_fade must be True or False")
    if not isinstance(smooth_speed, bool):
        raise ValueError("smooth_speed must be True or False")
    if form != "submerged" and (sway_damping_fade or smooth_speed):
        raise ValueError("sway_damping_fade and smooth_speed belong to the submerged form only")


def linear_damping_casadi(*, form, sway_damping_fade=False, smooth_speed=False):
    """The block: ``(nu_r, <couplings>, <parameters>) -> (D, tau, <coefficients>)``.

    Contract
    --------
    Keywords choose the graph (selectors): ``form`` one of
    ``LINEAR_DAMPING_FORMS`` (required); for ``"submerged"`` only,
    ``sway_damping_fade`` (bool) and ``smooth_speed`` (bool; module
    docstring). An unknown value raises ``ValueError``. Inputs: ``nu_r`` (6x1,
    velocity relative to the water, body axes, SI), the couplings of
    ``linear_damping_couplings(form)`` and the parameters of
    ``linear_damping_parameters(form, smooth_speed=smooth_speed)`` (with
    ``smooth_speed=True``: ``smooth_speed_epsilon``), each under its own name
    and declared shape. Outputs: the names of
    ``LINEAR_DAMPING_OUTPUTS[form]`` — ``D`` (6x6, positive), ``tau`` (6x1,
    the force on the vehicle) and ``damping_coefficients`` (the ``Dmtrx.m``
    diagonal before any fade) or ``damping_derivatives`` (``[Xu Yv Zw Kp Mq
    Nr]``, negative). No value is checked inside the graph
    (``check_linear_damping_values`` does that).
    """
    _check_selectors(form, sway_damping_fade, smooth_speed)
    couplings, declared = _COUPLINGS[form], linear_damping_parameters(form, smooth_speed=smooth_speed)
    c, p = symbols(couplings), symbols(declared)
    nu_r = ca.SX.sym("nu_r", 6)
    if form == "submerged":
        outputs = _submerged(nu_r, c, p, sway_damping_fade, p.get("smooth_speed_epsilon"))
    elif form == "floating":
        outputs = _floating(nu_r, c, p)
    else:
        outputs = _surface(nu_r, c, p)
    names = [d.name for d in couplings] + [d.name for d in declared]
    symbols_in = {**c, **p}
    return ca.Function(
        f"{form}_linear_damping",
        [nu_r, *[symbols_in[name] for name in names]],
        [outputs[name] for name in LINEAR_DAMPING_OUTPUTS[form]],
        ["nu_r", *names],
        list(LINEAR_DAMPING_OUTPUTS[form]),
    )


def check_linear_damping_values(values, *, form, smooth_speed=False):
    """Numbers for one form, couplings included, checked: ``{name: ca.DM}``.

    Contract
    --------
    ``values`` maps every coupling and parameter name of ``form`` (and
    ``smooth_speed``) to a number, a (nested) list or an array. Raises ``ValueError`` naming the
    input when a name is missing or unknown, a shape differs, an entry is
    not finite or outside its range (``check_values``); when the diagonal of
    the mass matrix is not positive; when ``G33``, ``G44`` or ``G55`` of
    ``restoring_matrix`` is negative (``sqrt`` of ``Dmtrx.m`` 58-60,
    ``otter.m`` 197-199); and, submerged, when the centre of gravity is not
    below the centre of buoyancy (``sqrt`` of ``Dmtrx.m`` 45-46, z down).
    """
    declared = linear_damping_parameters(form, smooth_speed=smooth_speed)
    numbers = check_values(_COUPLINGS[form] + declared, values)
    if form == "surface":
        mass, label = numbers["mass_matrix"], "mass_matrix"
    else:
        mass = numbers["rigid_body_mass_matrix"] + numbers["added_mass_matrix"]  # M = MRB + MA (Dmtrx.m 30)
        label = "rigid_body_mass_matrix + added_mass_matrix"
    if not all(float(mass[i, i]) > 0.0 for i in range(6)):
        raise ValueError(f"the diagonal of {label} must be positive")
    if "restoring_matrix" in numbers:
        restoring = numbers["restoring_matrix"]
        if not all(float(restoring[i, i]) >= 0.0 for i in (2, 3, 4)):
            raise ValueError("G33, G44 and G55 of restoring_matrix must be non-negative")
    if form == "submerged":
        z_g, z_b = float(numbers["center_of_gravity"][2]), float(numbers["center_of_buoyancy"][2])
        if not z_g - z_b > 0.0:
            raise ValueError(
                "center_of_gravity must lie below center_of_buoyancy (z down), "
                f"got z_g = {z_g:g} m, z_b = {z_b:g} m"
            )
    return numbers
