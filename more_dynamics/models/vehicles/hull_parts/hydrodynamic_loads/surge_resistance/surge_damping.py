"""Surge damping: linear at low speed, quadratic at high speed, blended (CasADi).

Three forms of one block, chosen by the selector ``surge_form``.
``surge_damping_parameters(surge_form, ittc_reynolds_bound=...)`` declares
the numbers of a form (name, shape, SI unit, meaning, admissible range);
``surge_damping_couplings(surge_form)`` declares what it takes from the
rigid-body block (``mass``, ``wetted_surface`` or ``M``), wired by name inside
a vehicle class; ``surge_damping_casadi(surge_form=...)`` builds ``(nu_r,
<couplings>, <parameters>) -> (tau, <coefficients>)``. The keywords
``surge_form``, ``surge_blend`` and ``ittc_reynolds_bound`` are selectors
(they decide the structure of the graph); every number in an equation that a
composition may set is a declared parameter and a named input, the Reynolds
floor included.

Equations (keys in References)::

    A11   = c rho (m / rho)^(5/3) / L^2        (addedMassSurge.m 32-33; c is the
                                                declared parameter
                                                surge_added_mass_factor)
    Xudot = -A11                               (forceSurgeDamping.m 60 @ ac77394)
    Xu    = -(m - Xudot) / T1 = -(m + A11) / T1
            (Fossen 2011, eqs. 6.71, p. 124, and 6.76, p. 125;
             forceSurgeDamping.m 61 @ ac77394)
    X     = sigma Xu u_r + (1 - sigma) Xuu |u_r| u_r
                                               (forceSurgeDamping.m 79-82 @ ac77394)

MSS attributes the surge added mass to Söding (1982), not read here.

``surge_form`` (the calibration of the quadratic coefficient):

* ``"max_thrust"``: ``Xuu = -max_thrust / max_speed^2`` (the steady state at
  top speed, ``forceSurgeDamping.m`` 64-66 at ``ac77394``);
* ``"ittc"``: ``Xuu = -1/2 rho S (1 + k) C_f`` with the ITTC-1957 line
  ``C_f = 0.075 / (log10 Rn - 2)^2``, ``Rn = L |u_r| / nu``, ``nu`` the
  declared parameter ``kinematic_viscosity`` (Fossen 2011, eqs. 6.82-6.85,
  p. 125, with ``C_R = 0``; ``forceSurgeDamping.m`` 69-74 at ``ac77394``;
  ``XuuITTC.m`` 32-39), evaluated at the current speed;
* ``"exp_ittc"``: the ship law of MSS release 2.0.2 (``osv.m`` 182-185,
  ``hydroVessel.m`` 88-89 at ``cc07579``), ``D_nl,11 = exp(-k_u |u_r|) D11 -
  Xuu |u_r|``, ``X = -D_nl,11 u_r`` (``osv.m`` 204), with ``D11 = M(1,1) /
  T1`` (``Dmtrx.m`` 62; Fossen 2011, eq. 6.76, p. 125) and ``Xuu`` the ITTC
  line above. The linear term fades exponentially instead of through the
  ``tanh`` blend; ``k_u`` (``exponential_decay_rate``, s/m) is a parameter
  (MSS fixes 3, ``osv.m`` 183). ``XuuITTC.m`` 38 estimates the wetted
  surface with the Mumford formula from ``L, B, T, C_B``; here the wetted
  surface is an input, as in ``"ittc"``.

MSS hard-codes the form factor ``k = 0.1`` (``forceSurgeDamping.m`` 70,
``XuuITTC.m`` 33) and the crossover speed ``u_cross = 2`` m/s
(``forceSurgeDamping.m`` 57); here both are parameters (Fossen 2011, p. 125:
``k`` typically 0.1 in transit, 0.25 in DP).

Declared parameters, not fixed constants (rule 16, E-65): the surge
added-mass coefficient ``c`` (``surge_added_mass_factor``, dimensionless,
``"ittc"`` and ``"max_thrust"``) of ``A11 = c rho (m/rho)^(5/3) / L^2``
(``addedMassSurge.m`` 33; MSS attributes it to Söding 1982, not read here),
and the kinematic viscosity ``nu`` (``kinematic_viscosity``, m^2/s, the ITTC
forms) of ``Rn = L |u_r| / nu``. ``nu`` is a property of the water, not a
fixed number: Fossen 2011 (p. 125, below eq. 6.85, read) gives ``nu = 1e-6``
m^2/s as its value *at 20 degC* (``XuuITTC.m`` 32; ``forceSurgeDamping.m`` 69
at ``ac77394`` fix the same number); a composition raises it from a
measurement of the water on the day of the run (e.g. a CTD cast) in colder or
more saline water.

Fixed constants of the cited lines, still not parameters: the ITTC line's
0.075 and 2 (Fossen 2011, eq. 6.83, p. 125). These belong to the published
regression itself, not to a measurable property of the water or the vehicle.

Revision: ``LIBRARY/modeling/forceSurgeDamping.m`` is cited at MSS
``ac77394``, the last revision that holds it; it was deleted in MSS
``108ceda`` (release 2.0, 2026-10-06). ``addedMassSurge.m``, ``XuuITTC.m``,
``Dmtrx.m``, ``osv.m`` and ``hydroVessel.m`` are cited at ``cc07579``
(``osv.m`` 188 at ``72656d1``, the line that held the surge term at zero).

Deviations from ``forceSurgeDamping.m`` (each with a selector that gives
MSS's line back):

* ``Xudot = -A11`` as MSS (the linearised mode ``(m + A11) du/dt = Xu u``
  decays with the time constant ``T1``, Fossen 2011, eqs. 6.71 and 6.76).
* ``surge_blend="symmetric"`` (default): ``sigma = 1 - tanh(|u_r| / u_cross)``.
  ``"mss_tanh"``: ``sigma = 1 - tanh(u_r / u_cross)`` (line 79), which exceeds
  1 for ``u_r < 0``, so the quadratic term changes sign and pushes the vehicle
  astern, against the dissipative property of damping (Fossen 2011,
  Property 6.3, p. 123); the file's own header (lines 8-12) describes a
  symmetric blend. The two are equal for ``u_r >= 0``. Not used by
  ``"exp_ittc"`` (no blend).
* ``ittc_reynolds_bound="floor"`` (default, the ITTC forms):
  ``Rn = max(L |u_r| / nu, Rn_min)``, the line of MSS's own ``XuuITTC.m``
  34-36 (Fossen 2011, p. 125: a minimum ``Rn`` should be used, ``C_F`` blows
  up at low speed). ``Rn_min`` is then the declared parameter
  ``ittc_reynolds_floor`` (dimensionless, ``> 100``, the pole of the ITTC
  line; ``XuuITTC.m`` 34 sets ``1e5``): a numerical bound chosen by the
  composition author, not a property of the vehicle.
  ``"mss_offset"``: the unbounded ``log10(Rn + 1e-10)`` of
  ``forceSurgeDamping.m`` 71-73, which has a pole at ``Rn = 100``; no floor
  is declared. The two are equal for ``Rn >= Rn_min`` up to the ``1e-10``.
* Inputs outside their domain are refused by the declared ranges
  (``check_surge_damping_values``), naming the input.

Conventions: ``nu_r = [u v w p q r]`` relative to the water; only ``u_r``
enters. ``tau`` is the force on the vehicle (added to the right-hand side).
A vehicle that uses this block sets the surge term of its linear damping
matrix to zero (MSS ``osv.m`` 188 at ``72656d1``; ``osv.m`` 185 replaces it
at ``cc07579``), or surge damping is counted twice.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 6,
    Property 6.3, p. 123; eqs. 6.71-6.85, pp. 124-125.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``LIBRARY/modeling/addedMassSurge.m``
    32-33, ``XuuITTC.m`` 32-39, ``Dmtrx.m`` 62; ``CRAFT/SHIP/models/osv.m``
    182-185, 204; ``CRAFT/hydroVessel.m`` 88-89; revision ``72656d1``:
    ``CRAFT/SHIP/models/osv.m`` 188.
[MSS ac77394] Fossen, T. I. *Marine Systems Simulator (MSS)*,
    https://github.com/cybergalactic/MSS, MIT licence, revision ``ac77394``
    (2026-10-05): ``LIBRARY/modeling/forceSurgeDamping.m`` 57-82 (revision
    of 2025-09-23; the file is deleted in later revisions).

Author:    Enio Krizman
Date:      2026-10-07
"""

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, check_values, symbols

# C_f = 0.075 / (log10 Rn - 2)^2 (Fossen 2011, eq. 6.83, p. 125;
# forceSurgeDamping.m 73 at ac77394; XuuITTC.m 36).
ITTC_FRICTION_FACTOR = 0.075
ITTC_LOG10_OFFSET = 2.0

# forceSurgeDamping.m 71 (ac77394): added to Rn inside the logarithm (MSS line only).
MSS_REYNOLDS_OFFSET = 1e-10

SURGE_FORMS = ("ittc", "max_thrust", "exp_ittc")
SURGE_BLENDS = ("symmetric", "mss_tanh")
ITTC_REYNOLDS_BOUNDS = ("floor", "mss_offset")
_ITTC_FORMS = ("ittc", "exp_ittc")

# Named outputs per form, in order.
SURGE_DAMPING_OUTPUTS = {
    "ittc": ("tau", "added_mass", "linear_coefficient", "quadratic_coefficient"),
    "max_thrust": ("tau", "added_mass", "linear_coefficient", "quadratic_coefficient"),
    "exp_ittc": ("tau", "linear_coefficient", "quadratic_coefficient", "surge_damping_coefficient"),
}


def _positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


_MASS = _positive("mass", "kg", "rigid-body mass m (the rigid-body block's mass)")
_WETTED_SURFACE = _positive("wetted_surface", "m^2", "wetted surface S (the rigid-body block's wetted_surface)")
_MASS_MATRIX = Parameter("mass_matrix", (6, 6), "kg, kg*m, kg*m^2 (6x6 blocks)",
                         "M = M_RB + M_A of the rigid-body block (only M(1,1) used)")
_LENGTH = _positive("length", "m", "length L (added mass, Reynolds number)")
_WATER_DENSITY = _positive("water_density", "kg/m^3", "water density rho")
_TIME_CONSTANT = _positive("time_constant", "s", "surge time constant T1")
_CROSSOVER_SPEED = _positive("crossover_speed", "m/s", "u_cross of the tanh blend")
_FORM_FACTOR = Parameter("form_factor", (1, 1), "1", "hull form factor k of (1 + k) C_f", 0.0)
_KINEMATIC_VISCOSITY = _positive(
    "kinematic_viscosity", "m^2/s",
    "kinematic viscosity nu of the water, Rn = L |u_r| / nu (Fossen 2011, eq. 6.85, "
    "p. 125); a water property, not a fixed number — nu = 1e-6 at 20 degC (Fossen "
    "2011, p. 125, below eq. 6.85; XuuITTC.m 32; forceSurgeDamping.m 69 @ ac77394); "
    "raised from that reference by a measurement of the water on the day of the run",
)
_SURGE_ADDED_MASS_FACTOR = _positive(
    "surge_added_mass_factor", "1",
    "coefficient c of A11 = c rho (m/rho)^(5/3) / L^2 (addedMassSurge.m 33; MSS "
    "attributes it to Soding 1982, not read here); an empirical constant of the "
    "published formula, overridable with an identified or published value",
)
# The bound is strict: log10 Rn - 2 = 0 at Rn = 100 (Fossen 2011, eq. 6.83, p. 125).
_REYNOLDS_FLOOR = Parameter(
    "ittc_reynolds_floor", (1, 1), "1",
    "minimum Reynolds number Rn_min of the ITTC line (ittc_reynolds_bound='floor'); "
    "numerical bound, not a vehicle quantity",
    10.0**ITTC_LOG10_OFFSET, minimum_exclusive=True,
)

_COUPLINGS = {
    "ittc": (_MASS, _WETTED_SURFACE),
    "max_thrust": (_MASS,),
    "exp_ittc": (_MASS_MATRIX, _WETTED_SURFACE),
}

_DECLARATIONS = {
    "ittc": (_LENGTH, _WATER_DENSITY, _TIME_CONSTANT, _FORM_FACTOR, _CROSSOVER_SPEED,
             _KINEMATIC_VISCOSITY, _SURGE_ADDED_MASS_FACTOR),
    "max_thrust": (
        _LENGTH,
        _WATER_DENSITY,
        _TIME_CONSTANT,
        _positive("max_speed", "m/s", "top speed u_max"),
        Parameter("max_thrust", (1, 1), "N", "surge thrust at the top speed", 0.0),
        _CROSSOVER_SPEED,
        _SURGE_ADDED_MASS_FACTOR,
    ),
    "exp_ittc": (
        _LENGTH,
        _WATER_DENSITY,
        _TIME_CONSTANT,
        _FORM_FACTOR,
        Parameter("exponential_decay_rate", (1, 1), "s/m", "k_u of exp(-k_u |u_r|)", 0.0),
        _KINEMATIC_VISCOSITY,
    ),
}


def _check_form(surge_form):
    if surge_form not in SURGE_FORMS:
        raise ValueError(f"surge_form must be one of {SURGE_FORMS}, got {surge_form!r}")


def surge_damping_parameters(surge_form, *, ittc_reynolds_bound="floor"):
    """The declared parameters of one form: a tuple of ``Parameter`` (name,
    shape, SI unit, meaning, admissible range), in the order of the block's
    inputs after the couplings; ``kinematic_viscosity`` (the ITTC forms,
    ``"ittc"`` and ``"exp_ittc"``) and ``surge_added_mass_factor`` (``"ittc"``
    and ``"max_thrust"``) are declared parameters too (rule 16); on the ITTC
    forms with ``ittc_reynolds_bound="floor"`` ``ittc_reynolds_floor`` is the
    last one."""
    _check_bound(surge_form, ittc_reynolds_bound)
    floor = surge_form in _ITTC_FORMS and ittc_reynolds_bound == "floor"
    return _DECLARATIONS[surge_form] + ((_REYNOLDS_FLOOR,) if floor else ())


def surge_damping_couplings(surge_form):
    """The inputs of one form that are outputs of the rigid-body block
    (``mass``, ``wetted_surface``; ``M`` for ``"exp_ittc"``), declared as
    ``Parameter`` in the order of the block's inputs after ``nu_r``."""
    _check_form(surge_form)
    return _COUPLINGS[surge_form]


def _check_bound(surge_form, ittc_reynolds_bound):
    _check_form(surge_form)
    if ittc_reynolds_bound not in ITTC_REYNOLDS_BOUNDS:
        raise ValueError(
            f"ittc_reynolds_bound must be one of {ITTC_REYNOLDS_BOUNDS}, got {ittc_reynolds_bound!r}"
        )
    if surge_form not in _ITTC_FORMS and ittc_reynolds_bound != "floor":
        raise ValueError(f"ittc_reynolds_bound does not apply to surge_form={surge_form!r} (no ITTC line)")


def _check_selectors(surge_form, surge_blend, ittc_reynolds_bound):
    _check_bound(surge_form, ittc_reynolds_bound)
    if surge_blend not in SURGE_BLENDS:
        raise ValueError(f"surge_blend must be one of {SURGE_BLENDS}, got {surge_blend!r}")
    if surge_form == "exp_ittc" and surge_blend != "symmetric":
        raise ValueError("surge_blend does not apply to surge_form='exp_ittc' (no tanh blend)")


def _ittc_quadratic_coefficient(p, c, u_r):
    """``Xuu(u_r)``: ``XuuITTC.m`` 35-39 with the floor (``p`` holds
    ``ittc_reynolds_floor``), else ``forceSurgeDamping.m`` 72-74
    (``ac77394``); Fossen 2011, eqs. 6.82-6.85, p. 125."""
    reynolds = p["length"] / p["kinematic_viscosity"] * ca.fabs(u_r)  # Rn = L |u_r| / nu (Fossen 2011, eq. 6.85, p. 125; forceSurgeDamping.m 72 @ ac77394)
    if "ittc_reynolds_floor" in p:
        reynolds = ca.fmax(reynolds, p["ittc_reynolds_floor"])  # Rn = max(Rn, Re_min) (XuuITTC.m 35)
    else:
        reynolds = reynolds + MSS_REYNOLDS_OFFSET  # (forceSurgeDamping.m 73 @ ac77394)
    friction = ITTC_FRICTION_FACTOR / (ca.log10(reynolds) - ITTC_LOG10_OFFSET) ** 2  # (Fossen 2011, eq. 6.83, p. 125; XuuITTC.m 36)
    # -1/2 rho S (1 + k) C_f (Fossen 2011, eq. 6.82, p. 125; XuuITTC.m 39; forceSurgeDamping.m 74 @ ac77394)
    return -0.5 * p["water_density"] * c["wetted_surface"] * (1 + p["form_factor"]) * friction


def _blended(nu_r, c, p, surge_form, surge_blend):
    """``"ittc"`` and ``"max_thrust"``: ``forceSurgeDamping.m`` 57-82 at ``ac77394``."""
    u_r = nu_r[0]
    mass, length, rho = c["mass"], p["length"], p["water_density"]
    displaced_volume = mass / rho  # (addedMassSurge.m 32)
    added_mass = p["surge_added_mass_factor"] * rho * displaced_volume ** (5 / 3) / length**2  # (addedMassSurge.m 33)
    surge_acceleration_derivative = -added_mass  # (forceSurgeDamping.m 60 @ ac77394)
    # (Fossen 2011, eqs. 6.71, p. 124, and 6.76, p. 125; forceSurgeDamping.m 61 @ ac77394)
    linear_coefficient = -(mass - surge_acceleration_derivative) / p["time_constant"]
    if surge_form == "ittc":
        quadratic_coefficient = _ittc_quadratic_coefficient(p, c, u_r)
    else:
        quadratic_coefficient = -p["max_thrust"] / p["max_speed"] ** 2  # (forceSurgeDamping.m 66 @ ac77394)
    # symmetric |u_r| by default; u_r as forceSurgeDamping.m 79 (ac77394) with "mss_tanh"
    blend_speed = ca.fabs(u_r) if surge_blend == "symmetric" else u_r
    sigma = 1 - ca.tanh(blend_speed / p["crossover_speed"])  # (forceSurgeDamping.m 79 @ ac77394; |u_r| by default)
    surge = (  # (forceSurgeDamping.m 82 @ ac77394)
        sigma * linear_coefficient * u_r
        + (1 - sigma) * quadratic_coefficient * ca.fabs(u_r) * u_r
    )
    return {
        "tau": ca.vertcat(surge, 0.0, 0.0, 0.0, 0.0, 0.0),
        "added_mass": added_mass,
        "linear_coefficient": linear_coefficient,
        "quadratic_coefficient": quadratic_coefficient,
    }


def _exponential(nu_r, c, p):
    """``"exp_ittc"``: MSS ``osv.m`` 182-185, 204 (``hydroVessel.m`` 88-89)."""
    u_r = nu_r[0]
    linear_damping = c["mass_matrix"][0, 0] / p["time_constant"]  # D11 = M(1,1) / T1 (Fossen 2011, eq. 6.76, p. 125; Dmtrx.m 62)
    quadratic_coefficient = _ittc_quadratic_coefficient(p, c, u_r)  # (XuuITTC.m 35-39; osv.m 184)
    # D_nl,11 = exp(-k_u |u_r|) D11 - Xuu |u_r| (osv.m 185; hydroVessel.m 89)
    damping = ca.exp(-p["exponential_decay_rate"] * ca.fabs(u_r)) * linear_damping - quadratic_coefficient * ca.fabs(u_r)
    return {
        "tau": ca.vertcat(-damping * u_r, 0.0, 0.0, 0.0, 0.0, 0.0),  # X = -D_nl,11 u_r (osv.m 204)
        "linear_coefficient": -linear_damping,  # Xu = -D11 (Fossen 2011, eq. 6.62, p. 123)
        "quadratic_coefficient": quadratic_coefficient,
        "surge_damping_coefficient": damping,
    }


def surge_damping_casadi(*, surge_form, surge_blend="symmetric", ittc_reynolds_bound="floor"):
    """The block: ``(nu_r, <couplings>, <parameters>) -> (tau = [X 0 0 0 0 0],
    <coefficients>)``.

    Contract
    --------
    Keywords choose the graph (selectors): ``surge_form`` in ``SURGE_FORMS``
    (required); ``surge_blend`` in ``SURGE_BLENDS`` (``"ittc"``,
    ``"max_thrust"`` only); ``ittc_reynolds_bound`` in
    ``ITTC_REYNOLDS_BOUNDS`` (the ITTC forms; module docstring). An unknown
    value raises ``ValueError``. Inputs: ``nu_r`` (6x1, relative to the
    water, SI), the couplings of ``surge_damping_couplings(surge_form)`` and
    the parameters of ``surge_damping_parameters(surge_form,
    ittc_reynolds_bound=...)`` (``ittc_reynolds_floor`` with ``"floor"`` on
    the ITTC forms), each by name. Outputs: the names
    of ``SURGE_DAMPING_OUTPUTS[surge_form]``: ``tau`` (the force on the
    vehicle), ``added_mass`` (A11 > 0), ``linear_coefficient`` (Xu < 0),
    ``quadratic_coefficient`` (Xuu < 0, a function of ``u_r`` on the ITTC
    forms), ``surge_damping_coefficient`` (``D_nl,11``, ``"exp_ittc"``). No
    value is checked inside the graph (``check_surge_damping_values``).
    """
    _check_selectors(surge_form, surge_blend, ittc_reynolds_bound)
    couplings = _COUPLINGS[surge_form]
    declared = surge_damping_parameters(surge_form, ittc_reynolds_bound=ittc_reynolds_bound)
    c, p = symbols(couplings), symbols(declared)
    nu_r = ca.SX.sym("nu_r", 6)
    if surge_form == "exp_ittc":
        outputs = _exponential(nu_r, c, p)
    else:
        outputs = _blended(nu_r, c, p, surge_form, surge_blend)
    names = [d.name for d in couplings] + [d.name for d in declared]
    symbols_in = {**c, **p}
    return ca.Function(
        "surge_damping",
        [nu_r, *[symbols_in[name] for name in names]],
        [outputs[name] for name in SURGE_DAMPING_OUTPUTS[surge_form]],
        ["nu_r", *names],
        list(SURGE_DAMPING_OUTPUTS[surge_form]),
    )


def check_surge_damping_values(values, *, surge_form, ittc_reynolds_bound="floor"):
    """Numbers for one form, couplings included, checked: ``{name: ca.DM}``
    (``check_values``: names, shapes, finite, ranges — ``ittc_reynolds_floor
    > 100`` included); for ``"exp_ittc"`` ``M(1,1)`` of ``mass_matrix`` must
    be positive."""
    declared = surge_damping_parameters(surge_form, ittc_reynolds_bound=ittc_reynolds_bound)
    numbers = check_values(_COUPLINGS[surge_form] + declared, values)
    if surge_form == "exp_ittc" and not float(numbers["mass_matrix"][0, 0]) > 0.0:
        raise ValueError("M(1,1) of mass_matrix must be positive")
    return numbers
