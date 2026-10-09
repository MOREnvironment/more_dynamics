"""Single-screw propeller with Wageningen B-series coefficients (CasADi).

``propeller_parameters(...)`` declares the numbers (name, shape, SI unit,
meaning, admissible range); ``propeller_casadi(open_water_model=...)``
builds ``(n, nu_r, <parameters by name>) -> (tau, shaft_axis,
coefficients_in_use)``.

Equations (keys in References), ``n`` in rpm saturated at ``n_max``,
``n_p = n / 60`` (MSS ``remus100.m`` 115):

* ``Va = (1 - w) U_r`` (``remus100.m`` 150, there ``Va = 0.944 U_r``);
* ``"bollard"``: ``X = rho D^4 K_T0 n_p |n_p|``, ``K = rho D^5 K_Q0 n_p |n_p|``
  (``remus100.m`` 175-176; Fossen 2011, eq. 12.265, p. 411);
* ``"linearized"``: for ``n_p > 0``, ``X = rho D^4 (K_T0 n_p |n_p| +
  (K_Tmax - K_T0) / Ja_max (Va / D) |n_p|)`` and the same for ``K`` with
  ``D^5`` and ``K_Q`` (``remus100.m`` 167-172); for ``n_p <= 0`` the bollard
  form (``remus100.m`` 173-176);
* ``"full"``: ``J = Va / (|n_p| D)``, ``X = rho D^4 K_T(J) n_p |n_p|``,
  ``K = rho D^5 K_Q(J) n_p |n_p|`` with the Wageningen polynomial
  (``wageningen_kt_kq.py``); ``J`` clipped to ``[0, 2 Ja_max]`` and the
  denominator bounded below at 1e-6 (guards of this module); ``K_T0``,
  ``K_Q0`` for ``n_p <= 0``;
* body wrench: thrust ``(1 - t) X`` along the shaft axis at ``position``
  (``remus100.m`` 249; moment ``r x f``, Fossen 2011, eq. 12.226, p. 400),
  plus a torque ``scale K`` about the shaft axis (``remus100.m`` 252, where
  the scale is ``1/10``; a pure couple, no ``r x f`` part).

The ``"linearized"`` model equals MSS ``remus100.m`` 148-176, 249 and 252.
``open_water_coefficients`` (a selector): ``"polynomial"`` (default) takes
``K_T0, K_Q0, K_Tmax, K_Qmax`` from the Wageningen polynomial at ``J = 0`` and
``J = Ja_max``; ``"given"`` declares them as the parameter
``thrust_torque_coefficients`` (``remus100.m`` 157-161 prints them, rounded
to four digits). Deviations from MSS: the shaft axis and position are
parameters (``remus100.m`` puts the thrust on the x axis at the CO).

Domain of the advance number. ``J = Va / (n_p D)`` is ``>= 0`` on the
forward branch, where both models use it: ``Va = (1 - w) U_r >= 0``
(``remus100.m`` 150) and ``n_p > 0`` (``remus100.m`` 167). ``Ja_max`` is the
``J`` at the top speed and the top shaft speed, ``Ja_max = (1 - w) U_max /
(D n_max / 60)`` (``remus100.m`` 152-153); it is declared ``> 0`` (the
linearised slope divides by it). The upper end of the range on which the
Wageningen regression was fitted is not stated in MSS (``wageningen.m``
1-31) and its primary source (Bernitsas, Ray and Kinley 1981) is not read
here, so no upper bound is declared. ``"full"``'s clip of ``J`` to ``[0, 2
Ja_max]`` and the ``"linearized"`` form are comparison methods of this
module and of ``remus100.m``, not a statement of that range.

``thrust_deduction`` and ``wake_fraction`` are fractions of the thrust and
of the inflow speed, declared in ``[0, 1]`` (the same range as the outboard
motor's ``thrust_deduction`` and ``propwash_factor``; the equations need no
exclusive endpoint).

``roll_moment_scale`` is an identified coefficient or a comparison form, not a
physics value: ``remus100.m`` 252 divides the shaft torque by 10 "to match
exp. results"; the physical shaft reaction torque is the scale 1.

Frames: the shaft axis is ``R_zyx(orientation)[:, 0]`` in BODY
(``more_transformations.more_casadi_transformations``
``MatrixTransforms.Rzyx_explicit``; Fossen 2011, eq. 2.18, p. 22).

Motor lag and the thrust-command input type are states / interfaces of the
assembly, not part of this map.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 2,
    eq. 2.18, p. 22; Ch. 12, eq. 12.226, p. 400; eq. 12.265, p. 411.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 110, 115,
    127, 148-176, 249, 252.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter, symbols

from more_dynamics.models.shared.force_producer_common import safe_norm, saturate, wrench
from .wageningen_kt_kq import wageningen_expression, wageningen_parameters

OPEN_WATER_MODELS = ("linearized", "full", "bollard")
OPEN_WATER_COEFFICIENTS = ("polynomial", "given")
MIN_ADVANCE_DENOMINATOR = 1e-6  # m/s: lower bound of |n_p| D in J (guard of this module)

PROPELLER_OUTPUTS = ("tau", "shaft_axis", "coefficients_in_use")


def _positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


_COMMON = (
    _positive("propeller_diameter", "m", "propeller diameter D"),
    _positive("max_shaft_speed", "rpm", "shaft speed limit n_max"),
    Parameter("thrust_deduction", (1, 1), "1", "thrust deduction t: thrust (1 - t) X", 0.0, 1.0),
    Parameter("wake_fraction", (1, 1), "1", "wake fraction w: Va = (1 - w) U_r", 0.0, 1.0),
    *wageningen_parameters(),
    _positive("max_advance_number", "1",
              "Ja_max = (1 - w) U_max / (D n_max / 60), J >= 0 on the forward branch; "
              "linearised slope and J clip (module docstring: domain)"),
    Parameter("roll_moment_scale", (1, 1), "1",
              "scale of the shaft torque K applied about the shaft axis; identified coefficient "
              "or comparison form, not a physics value (physical: 1; module docstring)"),
    _positive("water_density", "kg/m^3", "water density rho"),
    Parameter("position", (3, 1), "m", "CO -> propeller, body axes (FRD)"),
    Parameter("orientation", (3, 1), "rad", "shaft orientation (roll, pitch, yaw), zyx Euler angles"),
)
_GIVEN = Parameter("thrust_torque_coefficients", (4, 1), "1", "[K_T0, K_Q0, K_Tmax, K_Qmax]")


def _check_selectors(open_water_model, open_water_coefficients, clip_advance_ratio):
    if open_water_model not in OPEN_WATER_MODELS:
        raise ValueError(f"open_water_model must be one of {OPEN_WATER_MODELS}, got {open_water_model!r}")
    if open_water_coefficients not in OPEN_WATER_COEFFICIENTS:
        raise ValueError(f"open_water_coefficients must be one of {OPEN_WATER_COEFFICIENTS}, "
                         f"got {open_water_coefficients!r}")
    if not isinstance(clip_advance_ratio, bool):
        raise ValueError("clip_advance_ratio must be True or False")


def propeller_parameters(*, open_water_coefficients="polynomial"):
    """The declared parameter set: a tuple of ``Parameter`` in the order of
    the block's inputs; ``thrust_torque_coefficients`` (4x1) only with
    ``open_water_coefficients="given"``."""
    _check_selectors("linearized", open_water_coefficients, False)
    return _COMMON + ((_GIVEN,) if open_water_coefficients == "given" else ())


def _open_water(p, coefficients, model, clip_advance_ratio, n_rps, advance_speed):
    """Shaft thrust ``X_prop`` and torque ``K_prop`` of the chosen model."""
    kt0, kq0, kt_max, kq_max = coefficients[0], coefficients[1], coefficients[2], coefficients[3]
    rho, d = p["water_density"], p["propeller_diameter"]
    quadratic = ca.fabs(n_rps) * n_rps  # n_p |n_p| (Fossen 2011, eq. 12.265, p. 411)

    if model == "bollard":  # (remus100.m 175-176)
        return rho * d**4 * kt0 * quadratic, rho * d**5 * kq0 * quadratic

    if model == "linearized":  # (remus100.m 167-176)
        kt_slope = (kt_max - kt0) / p["max_advance_number"]  # (K_Tmax - K_T0) / Ja_max (remus100.m 170)
        kq_slope = (kq_max - kq0) / p["max_advance_number"]  # (remus100.m 172)
        inflow = (advance_speed / d) * ca.fabs(n_rps)  # (Va / D) |n_p| (remus100.m 170, 172)
        thrust = ca.if_else(  # (remus100.m 169-170, 175)
            n_rps > 0.0,
            rho * d**4 * (kt0 * quadratic + kt_slope * inflow),
            rho * d**4 * kt0 * quadratic,
        )
        torque = ca.if_else(  # (remus100.m 171-172, 176)
            n_rps > 0.0,
            rho * d**5 * (kq0 * quadratic + kq_slope * inflow),
            rho * d**5 * kq0 * quadratic,
        )
        return thrust, torque

    # "full": J = Va / (|n_p| D) clipped to [0, 2 Ja_max] (guards of this module),
    # then the polynomial's own clip only with ``clip_advance_ratio``;
    # K_T(0), K_Q(0) when n <= 0; zero force at n = 0 through |n| n.
    advance_number = advance_speed / ca.fmax(ca.fabs(n_rps) * d, MIN_ADVANCE_DENOMINATOR)  # J = Va / (n D) (wageningen.m 6; remus100.m 147); denominator bounded below (guard)
    advance_number = saturate(advance_number, 0.0, 2.0 * p["max_advance_number"])  # guard of this module
    kt_j, kq_j = wageningen_expression(advance_number, p["pitch_diameter_ratio"], p["blade_area_ratio"],
                                       p["blade_count"], clip_advance_ratio)
    kt = ca.if_else(n_rps <= 0.0, kt0, kt_j)  # K_T0 astern (remus100.m 175, bollard form)
    kq = ca.if_else(n_rps <= 0.0, kq0, kq_j)  # K_Q0 astern (remus100.m 176)
    return rho * d**4 * kt * quadratic, rho * d**5 * kq * quadratic  # (Fossen 2011, eq. 12.265, p. 411)


def propeller_casadi(*, open_water_model, open_water_coefficients="polynomial", clip_advance_ratio=False):
    """The block: ``(n, nu_r, <parameters by name>) -> (tau, shaft_axis,
    coefficients_in_use)``.

    Contract
    --------
    Keywords choose the graph: ``open_water_model`` in ``OPEN_WATER_MODELS``
    (required), ``open_water_coefficients`` in ``OPEN_WATER_COEFFICIENTS``,
    ``clip_advance_ratio`` a bool (passed to the Wageningen polynomial). An
    unknown value raises ``ValueError``. Inputs: ``n`` (1x1, rpm, saturated
    at ``max_shaft_speed``), ``nu_r`` (6x1; advance speed ``(1 - w)
    |nu_r[0:3]|``) and the parameters of ``propeller_parameters(...)`` by
    name. Outputs: ``tau`` (6x1, N and N m, BODY, about the CO),
    ``shaft_axis`` (3x1), ``coefficients_in_use`` (4x1, ``[K_T0, K_Q0,
    K_Tmax, K_Qmax]``).
    """
    _check_selectors(open_water_model, open_water_coefficients, clip_advance_ratio)
    declared = propeller_parameters(open_water_coefficients=open_water_coefficients)
    p = symbols(declared)
    n = ca.SX.sym("n", 1)
    nu_r = ca.SX.sym("nu_r", 6)
    geometry = (p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"])
    if open_water_coefficients == "given":
        coefficients = p["thrust_torque_coefficients"]
    else:
        kt0, kq0 = wageningen_expression(0.0, *geometry, clip_advance_ratio)  # K_T, K_Q at J = 0 (remus100.m 156)
        kt_max, kq_max = wageningen_expression(p["max_advance_number"], *geometry, clip_advance_ratio)  # at J = Ja_max (remus100.m 159)
        coefficients = ca.vertcat(kt0, kq0, kt_max, kq_max)

    n_rps = saturate(n, -p["max_shaft_speed"], p["max_shaft_speed"]) / 60.0  # n_p = sat(n, n_max) / 60 (remus100.m 115)
    advance_speed = (1.0 - p["wake_fraction"]) * safe_norm(nu_r[0:3])  # Va = (1 - w) U_r (remus100.m 127, 150)
    thrust, torque = _open_water(p, coefficients, open_water_model, clip_advance_ratio, n_rps, advance_speed)

    axis = MatrixTransforms.Rzyx_explicit(p["orientation"])[:, 0]  # shaft axis R_zyx[:, 0] (Fossen 2011, eq. 2.18, p. 22)
    force = axis * ((1.0 - p["thrust_deduction"]) * thrust)  # (1 - t) X along the shaft (remus100.m 249)
    shaft_torque = axis * (p["roll_moment_scale"] * torque)  # scale K about the shaft axis (remus100.m 252)
    # [f; r x f] of the thrust (Fossen 2011, eq. 12.226, p. 400) plus the shaft torque as a couple (remus100.m 252)
    tau = wrench(force, p["position"]) + ca.vertcat(ca.DM.zeros(3), shaft_torque)

    names = [d.name for d in declared]
    outputs = {"tau": tau, "shaft_axis": axis, "coefficients_in_use": coefficients}
    return ca.Function(
        "propeller",
        [n, nu_r, *[p[name] for name in names]],
        [outputs[name] for name in PROPELLER_OUTPUTS],
        ["n", "nu_r", *names],
        list(PROPELLER_OUTPUTS),
    )
