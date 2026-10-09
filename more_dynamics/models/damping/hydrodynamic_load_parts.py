"""Members of the ``hydrodynamic_loads`` list slot (one list, not four
single sockets, so a class's own set of loads needs no ``none`` parts): the
time-constant damping forms, slender-body lift, and the ITTC surge
resistance line.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Through ``more_dynamics.models.damping``, ``.lift_drag`` and ``.surge_resistance``,
    which cites it line by line (sec. 6.4, eqs. 6.82-6.85, p. 125).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/Dmtrx.m;
    CRAFT/AUV/models/remus100.m 218 (submerged time constants);
    CRAFT/USV/models/otter.m 195-240 (surface calibration);
    LIBRARY/modeling/forceLiftDrag.m.

Author:    Enio Krizman
Date:      2026-10-08
"""

import math

import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.shared.wiring import function_from
from more_dynamics.models.damping.linear_damping import linear_damping_casadi, linear_damping_parameters
from more_dynamics.models.lift_drag.lift_drag import lift_drag_casadi, lift_drag_parameters
from more_dynamics.models.surge_resistance.surge_damping import surge_damping_casadi, surge_damping_parameters

# Declared (own) names, by hand against each block's actual remaining
# inputs: a name that also names a vehicle quantity (``weight``,
# ``center_of_gravity``, ``center_of_buoyancy``, ``mass_matrix``,
# ``restoring_matrix``, ``span``, ``water_density``, ``length``,
# ``kinematic_viscosity``, ``rigid_body_mass_matrix``, ``added_mass_matrix``)
# is never frozen here: the vehicle connects it by name.
_SUBMERGED_OWN = ("time_constants", "damping_ratios")
_SURFACE_OWN = ("max_forward_thrust", "max_speed", "time_constants", "damping_ratios", "yaw_damping_nonlinearity")
_LIFT_DRAG_OWN = ("planform_area", "parasitic_drag_coefficient", "oswald_efficiency")
_ITTC_OWN = ("time_constant", "form_factor", "crossover_speed", "surge_added_mass_factor")
_SMOOTH_SPEED_EPSILON = "smooth_speed_epsilon"  # numerical regularisation of the two guards (one value for both)


def _subset(declared, names):
    by_name = {d.name: d for d in declared}
    return tuple(by_name[n] for n in names)


def time_constant_damping_submerged(*, smooth_speed=False):
    """Time-constant linear damping, submerged form (Dmtrx.m; remus100.m
    218); couplings ``rigid_body_mass_matrix``, ``added_mass_matrix``.
    ``smooth_speed=True`` is the numerical guard of ``models/damping``: the
    speed of the surge fade has a finite derivative at rest."""
    return linear_damping_casadi(form="submerged", smooth_speed=smooth_speed)


def time_constant_damping_submerged_parameters(*, smooth_speed=False):
    own = _SUBMERGED_OWN + (_SMOOTH_SPEED_EPSILON,) if smooth_speed else _SUBMERGED_OWN
    return _subset(linear_damping_parameters("submerged", smooth_speed=smooth_speed), own)


def time_constant_damping_surface():
    """Time-constant linear damping, surface form — an MSS calibration from
    top speed and thrust (otter.m 195-240); couplings ``mass_matrix``,
    ``restoring_matrix``."""
    return linear_damping_casadi(form="surface")


def time_constant_damping_surface_parameters():
    return _subset(linear_damping_parameters("surface"), _SURFACE_OWN)


def hull_lift_drag(*, smooth_speed=False):
    """Slender-body lift and induced drag (forceLiftDrag.m).
    ``smooth_speed=True`` is the numerical guard of ``models/lift_drag``: the
    angle of attack and the flow direction have finite derivatives at rest."""
    return lift_drag_casadi(smooth_speed=smooth_speed)


def hull_lift_drag_parameters(*, smooth_speed=False):
    own = _LIFT_DRAG_OWN + (_SMOOTH_SPEED_EPSILON,) if smooth_speed else _LIFT_DRAG_OWN
    return _subset(lift_drag_parameters(smooth_speed=smooth_speed), own)


def surge_resistance_ittc():
    """The ITTC 1957 friction line plus a form factor (Fossen 2011, eqs.
    6.82-6.85, p. 125, through the surge-damping block)."""
    return surge_damping_casadi(surge_form="ittc")


def surge_resistance_ittc_parameters():
    """``ittc_reynolds_floor`` is declared conditionally (default bound
    ``"floor"``); included so a composition can set the floor."""
    declared = surge_damping_parameters("ittc", ittc_reynolds_bound="floor")
    return _subset(declared, _ITTC_OWN + ("ittc_reynolds_floor",))


def surge_resistance_ittc_residual():
    """The ITTC 1957 friction line plus the RA14 residual-resistance prior
    (Fossen 2011, eq. 6.82's own decomposition, p. 125; Radojcic et al. 2014
    Appendix 1, p. 24, through the surge-damping block)."""
    return surge_damping_casadi(surge_form="ittc_residual")


def surge_resistance_ittc_residual_parameters():
    """Own declared parameters -- the same set as ``"ittc"``: the RA14 table
    is a module constant of the published regression
    (``models/surge_resistance/ra14_residual_resistance.py``), not a composition value."""
    declared = surge_damping_parameters("ittc_residual", ittc_reynolds_bound="floor")
    return _subset(declared, _ITTC_OWN + ("ittc_reynolds_floor",))


MANOEUVRING_DAMPING_FORMS = ("time_constants", "linear_coupled", "linear_coupled_modulus")

SWAY_YAW_DERIVATIVES_PARAMETER = Parameter(
    "sway_yaw_derivatives", (4, 1), "N*s/m, N*s, N*m*s/m, N*m*s",
    "[Y_v, Y_r, N_v, N_r] of the coupled linear sway-yaw law tau_Y = Y_v v_r + Y_r r_r, "
    "tau_N = N_v v_r + N_r r_r (manoeuvring_damping='linear_coupled'); given, identified, or from "
    "the time-constant formula Y_v = -M22/T_sway, N_r = -M66/T_yaw with Y_r = N_v = 0 (unidentified, "
    "E-101 Q2 a)")
MODULUS_COEFFICIENTS_PARAMETER = Parameter(
    "modulus_coefficients", (2, 1), "N*s^2/m^2, N*m*s^2/m^2",
    "[Y_vv, N_rr] of a second-order sway/yaw term Y_vv |v_r| v_r, N_rr |r_r| r_r "
    "(manoeuvring_damping='linear_coupled_modulus'); zero (the linear_coupled law alone) until "
    "identified -- the exact modulus form (A-64 Section 3.2 row R7) is not read in this job and is "
    "not used by any default (rule 21)")


def coupled_manoeuvring_damping_surface(*, modulus=False):
    """Surge/heave/roll/pitch as the surface calibration (otter.m 195-207:
    ``Xu`` from the top speed and thrust, ``Zw``/``Kp``/``Mq`` from damping
    ratios and natural frequencies, unchanged from ``time_constant_damping_surface``);
    sway and yaw as a coupled linear law from ``sway_yaw_derivatives`` (E-101
    Q2 a: ``Y_r = N_v = 0`` until identified), optionally plus a second-order
    ``modulus_coefficients`` term (zero, unsourced placeholder, module
    docstring, ``modulus=True``)."""
    mass = ca.SX.sym("mass_matrix", 6, 6)
    restoring = ca.SX.sym("restoring_matrix", 6, 6)
    nu_r = ca.SX.sym("nu_r", 6)
    max_forward_thrust = ca.SX.sym("max_forward_thrust")
    max_speed = ca.SX.sym("max_speed")
    damping_ratios = ca.SX.sym("damping_ratios", 3)
    sway_yaw = ca.SX.sym("sway_yaw_derivatives", 4)
    y_v, y_r, n_v, n_r = sway_yaw[0], sway_yaw[1], sway_yaw[2], sway_yaw[3]
    zeta3, zeta4, zeta5 = damping_ratios[0], damping_ratios[1], damping_ratios[2]
    w3 = ca.sqrt(restoring[2, 2] / mass[2, 2])  # (Fossen 2011, eqs. 4.51-4.53, p. 68; otter.m 197)
    w4 = ca.sqrt(restoring[3, 3] / mass[3, 3])  # (otter.m 198)
    w5 = ca.sqrt(restoring[4, 4] / mass[4, 4])  # (otter.m 199)
    x_u = -max_forward_thrust / max_speed  # (otter.m 202)
    z_w = -2 * zeta3 * w3 * mass[2, 2]  # (Fossen 2011, eq. 6.78, p. 125; otter.m 204)
    k_p = -2 * zeta4 * w4 * mass[3, 3]  # (eq. 6.79; otter.m 205)
    m_q = -2 * zeta5 * w5 * mass[4, 4]  # (eq. 6.80; otter.m 206)
    zero = ca.SX(0.0)
    d = ca.vertcat(
        ca.horzcat(-x_u, zero, zero, zero, zero, zero),
        ca.horzcat(zero, -y_v, zero, zero, zero, -y_r),
        ca.horzcat(zero, zero, -z_w, zero, zero, zero),
        ca.horzcat(zero, zero, zero, -k_p, zero, zero),
        ca.horzcat(zero, zero, zero, zero, -m_q, zero),
        ca.horzcat(zero, -n_v, zero, zero, zero, -n_r),
    )
    tau = -d @ nu_r  # (E-101 Q2 a: the coupled sway/yaw law, Xu/Zw/Kp/Mq as otter.m 202-206)
    damping_derivatives = ca.vertcat(x_u, y_v, z_w, k_p, m_q, n_r)  # the diagonal part only
    inputs = {"nu_r": nu_r, "mass_matrix": mass, "restoring_matrix": restoring,
             "max_forward_thrust": max_forward_thrust, "max_speed": max_speed,
             "damping_ratios": damping_ratios, "sway_yaw_derivatives": sway_yaw}
    if modulus:
        modulus_coefficients = ca.SX.sym("modulus_coefficients", 2)
        # placeholder, zero by default (module docstring: not sourced, not used by any default)
        tau = tau + ca.vertcat(0.0, modulus_coefficients[0] * ca.fabs(nu_r[1]) * nu_r[1], 0.0, 0.0, 0.0,
                              modulus_coefficients[1] * ca.fabs(nu_r[5]) * nu_r[5])
        inputs["modulus_coefficients"] = modulus_coefficients
    return function_from("coupled_manoeuvring_damping_surface", inputs,
                         {"D": d, "tau": tau, "damping_derivatives": damping_derivatives})


def coupled_manoeuvring_damping_surface_parameters(*, modulus=False):
    """``max_forward_thrust``, ``max_speed``, ``damping_ratios`` (the surface
    form's own, unchanged), ``sway_yaw_derivatives``, and
    ``modulus_coefficients`` when ``modulus=True``."""
    surface_own = _subset(linear_damping_parameters("surface"), ("max_forward_thrust", "max_speed", "damping_ratios"))
    declared = surface_own + (SWAY_YAW_DERIVATIVES_PARAMETER,)
    return declared + ((MODULUS_COEFFICIENTS_PARAMETER,) if modulus else ())


HULL_PLANFORM_PARAMETERS = (
    Parameter("planform_fraction", (1, 1), "1", "planform area as a fraction of the rectangle length x diameter",
              0.0, 1.0, minimum_exclusive=True),
)
HULL_PARASITIC_DRAG_PARAMETERS = (
    Parameter("hull_frontal_drag_coefficient", (1, 1), "1", "drag coefficient Cd of the hull on its frontal area",
              0.0, minimum_exclusive=True),
)


def hull_planform_area():
    """``(length, span, planform_fraction) -> planform_area``: ``S = f L D``, the planform of a slender body as a
    fraction ``f`` of the rectangle length x diameter (``remus100.m`` 133: ``S = 0.7 L D``; ``span`` is the
    body's diameter ``D``)."""
    length, span, fraction = ca.SX.sym("length"), ca.SX.sym("span"), ca.SX.sym("planform_fraction")
    return function_from("hull_planform_area", {"length": length, "span": span, "planform_fraction": fraction},
                         {"planform_area": fraction * length * span})  # S = f L D (remus100.m 133)


def hull_parasitic_drag():
    """``(span, planform_area, hull_frontal_drag_coefficient) -> parasitic_drag_coefficient``: ``C_D0 = Cd pi (D/2)^2
    / S``, the drag of the frontal circle ``pi (D/2)^2`` with the coefficient ``Cd`` re-expressed on the planform area
    (``remus100.m`` 142-144, ``F_drag = 1/2 rho Cd pi b^2 = 1/2 rho C_D0 S`` with ``b = D/2``)."""
    span, planform = ca.SX.sym("span"), ca.SX.sym("planform_area")
    drag = ca.SX.sym("hull_frontal_drag_coefficient")
    return function_from("hull_parasitic_drag", {"span": span, "planform_area": planform,
                                                 "hull_frontal_drag_coefficient": drag},
                         {"parasitic_drag_coefficient": drag * math.pi * (span / 2) ** 2 / planform})  # (remus100.m 144)
