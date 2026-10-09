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
