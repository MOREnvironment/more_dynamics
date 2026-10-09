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


def _subset(declared, names):
    by_name = {d.name: d for d in declared}
    return tuple(by_name[n] for n in names)


def time_constant_damping_submerged():
    """Time-constant linear damping, submerged form (Dmtrx.m; remus100.m
    218); couplings ``rigid_body_mass_matrix``, ``added_mass_matrix``."""
    return linear_damping_casadi(form="submerged")


def time_constant_damping_submerged_parameters():
    return _subset(linear_damping_parameters("submerged"), _SUBMERGED_OWN)


def time_constant_damping_surface():
    """Time-constant linear damping, surface form — an MSS calibration from
    top speed and thrust (otter.m 195-240); couplings ``mass_matrix``,
    ``restoring_matrix``."""
    return linear_damping_casadi(form="surface")


def time_constant_damping_surface_parameters():
    return _subset(linear_damping_parameters("surface"), _SURFACE_OWN)


def hull_lift_drag():
    """Slender-body lift and induced drag (forceLiftDrag.m)."""
    return lift_drag_casadi()


def hull_lift_drag_parameters():
    return _subset(lift_drag_parameters(), _LIFT_DRAG_OWN)


def surge_resistance_ittc():
    """The ITTC 1957 friction line plus a form factor (Fossen 2011, eqs.
    6.82-6.85, p. 125, through the surge-damping block)."""
    return surge_damping_casadi(surge_form="ittc")


def surge_resistance_ittc_parameters():
    """``ittc_reynolds_floor`` is declared conditionally (default bound
    ``"floor"``); included so a composition can set the floor."""
    declared = surge_damping_parameters("ittc", ittc_reynolds_bound="floor")
    return _subset(declared, _ITTC_OWN + ("ittc_reynolds_floor",))
