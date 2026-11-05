from __future__ import annotations

# REMUS 100 (estimated): hull, mass, forms and the site. Every added-mass Coriolis term is kept; the body lift acts in both
# planes as one point force; the fins take their lift slope from their geometry. Values of the one-value REMUS 100 set.
# Source of every value: the comment beside it, SOURCE.md gives its kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    integration_max_step = 0.05  # s, the step of the vehicle simulations
    length = 1.6  # m, the one-value REMUS 100 set (the same hull as the published reference part)
    diameter = 0.19  # m, the one-value REMUS 100 set (the same hull as the published reference part)
    body_density = 1054.7872613500267  # kg/m^3, the mean density of the 31.9 kg spheroid: mass / (4/3 pi a b^2) (derived)
    body_center_of_gravity = [0, 0, 0.02]  # m, CO -> CG
    roll_added_inertia_ratio = 0.3  # roll added inertia as a fraction of the roll inertia of the displaced fluid
    added_mass_form = "lamb_spheroid"  # ideal-fluid prolate spheroid
    coriolis_form = "kirchhoff_full"  # every term of the added-mass Coriolis matrix (Fossen 2011, eq. 6.43, p. 120); the Munk moment is kept and balanced by the body lift and the fins
    current_form = "none"  # still water
    body_mass = 31.9  # kg
