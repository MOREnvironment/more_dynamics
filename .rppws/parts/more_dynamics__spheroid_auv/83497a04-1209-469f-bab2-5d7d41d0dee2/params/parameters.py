from __future__ import annotations

# REMUS 100 (Prestero 2001; MSS remus100.m @ 99fdf4d): hull, mass, forms and the site.
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    integration_max_step = 0.05  # SIMremus100.m 60 (h = 0.05 s), the step of the MSS vehicle simulations
    length = 1.6  # remus100.m:132 L_auv, passed as L at :222
    diameter = 0.19  # remus100.m:133 D_auv, passed as B at :222
    body_density = 1054.7872613500267  # mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: 1025 at spheroid.m:35-36
    body_center_of_gravity = [0, 0, 0.02]  # remus100.m:138 r_bG
    roll_added_inertia_ratio = 0.3  # remus100.m:137 r44
    added_mass_form = "lamb_spheroid"  # imlay61.m 31-60
    coriolis_form = "munk_couplings_removed"  # remus100.m 208-211 (MSS omits the Munk-moment couplings of the added-mass Coriolis: its reduced model keeps no measured manoeuvring derivatives to balance them; with every term kept this part spins open-loop, SOURCE.md)
    current_form = "none"  # still water
    body_mass = 31.9  # document: remus100.m:3
