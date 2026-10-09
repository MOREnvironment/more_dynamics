from __future__ import annotations

# REMUS 100 (Prestero 2001; MSS remus100.m @ cc07579): hull, mass, forms and the site.
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    integration_max_step = 0.05  # SIMremus100.m 60 (h = 0.05 s), the step of the MSS vehicle simulations
    latitude = 1.1073560310932362  # remus100.m 96-97, mu = deg2rad(63.446827)
    water_density = 1026  # remus100.m:98 rho = 1026, the one value (MSS: 1026 at imlay61.m:31 and forceLiftDrag.m:26, 1025 at crossFlowDrag.m:36)
    kinematic_viscosity = 1e-06  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
    length = 1.6  # remus100.m:131 L_auv, passed as L at :221
    diameter = 0.19  # remus100.m:132 D_auv, passed as B at :221
    body_density = 1054.7872613500267  # mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: 1025 at spheroid.m:35
    body_center_of_gravity = [0, 0, 0.02]  # remus100.m:137 r_bG
    roll_added_inertia_ratio = 0.3  # remus100.m:136 r44
    added_mass_form = "lamb_spheroid"  # imlay61.m 31-59
    coriolis_form = "kirchhoff_full"  # m2c.m 33-48 (every term kept)
    current_form = "none"  # still water
