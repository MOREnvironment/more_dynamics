from __future__ import annotations

# Otter USV, catamaran (MSS otter.m @ cc07579): the two pontoons, mass, added-mass derivatives and the site.
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    integration_max_step = 0.05  # SIMremus100.m 60 (h = 0.05 s), the step of the MSS vehicle simulations
    gravity = 9.81  # otter.m 90-91
    water_density = 1025.0  # otter.m 90-91
    kinematic_viscosity = 1e-06  # cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s; Fossen 2011, p. 125, below eq. 6.85 (20 degC)
    length = 2.0  # otter.m 92-93, 104-107
    beam = 1.08  # otter.m 92-93, 104-107
    pontoon_beam = 0.25  # otter.m 104-107
    pontoon_lateral_offset = 0.395  # otter.m 104-107
    pontoon_block_coefficient = 0.4  # otter.m 104-107
    pontoon_waterplane_coefficient = 0.75  # otter.m 104-107
    hull_mass = 55.0  # otter.m 94-98
    payload_mass = 25.0  # generate_catamaran_mss.m 61-62 (payload mass)
    hull_center_of_gravity = [0.2, 0.0, -0.2]  # otter.m 94-98
    payload_position = [0.05, 0.0, -0.35]  # generate_catamaran_mss.m 61-62 (payload position)
    radii_of_gyration = [0.4, 0.25, 0.25]  # otter.m 94-98
    added_mass_coefficients = [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]  # otter.m 152-157
    added_mass_form = "scaled_derivatives"  # otter.m 152-159
    coriolis_form = "kirchhoff_full"  # m2c.m 33-48 (every term kept)
    current_form = "none"  # still water
