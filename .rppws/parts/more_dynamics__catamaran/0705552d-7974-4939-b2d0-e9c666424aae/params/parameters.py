from __future__ import annotations

# Otter USV, catamaran (MSS otter.m @ 99fdf4d): the two pontoons, mass, added-mass derivatives and the site.
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    integration_max_step = 0.05  # SIMremus100.m 60 (h = 0.05 s), the step of the MSS vehicle simulations
    length = 2.0  # otter.m 98-99, 111-114
    beam = 1.08  # otter.m 98-99, 111-114
    pontoon_beam = 0.25  # otter.m 111-114
    pontoon_lateral_offset = 0.395  # otter.m 111-114
    pontoon_block_coefficient = 0.4  # otter.m 111-114
    pontoon_waterplane_coefficient = 0.75  # otter.m 111-114
    hull_mass = 55.0  # otter.m 100
    payload_mass = 25.0  # generate_catamaran_mss.m 61-62 (payload mass)
    hull_center_of_gravity = [0.2, 0.0, -0.2]  # otter.m 101
    payload_position = [0.05, 0.0, -0.35]  # generate_catamaran_mss.m 61-62 (payload position)
    radii_of_gyration = [0.4, 0.25, 0.25]  # otter.m 102-104
    added_mass_coefficients = [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]  # otter.m 161-166
    added_mass_form = "scaled_derivatives_loaded_surge"  # otter.m 160-168 (surge from m + mP, line 162)
    coriolis_form = "kirchhoff_full"  # m2c.m 33-48 (every term kept)
    current_form = "none"  # still water
