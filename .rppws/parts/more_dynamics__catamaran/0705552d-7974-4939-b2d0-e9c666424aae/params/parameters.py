from __future__ import annotations

# Otter USV, catamaran (MSS otter_legacy.m @ 6e78e98): the two pontoons, mass, added-mass derivatives and the site.
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    integration_max_step = 0.05  # SIMremus100.m 60 (h = 0.05 s), the step of the MSS vehicle simulations
    length = 2.0  # otter_legacy.m 104-105, 117-120
    beam = 1.08  # otter_legacy.m 104-105, 117-120
    pontoon_beam = 0.25  # otter_legacy.m 117-120
    pontoon_lateral_offset = 0.395  # otter_legacy.m 117-120
    pontoon_block_coefficient = 0.4  # otter_legacy.m 117-120
    pontoon_waterplane_coefficient = 0.75  # otter_legacy.m 117-120
    hull_mass = 55.0  # otter_legacy.m 106
    payload_mass = 25.0  # generate_catamaran_mss.m 61-62 (payload mass)
    hull_center_of_gravity = [0.2, 0.0, -0.2]  # otter_legacy.m 107
    payload_position = [0.05, 0.0, -0.35]  # generate_catamaran_mss.m 61-62 (payload position)
    radii_of_gyration = [0.4, 0.25, 0.25]  # otter_legacy.m 108-110
    added_mass_coefficients = [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]  # otter_legacy.m 167-172
    added_mass_form = "scaled_derivatives_loaded_surge"  # otter_legacy.m 166-174 (surge from m + mP, line 162)
    coriolis_form = "kirchhoff_full"  # m2c.m 33-48 (every term kept)
    current_form = "none"  # still water
