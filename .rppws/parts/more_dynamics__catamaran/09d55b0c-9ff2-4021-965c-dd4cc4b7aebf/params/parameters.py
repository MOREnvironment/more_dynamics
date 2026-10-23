from __future__ import annotations


class ComponentParameters:
    integration_max_step = 0.05
    length = 2.0
    beam = 1.08
    pontoon_beam = 0.25
    pontoon_lateral_offset = 0.395
    pontoon_block_coefficient = 0.4
    pontoon_waterplane_coefficient = 0.75
    hull_mass = 55.0
    payload_mass = 25.0
    hull_center_of_gravity = [0.2, 0.0, -0.2]
    payload_position = [0.05, 0.0, -0.35]
    radii_of_gyration = [0.4, 0.25, 0.25]
    added_mass_coefficients = [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]
    added_mass_form = 'scaled_derivatives_loaded_surge'
    coriolis_form = 'kirchhoff_full'
    current_form = 'none'
