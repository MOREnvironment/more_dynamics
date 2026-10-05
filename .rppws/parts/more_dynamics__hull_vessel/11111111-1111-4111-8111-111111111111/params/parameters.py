from __future__ import annotations

# Grethe (Mariner 5) values. Source of every value, with page:
# scripts/vehicle_models/params/grethe_mariner5.yaml (key in each comment).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    length = 5.2  # hull.length
    beam = 2.15  # hull.beam
    draft = 0.3  # hull.draft
    block_coefficient = 0.233  # hydrostatics.block_coefficient
    radii_of_gyration = [0.35, 0.25, 0.25]  # mass.radii_of_gyration_scale
    added_mass_coefficients = [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]  # mass.added_mass_scales.workspace_rppws_value (OPEN: library -1.2 in yaw)
    center_of_gravity = [0.0, 0.0, 0.025]  # mass.center_of_gravity
    center_of_buoyancy = [0.0, 0.0, 0.0]  # mass.center_of_buoyancy
    water_density = 1025.0  # hydrostatics.water_density
