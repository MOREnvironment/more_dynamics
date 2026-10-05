from __future__ import annotations

# Grethe (Mariner 5) values. Source of every value, with page:
# scripts/vehicle_models/params/grethe_mariner5.yaml (key in each comment).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    length = 5.2  # hull.length
    beam = 2.15  # hull.beam
    draft = 0.3  # hull.draft
    block_coefficient = 0.233  # hydrostatics.block_coefficient
    waterplane_coefficient = 0.8  # hydrostatics.waterplane_coefficient
    water_density = 1025.0  # hydrostatics.water_density
    gravity = 9.81  # hydrostatics.gravity
    center_of_gravity = [0.0, 0.0, 0.025]  # mass.center_of_gravity
    longitudinal_center_of_flotation = -0.3  # hydrostatics.longitudinal_center_of_flotation
    coefficient_scales = [1.0, 0.7, 1.0]  # hydrostatics.restoring_scales.workspace_rppws_value (OPEN: library 1.0 in roll)
    reference_point = [0.0, 0.0, 0.0]  # hydrostatics.restoring_reference_point
