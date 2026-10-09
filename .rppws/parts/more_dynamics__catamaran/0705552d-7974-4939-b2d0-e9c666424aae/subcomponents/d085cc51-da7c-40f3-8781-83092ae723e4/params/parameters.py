from __future__ import annotations

# Otter restoring: two pontoons at the equilibrium draft (otter.m 121-122, 172-193).
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    hull_count = 2  # otter.m 172-193: two pontoons
    longitudinal_inertia_factor = 0.8  # otter.m 179 (I_L = 0.8 ...)
    longitudinal_center_of_flotation = -0.2  # otter.m 177, 192
    reference_point = [0.0, 0.0, 0.0]  # the body-frame origin (CO)
