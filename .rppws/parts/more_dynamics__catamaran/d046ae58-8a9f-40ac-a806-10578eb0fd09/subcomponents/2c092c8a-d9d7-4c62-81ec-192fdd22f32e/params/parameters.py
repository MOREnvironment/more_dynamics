from __future__ import annotations

# Otter restoring: two pontoons at the equilibrium draft (otter.m 128-129, 181-202).
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    hull_count = 2  # otter.m 181-202: two pontoons
    longitudinal_inertia_factor = 0.8  # otter.m 188 (I_L = 0.8 ...)
    longitudinal_center_of_flotation = -0.2  # otter.m 186, 201
    reference_point = [0.0, 0.0, 0.0]  # the body-frame origin (CO)
    payload_trim_form = "lever_arm_and_heave_shift"  # otter.m 259-269 as MSS has it: the payload weight in heave a second time and the kinematics at the shifted attitude (F-MSS-31, F-MSS-24); the library's other form is "lever_arm_moment"
