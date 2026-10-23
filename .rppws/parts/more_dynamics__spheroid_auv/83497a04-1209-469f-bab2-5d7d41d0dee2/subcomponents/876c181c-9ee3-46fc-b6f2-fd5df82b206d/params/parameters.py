from __future__ import annotations

# REMUS 100 rudder and stern plane, deflection-only form (remus100.m 229-246).
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    rudder_area = 0.0133  # remus100.m:184 A_r = 2 * S_fin (S_fin :180)
    stern_plane_area = 0.0133  # remus100.m:189 A_s = 2 * S_fin (S_fin :180)
    rudder_lift_coefficient = 0.5  # remus100.m:183 CL_delta_r
    stern_plane_lift_coefficient = 0.7  # remus100.m:188 CL_delta_s
    rudder_position = -0.8  # remus100.m:185 x_r = -a, a = L_auv / 2 (one geometry)
    stern_plane_position = -0.8  # remus100.m:190 x_s = -a, a = L_auv / 2 (one geometry)
    max_deflection = 0.3490658503988659  # remus100.m:110 delta_max = deg2rad(20)
