from __future__ import annotations

# REMUS 100 rudder and stern plane, deflection-only form (remus100.m 228-245).
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    rudder_area = 0.0133  # remus100.m:183 A_r = 2 * S_fin (S_fin :179)
    stern_plane_area = 0.0133  # remus100.m:188 A_s = 2 * S_fin (S_fin :179)
    rudder_lift_coefficient = 0.5  # remus100.m:182 CL_delta_r
    stern_plane_lift_coefficient = 0.7  # remus100.m:187 CL_delta_s
    rudder_position = -0.8  # remus100.m:184 x_r = -a, a = L_auv / 2 (one geometry)
    stern_plane_position = -0.8  # remus100.m:189 x_s = -a, a = L_auv / 2 (one geometry)
    max_deflection = 0.3490658503988659  # remus100.m:109 delta_max = deg2rad(20)
