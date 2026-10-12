from __future__ import annotations

# REMUS 100 hull loads: damping time constants, lift and drag (remus100.m).
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    time_constants = [20, 20, 1]  # remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217)
    damping_ratios = [0.3, 0.8]  # remus100.m:194, 195 [zeta4 zeta5]
    planform_area = 0.2128  # remus100.m:133 S = 0.7 * L_auv * D_auv
    parasitic_drag_coefficient = 0.05595961914206819  # remus100.m:143-144 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry)
    oswald_efficiency = 0.3  # coeffLiftDrag.m:54 e = 0.3
    planform_fraction = 0.7  # document: remus100.m:133 (S = 70% of the rectangle L_auv * D_auv)
    hull_frontal_drag_coefficient = 0.42  # document: remus100.m:142 (Cd), Allen et al. (2000) as MSS cites it, not read
