from __future__ import annotations

# REMUS 100 (estimated) hull loads: damping time constants, the body lift as one point force, the wing drag.
# Source of every value: the comment beside it, SOURCE.md of the vehicle gives its kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    time_constants = [20, 20, 1]  # s, [T1 T2 T6]: surge, sway and heave, yaw (design values of the published reference part)
    damping_ratios = [0.3, 0.8]  # [zeta4 zeta5]: roll and pitch (design values of the published reference part)
    planform_area = 0.2128  # m^2, 0.7 L D
    parasitic_drag_coefficient = 0.05595961914206819  # Cd pi (D/2)^2 / S with Cd 0.42
    oswald_efficiency = 0.3  # Oswald factor of the induced drag
    planform_fraction = 0.7  # S = 70 % of the rectangle L x D
    hull_frontal_drag_coefficient = 0.42  # Allen et al. (2000), cited second-hand, not read
    body_lift_form = "point_force_both_planes"  # Prestero 2001, eqs. 4.27-4.36, pp. 30-31
    body_lift_slope_method = "slender_body_hoerner"  # Prestero 2001, eqs. 4.30-4.32, p. 30
    body_lift_slope_per_degree = 0.003  # 1/deg, Prestero 2001, eq. 4.31, p. 30 (Hoerner, for 6.7 <= l/d <= 10)
    pressure_centre_method = "fraction_of_length_from_nose"  # Prestero 2001, eq. 4.35, p. 31
    pressure_centre_fraction = 0.65  # Prestero 2001, eq. 4.35, p. 31
    hull_drag_form = "wing_drag"  # the drag half of the wing block
    roll_damping_form = "linear"  # time-constant damping alone
