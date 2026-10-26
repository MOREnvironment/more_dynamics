from __future__ import annotations

# REMUS 100 (MSS) stern plane fin: one of the two fins of the stern plane (remus100.m 188-190: CL_delta_s 0.7, A_s = 2 S_fin: two fins of S_fin).
# Source of every value: the comment beside it, SOURCE.md gives its kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    fin_position = [-0.8, 0.0, 0.0]  # m, remus100.m:185, 190 x = -a, a = L_auv / 2 (one geometry: derived, see SOURCE.md)
    chord_axis = [1.0, 0.0, 0.0]  # unit chord direction, BODY FRD (convention; the drag acts along it, remus100.m:239-240)
    lift_axis = [0.0, 0.0, -1.0]  # unit positive-lift direction, BODY FRD (the force of a positive deflection, remus100.m:246: negative on its axis)
    fin_area = 0.00665  # m^2, remus100.m:180 S_fin (the pair has A = 2 * S_fin, remus100.m:184, 189)
    servo_dynamics = 'none'  # the deflection is the command (remus100.m:114-115)
    servo_rate_limit = False  # no rate limit in remus100.m
    servo_angle_limit = 'on_command'  # the command saturated at delta_max (remus100.m:114-115)
    max_deflection = 0.3490658503988659  # rad, remus100.m:110 delta_max = deg2rad(20)
    inflow = 'translational'  # the fin sees the vehicle's relative translation (remus100.m:235-236)
    flow_angle = 'none'  # no flow-angle correction (remus100.m:235-246)
    interference = 'none'  # no interference factors
    section = 'quadratic_drag'  # lift slope * delta, drag = lift slope * delta^2 (remus100.m:239-240, 243, 246)
    lift_slope = 0.7  # 1/rad, remus100.m:188 CL_delta_s
