from __future__ import annotations

# REMUS 100 (estimated) propulsor: the propeller of the REMUS 100 as one Propulsor (fixed shaft, speed command in rpm, linearised open-water law), with the estimated roll reaction.
# Source of every value: the comment beside it, SOURCE.md gives its kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    drive = 'speed_command_rpm'  # the shaft-speed command in rpm, saturated (remus100.m:116)
    steering = 'fixed'  # a fixed shaft
    mounting = 'fixed_euler'  # shaft along x_b through the CO (remus100.m:250-253)
    inflow = 'speed_magnitude'  # Va = (1 - w) U_r, thrust (1 - t) X (remus100.m:151, 250)
    law = 'open_water_linearized_given'  # linearised open-water law with the four printed coefficients (remus100.m:168-178)
    housing = 'none'  # no housing load
    min_command = -1525  # rpm: -n_max, remus100.m:111, 116 (sat symmetric)
    max_command = 1525  # rpm: n_max, remus100.m:111
    thrust_deduction = 0.1  # remus100.m:150 t_prop
    wake_fraction = 0.05600000000000005  # remus100.m:151 Va = 0.944 * U_r, w = 1 - 0.944
    position = [0, 0, 0]  # m, remus100.m:250-253 thrust on x_b through the CO (no moment arm)
    orientation = [0, 0, 0]  # rad, remus100.m:250, 253 shaft along x_b
    diameter = 0.14  # m, remus100.m:149 D_prop
    pitch_diameter_ratio = 1  # remus100.m:157 wageningen(0,1,0.718,3), second argument
    blade_area_ratio = 0.718  # remus100.m:156-157 blade-area ratio 0.718
    blade_count = 3  # remus100.m:156-157 3 blades
    max_advance_number = 0.6632  # remus100.m:154 Ja_max
    thrust_torque_coefficients = [0.4566, 0.07, 0.1798, 0.0312]  # remus100.m:158, 159, 161, 162 [KT_0 KQ_0 KT_max KQ_max]
    shaft_torque_scale = 0.35  # scale of the shaft-torque reaction on the hull, Prestero 2001, Sec. 4.6.2, p. 33; calibrated to a steady roll of -5.3 deg in the generic-models REMUS actuator set, carried over here, not re-calibrated for this plant
