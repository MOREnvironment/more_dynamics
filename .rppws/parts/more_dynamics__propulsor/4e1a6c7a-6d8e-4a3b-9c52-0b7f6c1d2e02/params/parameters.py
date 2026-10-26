from __future__ import annotations

# Otter port thruster: one of the two fixed thrusters of the Otter as a Propulsor (port).
# Source of every value: the comment beside it, SOURCE.md gives its kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    drive = 'speed_command_rad_s'  # the shaft-speed command in rad/s, saturated (otter.m:230)
    steering = 'fixed'  # a fixed thruster
    mounting = 'fixed_direction'  # thrust along x_b at a lever arm (otter.m:139-140, 222-224)
    inflow = 'none'  # otter.m takes nothing from the water for the thrust
    law = 'quadratic_asymmetric'  # T = k n |n| with k_pos ahead and k_neg astern (otter.m:232-238)
    housing = 'none'  # no housing load
    min_command = -101.73666550381692  # rad/s: n_min = -sqrt(0.5 * 13.6 * g / k_neg), g = 9.81 (otter.m:146, mssConstants.m:14)
    max_command = 103.93086427398005  # rad/s: n_max = sqrt(0.5 * 24.4 * g / k_pos), g = 9.81 (otter.m:145, mssConstants.m:14)
    position = [0.0, -0.395, 0.0]  # m: lever arm l = -y_pont, y_pont = 0.395 (otter.m:112, 139-140), no x or z offset in the model
    direction = [1.0, 0.0, 0.0]  # thrust along the body x axis (otter.m:242)
    positive_thrust_coefficient = 0.01108  # N s^2/rad^2: otter.m:143 k_pos = 0.02216 / 2 (one thruster)
    negative_thrust_coefficient = 0.006445  # N s^2/rad^2: otter.m:144 k_neg = 0.01289 / 2 (one thruster)
    shaft_torque_scale = 0.0  # the model carries no thruster torque on the body (zero with its reason)
