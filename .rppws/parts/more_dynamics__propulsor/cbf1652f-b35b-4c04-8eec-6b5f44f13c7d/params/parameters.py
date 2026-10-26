from __future__ import annotations

# Grethe outboard: the one steerable electric outboard of Grethe as a Propulsor, the values of the model file grethe.m
# (MSS fork branch grethe-vehicle @ 679a6e8; lines beside each value). The kind of each value is in SOURCE.md.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    drive = 'thrust_command_lag'  # the thrust command (N) through a first-order lag, the actual thrust the model file takes as ui(1) (grethe.m 38-48; SIMgrethe.m 257)
    steering = 'servo_lag_rate_angle'  # the steering angle command through a lag with rate and angle limits (grethe.m 265-267; SIMgrethe.m 258-259)
    mounting = 'azimuth'  # thrust force T [cos d, sin d, 0] at the motor position (grethe.m 395-396)
    inflow = 'horizontal_speed_floored'  # the power limit is set by the horizontal speed, floored (grethe.m 391)
    law = 'thrust_power_limited'  # T = sat(command, +-min(T_max, P / U)) (grethe.m 391-393)
    housing = 'leg_normal_force'  # the normal force of the outboard leg (grethe.m 400-408)
    steering_convention = 'port_positive'  # the same formula the model file names starboard_positive (a positive angle = thrust force to starboard, bow to port; grethe.m 38-48, 276)
    shaft_torque_scale = 0.0  # the model file carries no shaft torque on the body (zero with its reason)
    thrust_time_constant = 0.3  # s: T_thrust (grethe.m 268)
    max_deflection = 0.6108652381980153  # rad: delta_max = deg2rad(35) (grethe.m 265)
    max_rate = 1.0471975511965976  # rad/s: ddelta_max = deg2rad(60) (grethe.m 266)
    time_constant = 0.5  # s: T_delta (grethe.m 267)
    position = [-2.24, 0.0, 0.55]  # m: r_P, the motor position from the CO, body axes (grethe.m 264)
    speed_floor = 0.1  # m/s: lowest horizontal speed of the power limit (grethe.m 391)
    max_thrust = 511.5  # N: T_max (grethe.m 262)
    delivered_power = 1120.0  # W: P_prop, power delivered to the water (grethe.m 263)
    leg_area = 0.04  # m^2: A_leg (grethe.m 269)
    leg_aspect_ratio = 1.5  # 1: AR_leg (grethe.m 270)
