from __future__ import annotations

# Marie stern plane fin: Marie as MGM sets her up (more_generic_models vehicles/LAUV_Marie @ e86985cd). Source of every value: the comment beside it, SOURCE.md of the vehicle part has the table.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    fin_position = [-0.83, 0.0, 0.0]  # MGM actuators.yaml:46 (fin_x), :45 (fin_radius), :43 (fin_roll_angle)
    chord_axis = [1.0, 0.0, 0.0]  # MGM lifting_fins.py: chord = BODY x
    lift_axis = [0.0, 0.0, -1.0]  # MGM actuators.yaml:43 (fin_roll_angle), :44 (deflection_sign)
    fin_area = 0.0065  # MGM actuators.yaml:47
    servo_dynamics = 'none'  # MGM header actuators.yaml:20: "fins follow their command (no servo state)"
    servo_rate_limit = False  # MGM header actuators.yaml:20: the rate limit sits on the command in control.yaml
    servo_angle_limit = 'on_command'  # MGM range delta_max
    max_deflection = 0.4363323129985824  # MGM actuators.yaml:58
    inflow = 'rigid_point'  # MGM actuators.yaml:55
    flow_angle = 'small_angle'  # MGM actuators.yaml:56
    interference = 'none'  # MGM actuators.yaml:57
    section = 'linear_section'  # MGM drag_law none, stall_law none: lift only
    lift_law = 'finite_wing'  # MGM actuators.yaml:51
    lift_slope_form = 'low_aspect_ratio'  # MGM actuators.yaml:52
    fin_span = 0.0875  # MGM actuators.yaml:48
    aspect_ratio_factor = 2.0  # MGM actuators.yaml:49
    section_efficiency = 0.9  # MGM actuators.yaml:50
    zero_lift_drag = 0.0  # MGM drag_law none
