from __future__ import annotations

# Marie propulsor: Marie as MGM sets her up (more_generic_models vehicles/LAUV_Marie @ e86985cd). Source of every value: the comment beside it, SOURCE.md of the vehicle part has the table.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    drive = 'speed_command_rpm'  # selector: MGM PowerLawPropeller is a fixed shaft, rpm command, X = k|n|^p, no inflow term
    steering = 'fixed'  # selector: MGM PowerLawPropeller is a fixed shaft, rpm command, X = k|n|^p, no inflow term
    mounting = 'fixed_direction'  # selector: MGM PowerLawPropeller is a fixed shaft, rpm command, X = k|n|^p, no inflow term
    inflow = 'none'  # selector: MGM PowerLawPropeller is a fixed shaft, rpm command, X = k|n|^p, no inflow term
    law = 'power_law'  # selector: MGM PowerLawPropeller is a fixed shaft, rpm command, X = k|n|^p, no inflow term
    housing = 'none'  # selector: MGM PowerLawPropeller is a fixed shaft, rpm command, X = k|n|^p, no inflow term
    steering_convention = 'port_positive'  # plugin default, no MGM counterpart
    shaft_torque_scale = 1.0  # MGM power_law.py wrench: reaction = -shaft_sense c X
    min_command = -2340.0  # MGM actuators.yaml:64
    max_command = 2340.0  # MGM actuators.yaml:64
    position = [0.0, 0.0, 0.0]  # MGM power_law/params.py default position
    direction = [1.0, 0.0, 0.0]  # MGM power_law/params.py default orientation_rpy 0
    power_law_coefficient = 0.00842703  # MGM actuators.yaml:65
    power_law_exponent = 2.2633  # MGM actuators.yaml:66
    astern_thrust_ratio = 1.0  # MGM actuators.yaml:67
    reaction_torque_per_thrust = 0.002852  # MGM actuators.yaml:68
