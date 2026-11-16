from __future__ import annotations

from pathlib import Path

import yaml

# rpp rewrites this file as plain values, without the comments and the yml read, when parameters are saved from the editor or by `parameter set`.
# Keel fin of Marie (MGM LAUV_Marie): the keel fin of MGM's LiftingFins, fin index 1 (MGM actuators.yaml:42-46). Fin parameters it shares with the other three
# fins are read from data/fins_common.yml of the vehicle part (part/subcomponents/<id>/params/parameters.py: four levels up); the table is in SOURCE.md.
with open(Path(__file__).resolve().parents[3] / "data" / "fins_common.yml", encoding="utf-8") as _file:
    _common = yaml.safe_load(_file)


class ComponentParameters:
    fin_position = [-0.83, 0.0, 0.0]  # MGM (fin_x, fin_radius cos(roll), fin_radius sin(roll)): actuators.yaml:46 (-0.83), :45 (fin_radius 0.0, placeholder M-0679), :43 (roll 90 deg); lifting_fins.py:92
    lift_axis = [0.0, -1.0, 0.0]  # MGM sign * R_x(roll) z, lifting_fins.py:78 and :179: roll 90 deg (actuators.yaml:43), deflection_sign +1 (:44) = (+1) (0, -sin 90, cos 90) = (+1)(0, -1, 0); command 0 (rudder_angle) = input of this fin (actuators.yaml:42, fin_command)
    chord_axis = _common["chord_axis"]
    fin_area = _common["fin_area"]
    servo_dynamics = _common["servo_dynamics"]
    servo_rate_limit = _common["servo_rate_limit"]
    servo_angle_limit = _common["servo_angle_limit"]
    max_deflection = _common["max_deflection"]
    inflow = _common["inflow"]
    flow_angle = _common["flow_angle"]
    interference = _common["interference"]
    section = _common["section"]
    lift_law = _common["lift_law"]
    lift_slope_form = _common["lift_slope_form"]
    fin_span = _common["fin_span"]
    aspect_ratio_factor = _common["aspect_ratio_factor"]
    section_efficiency = _common["section_efficiency"]
    zero_lift_drag = _common["zero_lift_drag"]
