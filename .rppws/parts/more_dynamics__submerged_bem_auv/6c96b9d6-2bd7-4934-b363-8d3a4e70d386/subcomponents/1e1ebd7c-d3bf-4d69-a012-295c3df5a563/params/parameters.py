from __future__ import annotations

# Marie hull loads: Marie as MGM sets her up (more_generic_models vehicles/LAUV_Marie @ e86985cd). Source of every value: the comment beside it, SOURCE.md of the vehicle part has the table.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    T_1236 = [50.0, 5.0, 5.0, 5.0]  # MGM plant.yaml:30
    delta_zeta_45 = [0.2, 0.2]  # MGM plant.yaml:31
    kappa_456 = [5.0, 5.0, 0.0]  # MGM plant.yaml:32
    surge_damping_blend = 3.0  # MGM plant.yaml:33
    p_factor = 0.7  # MGM plant.yaml:34
    c_d = 0.42  # MGM plant.yaml:35
    e_oswald = 0.3  # MGM plant.yaml:36
    cross_flow_model = 'cylinder'  # MGM plant.yaml:41
    smooth_speed = True  # plugin default, no MGM counterpart
    smooth_speed_epsilon = 1e-09  # plugin default, no MGM counterpart
