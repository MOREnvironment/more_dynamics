from __future__ import annotations


class ComponentParameters:
    integration_max_step = 0.05
    open_inputs = ["current.current_speed", "current.current_direction", "current.current_vertical_speed"]
    diagnostic_outputs = ["mass_matrix", "relative_velocity"]
