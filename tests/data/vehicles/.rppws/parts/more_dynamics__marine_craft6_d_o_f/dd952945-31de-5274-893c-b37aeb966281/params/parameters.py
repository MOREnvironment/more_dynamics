from __future__ import annotations


class ComponentParameters:
    integration_max_step = 0.05
    open_inputs = ["current.current_speed", "current.current_direction", "current.current_vertical_speed"]
    diagnostic_outputs = ["force_producers.1.generated_force", "force_producers.1.deflection", "force_producers.1.angle_of_attack"]
