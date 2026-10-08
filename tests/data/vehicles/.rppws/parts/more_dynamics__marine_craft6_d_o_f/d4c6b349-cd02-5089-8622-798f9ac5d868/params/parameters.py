from __future__ import annotations


class ComponentParameters:
    integration_max_step = 0.05
    open_inputs = []
    diagnostic_outputs = ["mass_matrix", "rigid_body_mass_matrix", "added_mass_matrix", "rigid_body_coriolis_matrix", "added_mass_coriolis_matrix", "restoring_matrix", "hydrodynamic_loads.0.damping_matrix", "hydrodynamic_loads.0.hydrodynamic_force", "hydrodynamic_loads.1.hydrodynamic_force"]
