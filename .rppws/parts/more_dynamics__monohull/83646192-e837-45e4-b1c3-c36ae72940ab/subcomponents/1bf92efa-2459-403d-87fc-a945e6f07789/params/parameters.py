from __future__ import annotations

# Grethe restoring: one hull at the equilibrium draft, CO at the waterline midpoint
# (E-101 Q4 a). Source of every value: the comment beside it, from grethe_parameters.json
# (R10, E-104). rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    hull_count = 1  # convention: A-43 M15, monohull
    longitudinal_inertia_factor = 1.0  # estimate: A-43 M17
    longitudinal_center_of_flotation = -0.3  # estimate: A-43 M18
    reference_point = [-2.24, 0.0, 0.0]  # owner decision: E-101 Q4 a, CO at the midpoint of the waterline length
    wetted_surface_method = "regression_table"  # owner decision: A-64 Section 3.3 R11
