from __future__ import annotations

# Marie site: Marie as MGM sets her up (more_generic_models vehicles/LAUV_Marie @ e86985cd). Source of every value: the comment beside it, SOURCE.md of the vehicle part has the table.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    site_form = 'latitude'  # MGM scenario.py:213-215
    latitude = 1.1073560310932362  # MGM globals.yaml:8
    water_density = 1025.0  # MGM globals.yaml:6
    kinematic_viscosity = 1e-06  # MGM dynamics/plant/matrices/hydrodynamics.py:221
