"""The CasADi functions of the parts of a ``MarineCraft6DOF`` vehicle: one
function per physical form of each slot (site, current, hull form, rigid
body, added mass, its Coriolis term, restoring, hydrodynamic loads, force
producers). The plugin layer wraps each as a plugin class; the vehicle's
equation of motion is ``more_dynamics.models.vehicles.marine_craft_6dof``.

Author:    Enio Krizman
Date:      2026-10-08
"""
