from __future__ import annotations

# Otter hull loads: damping calibrated from top speed and thrust, cross-flow (otter.m 195-240).
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    max_forward_thrust = 239.364  # otter.m, 24.4 kgf = 24.4 * 9.81 N
    max_speed = 3.0864  # otter.m, 6 knots = 6 * 0.5144 m/s
    time_constants = [1.0, 1.0]  # otter.m 99-101, 202-240
    damping_ratios = [0.3, 0.2, 0.4]  # otter.m 99-101, 202-240
    yaw_damping_nonlinearity = 10.0  # otter.m 99-101, 202-240
    surge_resistance = "none"
    time_constant = 100.0  # osv.m 128, vessel.T1 = 100 s (read with "ittc")
    form_factor = 0.1  # XuuITTC.m 33, k = 0.1; Fossen 2011, p. 125 (read with "ittc")
    crossover_speed = 2.0  # forceSurgeDamping.m 57 @ ac77394, u_cross = 2 m/s (read with "ittc")
    surge_added_mass_factor = 2.7  # addedMassSurge.m 33, A11 = 2.7 rho nabla^(5/3) / L^2 (read with "ittc")
    ittc_reynolds_floor = 100000.0  # XuuITTC.m 35, Re_min = 1e5 (read with "ittc")
