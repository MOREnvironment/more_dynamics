from __future__ import annotations

# Grethe hull loads: linear coupled sway/yaw damping from the time constants (D_ii = M_ii / T_i,
# Y_r = N_v = 0 unidentified, E-101 Q2 a), ITTC + RA14 residual surge resistance (E-101 Q1 a).
# Source of every value: the comment beside it, from grethe_parameters.json (R10, E-104).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    manoeuvring_damping = "linear_coupled"  # owner decision: A-64 Section 3.4 G0 set
    sway_yaw_derivatives = [-251.875, 0.0, 0.0, -229.86112500000002]  # [Y_v, Y_r, N_v, N_r]: Y_v = -M22/T_sway, N_r = -M66/T_yaw (time_constants [8, 16] s); Y_r = N_v = 0, unidentified (E-101 Q2 a)
    cross_flow = "hoerner_strips"  # A-64 Section 3.4
    surge_resistance = "ittc_residual"  # owner decision: E-101 Q1 a
    time_constant = 15.0  # estimate: A-43 M26 (= M20 T1)
    form_factor = 0.1  # published: Fossen 2011 p. 125; A-43 M27
    crossover_speed = 2.0  # published (MSS constant): A-43 M28
    surge_added_mass_factor = 2.7  # published: addedMassSurge.m 33; A-43
    ittc_reynolds_floor = 100000.0  # published + owner default: A-43 M29
