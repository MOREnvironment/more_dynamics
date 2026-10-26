from __future__ import annotations

# Grethe hull loads: linear coupled sway/yaw damping from the time constants (D_ii = M_ii / T_i,
# Y_r = N_v = 0 unidentified, E-101 Q2 a), ITTC + RA14 residual surge resistance (E-101 Q1 a).
# Source of every value: the comment beside it, from grethe_parameters.json (R10, E-104).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    manoeuvring_damping = "time_constants_coupled"  # Y_v = -M22 / T_sway and N_r = -M66 / T_yaw computed, Y_r and N_v given (grethe.m 361-369)
    time_constants = [8.0, 4.0]  # s: [T_sway, T_yaw] (grethe.m 251-252, estimates; T_yaw chosen with N_v for a yaw-stable open loop)
    sway_yaw_coupling = [0.0, 2500.0]  # [Y_r, N_v], N s/m and N m s/m (grethe.m 363-364: Y_r = 0 not identified; N_v = 2500 chosen, no source)
    cross_flow = "hoerner_strips"  # A-64 Section 3.4
    surge_resistance = "ittc_residual_ramped"  # ITTC friction plus the RA14 residual, its range edge softened by a smoothstep (grethe.m 386-388, 443-495; owner E-101 Q1 a)
    time_constant = 15.0  # estimate: A-43 M26 (= M20 T1)
    form_factor = 0.1  # published: Fossen 2011 p. 125; A-43 M27
    crossover_speed = 2.0  # published (MSS constant): A-43 M28
    surge_added_mass_factor = 2.7  # published: addedMassSurge.m 33; A-43
    ittc_reynolds_floor = 100000.0  # published + owner default: A-43 M29
    max_forward_thrust = 0  # zero on purpose: surge resistance is the ITTC + residual forms alone, a thrust here would count a second resistance (A-74 Q1, owner E-135); a coast-down-fitted resistance option follows
    max_speed = 2.4  # fitted from her 2026 control-station logs, through water, logged band 2.32-2.47 m/s (A-74 section 1, a74_top_speed.py); replace by X-9 G14 sea runs
    damping_ratios = [0.3, 0.05, 0.5]  # [heave, roll, pitch]: roll 0.05 from her logs, band 0.03-0.08 (A-74 section 3, owner E-135); heave and pitch approved guess 2026-10-12 (owner E-135), no source, until quay free decays (X-9 G12)
