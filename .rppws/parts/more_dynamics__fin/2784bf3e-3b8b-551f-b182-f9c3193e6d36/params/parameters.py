from __future__ import annotations

# REMUS 100 (estimated) rudder fin: one of the two rudder fins of the estimated REMUS 100 (the vehicle links it twice).
# Source of every value: the comment beside it, SOURCE.md gives its kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    fin_position = [-0.8, 0.0, 0.0]  # m, the fin post at the tail of the hull centred on the origin: x = -L/2 (derived)
    chord_axis = [1.0, 0.0, 0.0]  # unit chord direction, BODY FRD (convention)
    lift_axis = [0.0, -1.0, 0.0]  # unit positive-lift direction, BODY FRD: negative on y: a positive rudder turns to starboard
    fin_area = 0.00665  # m^2, Prestero 2001, Table 2.2, p. 17 (S_fin, one fin)
    servo_dynamics = 'none'  # the deflection is the command
    servo_rate_limit = False  # no rate limit
    servo_angle_limit = 'on_command'  # the command saturated at the deflection limit
    max_deflection = 0.3490658503988659  # rad, 20 deg
    inflow = 'rigid_point'  # velocity at the fin post, Prestero 2001, eq. 4.40, p. 32
    flow_angle = 'small_angle'  # Prestero 2001, eqs. 4.41-4.43, pp. 32-33
    interference = 'none'  # both factors one, Prestero 2001, eq. 4.41, p. 32
    section = 'lift_times_angle'  # C_L = a alpha, C_D = C_L alpha
    lift_law = 'finite_wing'  # the slope from the fin geometry
    lift_slope_form = 'low_aspect_ratio'  # Prestero 2001, eqs. 4.38-4.39, p. 31
    fin_span = 0.0857  # m, Prestero 2001, Table 2.2, p. 17 (b_fin)
    aspect_ratio_factor = 2.0  # AR_e = 2 b^2 / S, Prestero 2001, eq. 4.39, p. 31
    section_efficiency = 0.9  # the empirical factor a-bar, Prestero 2001, Table 2.2, p. 17 and eq. 4.38, p. 31
