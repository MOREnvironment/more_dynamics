from __future__ import annotations


class ComponentParameters:
    max_forward_thrust = 239.364
    max_speed = 3.0864
    time_constants = [1.0, 1.0]
    damping_ratios = [0.3, 0.2, 0.4]
    yaw_damping_nonlinearity = 10.0
    manoeuvring_damping = 'time_constants'
    sway_yaw_derivatives = [0.0, 0.0, 0.0, 0.0]
    modulus_coefficients = [0.0, 0.0]
    cross_flow = 'hoerner_strips'
    surge_resistance = 'none'
    time_constant = 100.0
    form_factor = 0.1
    crossover_speed = 2.0
    surge_added_mass_factor = 2.7
    ittc_reynolds_floor = 100000.0
