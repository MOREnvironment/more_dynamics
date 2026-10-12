from __future__ import annotations

# REMUS 100 fin: one lifting fin (lifting_fin with rigid-point inflow, small-angle flow angle, quadratic-drag
# section, first-order servo). Source of every value: the comment beside it, SOURCE.md gives its kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    fin_position = [-0.638, 0.0, 0.0]  # m, Prestero 2001, Table A.5, p. 103 (x_finpost)
    chord_axis = [1.0, 0.0, 0.0]  # unit chord direction, BODY FRD (convention)
    lift_axis = [0.0, 1.0, 0.0]  # unit positive-lift direction, BODY FRD (convention)
    fin_area = 0.00665  # m^2, Prestero 2001, Table A.5, p. 103 (S_fin)
    time_constant = 0.1  # s, Sarhadi 2026, Fig. 4, p. 4 (a demanding choice of that study, not a measurement)
    max_deflection = 0.3490658503988659  # rad = 20 deg, Sarhadi 2026, Fig. 4, p. 4
    max_rate = 0.5235987755982988  # rad/s = 30 deg/s, Sarhadi 2026, Fig. 4, p. 4
    lift_slope = 3.12  # 1/rad, Prestero 2001, Table A.5, p. 103 (c_L_alpha)
    zero_lift_drag = 0.0  # 1, zero with its reason: Prestero 2001, eq. 4.37, p. 31 carries no fin drag (read by linear_section only)
