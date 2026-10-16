from __future__ import annotations


class ComponentParameters:
    integration_max_step = 0.05  # numerical setting (plugin default)
    length = 5.3  # npsauv.m:105 L = 5.3 (via MSS, Healey & Lienard 1993 not read)
    body_mass = 5443.425076452599  # npsauv.m:109 mass = W / g (derived from the given weight) (derived, via MSS)
    inertia_diagonal = [2038, 13587, 13587]  # npsauv.m:110 [Ix Iy Iz] (via MSS, Healey & Lienard 1993 not read)
    inertia_products = [-13.58, -13.58, -13.58]  # npsauv.m:111 [Ixy Iyz Ixz] (full tensor, not diagonal) (via MSS, Healey & Lienard 1993 not read)
    body_center_of_gravity = [0, 0, 0.061]  # npsauv.m:106 [xG yG zG] (via MSS, Healey & Lienard 1993 not read)
    weight = 53400  # npsauv.m:108 W = 53400 (via MSS, Healey & Lienard 1993 not read)
    buoyancy = 53400  # npsauv.m:108 B = 53400 (neutral, W = B as given) (via MSS, Healey & Lienard 1993 not read)
    center_of_buoyancy = [0, 0, 0]  # npsauv.m:107 [xB yB zB] (via MSS, Healey & Lienard 1993 not read)
    derivative_scaling = 'nondimensional'  # npsauv.m 116-124: the published nondimensional form (Fossen 2011, App. D.2 prime-scaling)
    current_form = 'yaw_rate_terms'  # npsauv.m 83-90: the MSS form of the current
    nondim_Xudot = -0.0076  # npsauv.m: nondimensional hydrodynamic derivative Xudot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Yvdot = -0.055  # npsauv.m: nondimensional hydrodynamic derivative Yvdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Ypdot = 0.00012  # npsauv.m: nondimensional hydrodynamic derivative Ypdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Yrdot = 0.0012  # npsauv.m: nondimensional hydrodynamic derivative Yrdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Zwdot = -0.24  # npsauv.m: nondimensional hydrodynamic derivative Zwdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Zqdot = -0.0068  # npsauv.m: nondimensional hydrodynamic derivative Zqdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Kvdot = 0.00012  # npsauv.m: nondimensional hydrodynamic derivative Kvdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Kpdot = -0.001  # npsauv.m: nondimensional hydrodynamic derivative Kpdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Krdot = -3.4e-05  # npsauv.m: nondimensional hydrodynamic derivative Krdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Mwdot = -0.0068  # npsauv.m: nondimensional hydrodynamic derivative Mwdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Mqdot = -0.017  # npsauv.m: nondimensional hydrodynamic derivative Mqdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Nvdot = 0.0012  # npsauv.m: nondimensional hydrodynamic derivative Nvdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Npdot = -3.4e-05  # npsauv.m: nondimensional hydrodynamic derivative Npdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
    nondim_Nrdot = -0.0034  # npsauv.m: nondimensional hydrodynamic derivative Nrdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read)
