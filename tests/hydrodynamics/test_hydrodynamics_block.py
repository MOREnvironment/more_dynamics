"""Gate tests for the hydrodynamics blocks: linear damping, cross-flow drag, lift/drag.

Written before the blocks exist (job A-22, 2026-10-05); the blocks are ported
by job A-24. Every block test fails today with "block not ported yet".

Rewired by job A-30 (2026-10-05): section "G1-MSS" compares the default paths
with MATLAB running current MSS (job A-26, ``*_mss_current.csv``); where the
default deviates by register (``30_checks/README.md`` D-MSS-1, Reynolds number
on the diameter) the MATLAB file is compared with the MSS flag
(``cross_flow_reynolds_length="length"``) and the default with the
transcription plus the register's stated effect. D-MSS-2 (owner E-24): Hoerner
below its table raises ``ValueError``. The test-side transcriptions stay as the
second check and are themselves checked against the new files.

Completed by job U3a (verifier, 2026-10-06; ledger
``agents-more/30_checks/2026-10-06_U3_hydrodynamics_gates.md``):

* MSS ``ac77394`` (issue #81) puts the cross-flow Reynolds number on the
  diameter, so D-MSS-1 is resolved: the default path is G1 against MATLAB at
  ``ac77394`` (``spheroid_matlab_reference_mss_ac77394.csv``) and the
  ``cross_flow_reynolds_length`` flag goes (owner, 2026-10-06, E-20 "a"). The
  block tests that drove that flag are removed here; ``test_D_MSS_1_flag_is_gone``
  fails until the porter (U3b) removes it. The ``99bf0b3`` file stays for history.
* Every U3 row has a test: ``get_D_linear_surface_vessel`` is ``Dmtrx.m``'s
  surface-craft branch (a third function in ``linear_damping``), and
  ``force_surge_damping`` is ``forceSurgeDamping.m`` (new module
  ``surge_damping``); both G1 against MATLAB calling those MSS functions.
* Physical-sign tests (owner E-23: MSS is the reference, not the truth): the
  linear surge coefficient reproduces its time constant (MSS right, the numpy
  source wrong); MSS's surge blend is not dissipative in reverse and its ITTC
  branch blows up at Rn = 100 (documented, owner decision pending); cross-flow
  pitch/yaw terms equal the first-principles strip integral.
* A-27 finding 4: the Jacobian entries that are non-finite today are pinned.
* E-25: no local transform in our hydrodynamics modules.

Owner rulings these tests encode
--------------------------------
* E-18 / E-20 Q2 a: the **default** path equals current MSS; each template
  departure stays behind a named flag that reproduces the old source value.
* E-10 a: cross-flow drag on MSS's 20 strip midpoints with the full-precision
  ``cylinderDrag.m`` table (default); the 21-end-point grid behind a flag.
* E-20 Q3 b: the cross-flow Reynolds number uses the **diameter**, which is
  newest MSS since ``ac77394`` (``cylinderDrag.m`` 79-80); no flag (E-20, "a",
  2026-10-06).
* E-16: generic by construction (no vehicle name or number inside a block).
* E-28: Luka's ``linear_surface.py`` and its tests are not touched; the U3
  template audit of it is a reading, in the ledger.

Contract the porter must provide
--------------------------------
``more_dynamics.models.hydrodynamics.linear_damping``:

* ``preprocess_submerged_linear_damping(rigid_body_mass_matrix,
  added_mass_matrix, weight, center_of_gravity, center_of_buoyancy,
  time_constants, damping_ratios, sway_damping_fade=False)``. Matrices 6x6,
  weight N, centres (3,) m from the CO, ``time_constants = (T1, T2, T6)`` s,
  ``damping_ratios = (zeta4, zeta5)``. Constants expose
  ``damping_coefficients`` (6,): the ``Dmtrx.m`` diagonal before any fade.
* ``submerged_linear_damping_casadi(constants)``: input ``["nu_r"]`` (6x1),
  outputs ``["D", "tau"]``: ``D`` (6x6, positive diagonal as ``Dmtrx.m``) with
  ``D[0,0] *= exp(-3 U_r)`` (``remus100.m:216``), ``U_r = |nu_r[0:3]|``, and
  ``D[1,1]`` faded the same way only when ``sway_damping_fade`` is True (the
  template departure, A-18); ``tau = -D @ nu_r``.
* ``preprocess_surface_linear_damping(mass_matrix, restoring_matrix,
  max_forward_thrust, max_speed, time_constants, damping_ratios,
  yaw_damping_nonlinearity)``. ``mass_matrix`` = M_RB + M_A (6x6);
  ``restoring_matrix`` = G about the centre of flotation (6x6, only G33, G44,
  G55 are used); ``max_forward_thrust`` N; ``max_speed`` m/s;
  ``time_constants = (T_sway, T_yaw)``; ``damping_ratios = (zeta3, zeta4,
  zeta5)``; ``yaw_damping_nonlinearity`` = the 10 of ``otter.m:240``.
  Constants expose ``damping_derivatives`` (6,) = ``[Xu Yv Zw Kp Mq Nr]``,
  the negative numbers of ``otter.m`` 202-207.
* ``surface_linear_damping_casadi(constants)``: input ``["nu_r"]``, outputs
  ``["D", "tau"]``: ``D = -diag(damping_derivatives)`` (positive, the same
  sign as the submerged ``D``) and ``tau`` = ``otter.m`` 235-242 (the force
  on the vehicle, quadratic yaw term included).
* ``preprocess_floating_linear_damping(rigid_body_mass_matrix,
  added_mass_matrix, restoring_matrix, time_constants, damping_ratios)``
  (U3a): ``Dmtrx.m`` 51-63, the surface-craft branch, which is the numpy
  ``get_D_linear_surface_vessel`` with all six ``coeff`` = 1. Matrices 6x6
  (``restoring_matrix`` = G, only G33, G44, G55 used), ``time_constants =
  (T1, T2, T6)`` s, ``damping_ratios = (zeta3, zeta4, zeta5)``; MSS fixes
  ``zeta3 = 0.2`` (``Dmtrx.m`` 57), so it is a parameter like the densities
  and the MSS tests pass 0.2. The source's ``coeff`` multipliers are not in
  the contract: ``coeff_i`` scales ``D_ii``, which a parameter set gets by
  scaling ``T_i`` or ``zeta_i`` (G2 checks that mapping). Constants expose
  ``damping_coefficients`` (6,).
* ``floating_linear_damping_casadi(constants)``: input ``["nu_r"]``, outputs
  ``["D", "tau"]``: ``D`` = that diagonal (positive, ``Dmtrx.m`` sign),
  ``tau = -D @ nu_r``.

``more_dynamics.models.hydrodynamics.surge_damping`` (U3a; ``forceSurgeDamping.m``):

* ``preprocess_surge_damping_ittc(mass, length, wetted_surface,
  water_density, time_constant, form_factor, crossover_speed)``: the
  eight-argument call (``forceSurgeDamping.m`` 68-75), ``Xuu`` from the
  ITTC-1957 line at the current speed. ``form_factor`` = the ``k = 0.1`` of
  line 70, ``crossover_speed`` = the ``u_cross = 2`` of line 57; MSS
  hard-codes both, they are parameters here (E-16) and the tests pass MSS's.
* ``preprocess_surge_damping_max_thrust(mass, length, water_density,
  time_constant, max_speed, max_thrust, crossover_speed)``: the nine-argument
  call (lines 64-66), ``Xuu = -max_thrust / max_speed^2``.
* Both: ``A11 = 2.7 rho (m/rho)^(5/3) / L^2`` (``addedMassSurge.m`` 33-34),
  ``Xudot = -A11`` (line 60), ``Xu = -(m - Xudot) / T1 = -(m + A11) / T1``
  (line 61). Constants expose ``added_mass`` (A11 > 0), ``linear_coefficient``
  (Xu < 0) and, for ``max_thrust``, ``quadratic_coefficient`` (Xuu < 0).
* ``surge_damping_casadi(constants)``: input ``["nu_r"]`` (6x1), output
  ``["tau"]`` = ``[X 0 0 0 0 0]`` with ``sigma = 1 - tanh(u_r / u_cross)``
  and ``X = sigma Xu u_r + (1 - sigma) Xuu |u_r| u_r`` (lines 79-82),
  ``u_r = nu_r[0]``; the force on the vehicle (added, as ``osv.m`` 180-182).
  The vehicle that uses it must set its linear ``D[0,0] = 0`` (``osv.m`` 185).

``more_dynamics.models.hydrodynamics.cross_flow``:

* ``preprocess_cross_flow_drag(length, beam, draft, water_density,
  drag_model, strip_grid="midpoint")``; ``drag_model`` in {"cylinder",
  "hoerner"} (no default), ``strip_grid`` in {"midpoint", "endpoint"}. The
  Reynolds number is on the ``beam`` argument, the diameter (MSS calls
  ``crossFlowDrag(L, D, D, ...)``; ``cylinderDrag.m`` 79-80 at ``ac77394``);
  ``cross_flow_reynolds_length`` is gone (owner 2026-10-06).
* ``cross_flow_drag_casadi(constants)``: input ``["nu_r"]``, output ``["tau"]``.

``more_dynamics.models.hydrodynamics.lift_drag``:

* ``preprocess_lift_drag(span, planform_area, parasitic_drag_coefficient,
  oswald_efficiency, water_density)``.
* ``lift_drag_casadi(constants)``: input ``["nu_r"]``, output ``["tau"]``;
  ``alpha = atan2(w_r, u_r)``, ``U_r = |nu_r[0:3]|`` (``remus100.m`` 127-128).

Map rows covered, numpy source read in full (``more_generic_models/more_generic_models/``)
---------------------------------------------------------------------------------------
* ``dynamics/plant/matrices/linear_damping.py``: ``get_damping_coeff`` 16-84,
  ``D_linear`` 88-112 (fades surge **and sway**), ``D_linear_catamaran``
  165-197, ``get_tau_damping`` 200-226, ``get_D_linear_surface_vessel``
  115-160 (U3a; ``Dmtrx.m`` surface branch, not Luka's module, E-28),
  ``force_surge_damping`` 230-294 (U3a; on ``ASVHull.get_tau_drag``'s path,
  ``asv_hull.py`` 185-197, with ``D[0,0] = 0`` at 145). Its ``Xudot`` is
  ``+A11`` (``added_mass.py`` 338-357 returns A11), so its ``Xu`` is
  ``-(m - A11)/T1``: a sign error against ``forceSurgeDamping.m`` 60-61 and
  Fossen (2011) eq. 6.76; G2 runs the source with that one sign corrected.
* ``dynamics/plant/matrices/hydrodynamics.py``: ``get_lift_drag_coeff``
  62-75, ``force_lift_drag`` 78-88, ``_hoerner`` 95-100, ``_cylinder_drag``
  103-120 (Re column rounded to 4 digits, A-2), ``cross_flow_drag`` 127-178,
  ``_get_strips`` 181-190 (21 end points).
* ``dynamics/plant/matrices/drag_models.py``: ``hoerner`` 16-64,
  ``cylinder_drag`` 68-162 (same rounded table; not on any vehicle's path).
* ``dynamics/plant/auv_spheroid/auv_spheroid.py``: ``compute_constant_values``
  157-197, ``get_D`` 235-236, ``get_tau_lift_drag`` 269-294,
  ``get_tau_cross_flow`` 296-323.
* ``dynamics/plant/asv_catamaran/asv_catamaran.py``: ``compute_constant_values``
  121-186, ``get_tau_damping`` 232-233, ``get_tau_cross_flow`` 235-244.

Conventions found (hidden assumptions, see the ledger)
------------------------------------------------------
* Two damping sign conventions: ``Dmtrx.m`` positive D subtracted
  (``- D*nu_r``); ``otter.m`` negative derivatives added (``+ tau_damp``).
  The contract outputs both as a positive ``D`` and a force ``tau``.
* The cross-flow coefficient is one number per call, set by the mid-body
  cross-flow speed ``sqrt(v_r^2 + w_r^2)``; kinematic viscosity 1e-6 m^2/s
  is inside ``Re = U * L * 1e6``; strip height = ``draft``.
* ``crossFlowDrag.m`` hard-codes rho = 1025, ``forceLiftDrag.m`` rho = 1026
  and ``coeffLiftDrag.m`` e = 0.3; the blocks take both as parameters and the
  MSS tests pass the hard-coded values.

Frozen reference: ``tests/data/hydrodynamics/`` (``SOURCE.md``).
"""

import ast
import importlib
import inspect
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "hydrodynamics"
# Outside this repo only through environment variables (agents-more rule 9):
# MSS_DIR = an MSS checkout, MORE_GENERIC_MODELS_DIR = the more_generic_models
# repository root. Default: not set. The gates read the frozen CSVs and the
# snapshot of every cited MSS / generator line (``cited_lines_snapshot.json``,
# made by ``snapshot_cited_lines.py``); only the live re-read of those lines
# and the numpy-source tests need the variables, and they skip without them.
CITED_SNAPSHOT = DATA_DIR / "cited_lines_snapshot.json"
LINEAR_DAMPING = "more_dynamics.models.hydrodynamics.linear_damping"
CROSS_FLOW = "more_dynamics.models.hydrodynamics.cross_flow"
LIFT_DRAG = "more_dynamics.models.hydrodynamics.lift_drag"
SURGE_DAMPING = "more_dynamics.models.hydrodynamics.surge_damping"

G1_TOLERANCE = 1e-9   # 30_checks/README.md, gate G1
G2_TOLERANCE = 1e-10  # 30_checks/README.md, gate G2
G4_FACTOR = 10.0      # 30_checks/README.md, gate G4
SEED = 20261005
N_RANDOM_STATES = 1000
# A-2 ledger: the rounded Re table moves the spheroid cross-flow force by 2.3e-2 N
ROUNDED_TABLE_RESIDUAL_BOUND = 0.03

VEHICLE_NAMES = ("remus", "otter", "grethe", "marie", "hugin", "lauv",
                 "mariner", "torqeedo", "cybership", "prestero")

REMUS = "$MSS_DIR/CRAFT/AUV/models/remus100.m"
OTTER = "$MSS_DIR/CRAFT/USV/models/otter.m"
DMTRX = "$MSS_DIR/LIBRARY/modeling/Dmtrx.m"
FLD = "$MSS_DIR/LIBRARY/modeling/forceLiftDrag.m"
CLD = "$MSS_DIR/LIBRARY/modeling/coeffLiftDrag.m"
CFD = "$MSS_DIR/LIBRARY/modeling/crossFlowDrag.m"
CYL = "$MSS_DIR/HYDRO/cylinderDrag.m"
HOERNER = "$MSS_DIR/HYDRO/Hoerner.m"
GRAVITY = "$MSS_DIR/INS/functions/gravity.m"
SPHEROID_M = "$MSS_DIR/LIBRARY/modeling/spheroid.m"
FSD = "$MSS_DIR/LIBRARY/modeling/forceSurgeDamping.m"
AMS = "$MSS_DIR/LIBRARY/modeling/addedMassSurge.m"
XITTC = "$MSS_DIR/LIBRARY/modeling/XuuITTC.m"
OSV = "$MSS_DIR/CRAFT/SHIP/models/osv.m"
GEN_SPH = "$MORE_GENERIC_MODELS_DIR/more_generic_models/test/plant/auv_spheroid/test_mss_matlab/test_dynamics_consitency.m"
GEN_CAT = "$MORE_GENERIC_MODELS_DIR/more_generic_models/test/plant/asv_catamaran/test_mss_matlab/test_dynamics_consistency.m"

# (file relative to <more>, line) -> how the stripped line starts. Pinned by
# test_cited_lines_are_unchanged; every constant below is parsed from these.
CITED_LINES = {
    (REMUS, 97): "mu = deg2rad(63.446827);",
    (REMUS, 127): "alpha = atan2( nu_r(3), nu_r(1) );",
    (REMUS, 128): "U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );",
    (REMUS, 132): "L_auv = 1.6;",
    (REMUS, 133): "D_auv = 0.19;",
    (REMUS, 134): "S = 0.7 * L_auv * D_auv;",
    (REMUS, 136): "b = 1.0096 * D_auv/2;",
    (REMUS, 138): "r_bG = [ 0 0 0.02 ]';",
    (REMUS, 139): "r_bB = [ 0 0 0 ]';",
    (REMUS, 144): "Cd = 0.42;",
    (REMUS, 145): "CD_0 = Cd * pi * b^2 / S;",
    (REMUS, 193): "T1 = 20;",
    (REMUS, 194): "T2 = 20;",
    (REMUS, 195): "zeta4 = 0.3;",
    (REMUS, 196): "zeta5 = 0.8;",
    (REMUS, 197): "T6 = 1;",
    (REMUS, 215): "D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bG' r_bB']);",
    (REMUS, 216): "D(1,1) = D(1,1) * exp(-3 * U_r);",
    (REMUS, 218): "tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);",
    (REMUS, 219): "tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r,'cylinder');",
    (DMTRX, 30): "M = MRB + MA;",
    (DMTRX, 44): "T3 = T2;",
    (DMTRX, 45): "w4 = sqrt( W * (r_bg(3)-r_bb(3)) / M(4,4) );",
    (DMTRX, 46): "w5 = sqrt( W * (r_bg(3)-r_bb(3)) / M(5,5) );",
    (DMTRX, 48): "D = diag( [M(1,1)/T1 M(2,2)/T2 M(3,3)/T3...",
    (DMTRX, 49): "M(4,4)*2*zeta4*w4  M(5,5)*2*zeta5*w5 M(6,6)/T6 ] );",
    (FLD, 26): "rho = 1026;",
    (FLD, 28): "[CL,CD] = coeffLiftDrag(b,S,CD_0,alpha,0);",
    (FLD, 30): "F_drag = 1/2 * rho * U_r^2 * S * CD;",
    (FLD, 31): "F_lift = 1/2 * rho * U_r^2 * S * CL;",
    (FLD, 35): "cos(alpha) * (-F_drag) - sin(alpha) * (-F_lift)",
    (FLD, 37): "sin(alpha) * (-F_drag) + cos(alpha) * (-F_lift)",
    (CLD, 54): "e = 0.3;",
    (CLD, 57): "AR = b^2/S;",
    (CLD, 60): "CL_alpha = pi * AR / ( 1 + sqrt(1 + (AR/2)^2) );",
    (CLD, 61): "CL_linear = CL_alpha * alpha;",
    (CLD, 68): "CD = CD_0 + CL.^2 / (pi * e * AR);",
    (CFD, 24): "rho = 1025;",
    (CFD, 25): "nStrips = 20;",
    (CFD, 26): "dx = L/nStrips;",
    (CFD, 41): "xL = -L/2 + (i - 0.5) * dx;",
    (CFD, 46): "U_h = abs(v_r + xL * r) * (v_r + xL * r);",
    (CFD, 47): "U_v = abs(w_r + xL * q) * (w_r + xL * q);",
    (CFD, 48): "Yh = Yh - 0.5 * rho * T * Cd_2D * U_h * dx;",
    (CFD, 49): "Zh = Zh - 0.5 * rho * T * Cd_2D * U_v * dx;",
    (CFD, 50): "Mh = Mh - 0.5 * rho * T * Cd_2D * xL * U_v * dx;",
    (CFD, 51): "Nh = Nh - 0.5 * rho * T * Cd_2D * xL * U_h * dx;",
    # cylinderDrag.m at ac77394 (issue #81): the Reynolds number on the diameter
    (CYL, 78): "U_crossflow = sqrt(nu_r(2)^2+nu_r(3)^2);",
    (CYL, 79): "nu_water = 1e-6;",
    (CYL, 80): "Re = U_crossflow * B / nu_water;",
    (CYL, 92): "if Re < 2e5",
    (DMTRX, 57): "zeta3 = 0.2;",
    (DMTRX, 58): "w3 = sqrt( G33 / M(3,3) );",
    (DMTRX, 62): "D = diag( [M(1,1)/T1 M(2,2)/T2 M(3,3)*2*zeta3*w3...",
    (DMTRX, 63): "M(4,4)*2*zeta4*w4  M(5,5)*2*zeta5*w5 M(6,6)/T6 ] );",
    (FSD, 57): "u_cross = 2;",
    (FSD, 60): "Xudot = -addedMassSurge(m,L,rho);",
    (FSD, 61): "Xu = -(m - Xudot) / T1;",
    (FSD, 64): "if (nargin == 9)",
    (FSD, 66): "Xuu = - thrust_max / u_max^2;",
    (FSD, 69): "nu_kin = 1e-6;",
    (FSD, 70): "k = 0.1;",
    (FSD, 71): "eps = 1e-10;",
    (FSD, 72): "Rn = (L / nu_kin) * abs(u_r);",
    (FSD, 73): "Cf = 0.075 / (log10(Rn + eps) -  2)^2;",
    (FSD, 74): "Xuu = -0.5 * rho * S * (1+k) * Cf;",
    (FSD, 79): "sigma = 1 - tanh(u_r / u_cross);",
    (FSD, 82): "X = sigma .* Xu .* u_r + (1 - sigma) .* Xuu .* abs(u_r) .* u_r;",
    (AMS, 33): "nabla = m / rho;",
    (AMS, 34): "A11 = 2.7 * rho * nabla^(5/3) / L^2;",
    (XITTC, 29): "C_B = 0.65;",
    (XITTC, 32): "nu_kin = 1e-6;",
    (XITTC, 33): "k = 0.1;",
    (XITTC, 34): "Re_min = 1e5;",
    (XITTC, 35): "Re = max(L * abs(u_r) / nu_kin,Re_min);",
    (XITTC, 36): "Cf = 0.075 / (log10(Re) - 2)^2;",
    (XITTC, 38): "S = 1.025 * L * (C_B*B + 1.7*T);",
    (XITTC, 39): "Xuu = -0.5 * rho * S * (1+k) * Cf;",
    (OSV, 62): "vessel.L = 83;",
    (OSV, 63): "vessel.B = 18;",
    (OSV, 64): "vessel.T = 5;",
    (OSV, 65): "vessel.rho = 1025;",
    (OSV, 66): "vessel.Cb = 0.65;",
    (OSV, 67): "vessel.S = vessel.L * vessel.B + 2 * vessel.T * vessel.B;",
    (OSV, 70): "vessel.K_max = [300e3 300e3 420e3 655e3]';",
    (OSV, 76): "vessel.thrust_max = vessel.K_max(3)+vessel.K_max(4);",
    (OSV, 77): "vessel.U_max = 7.7;",
    (OSV, 79): "vessel.nabla = vessel.Cb * vessel.L * vessel.B * vessel.T;",
    (OSV, 80): "vessel.m = vessel.rho * vessel.nabla;",
    (OSV, 129): "vessel.T1 = 100;",
    (OSV, 180): "[X,Xuu,Xu] = forceSurgeDamping(flag,nu_r(1),vessel.m,vessel.S,vessel.L, ...",
    (OSV, 181): "vessel.T1,vessel.rho,vessel.U_max,vessel.thrust_max);",
    (OSV, 185): "vessel.D(1,1) = 0;",
    (HOERNER, 47): "if B/(2*T) <= 4.00309",
    (HOERNER, 48): "CY_2D = interp1(CD_DATA(:,1),CD_DATA(:,2),B/(2*T));",
    (HOERNER, 50): "CY_2D = 0.559315;",
    (GRAVITY, 11): "g = 9.7803253359 * ( 1 + 0.001931850400 * sin(mu)^2 ) /...",
    (GRAVITY, 12): "sqrt( 1 - 0.006694384442 * sin(mu)^2 );",
    (OTTER, 90): "g   = 9.81;",
    (OTTER, 91): "rho = 1025;",
    (OTTER, 92): "L = 2.0;",
    (OTTER, 94): "m = 55.0;",
    (OTTER, 99): "T_sway = 1;",
    (OTTER, 100): "T_yaw = 1;",
    (OTTER, 101): "Umax = 6 * 0.5144;",
    (OTTER, 104): "B_pont  = 0.25;",
    (OTTER, 107): "Cb_pont = 0.4;",
    (OTTER, 121): "nabla = (m+mp)/rho;",
    (OTTER, 122): "T = nabla / (2 * Cb_pont * B_pont * L);",
    (OTTER, 192): "LCF = -0.2;",
    (OTTER, 197): "w3 = sqrt( G33/M(3,3) );",
    (OTTER, 198): "w4 = sqrt( G44/M(4,4) );",
    (OTTER, 199): "w5 = sqrt( G55/M(5,5) );",
    (OTTER, 202): "Xu = -24.4 * g / Umax;",
    (OTTER, 203): "Yv = -M(2,2) / T_sway;",
    (OTTER, 204): "Zw = -2 * 0.3 * w3 * M(3,3);",
    (OTTER, 205): "Kp = -2 * 0.2 * w4 * M(4,4);",
    (OTTER, 206): "Mq = -2 * 0.4 * w5 * M(5,5);",
    (OTTER, 207): "Nr = -M(6,6) / T_yaw;",
    (OTTER, 235): "Xh = Xu * nu_r(1);",
    (OTTER, 236): "Yh = Yv * nu_r(2);",
    (OTTER, 237): "Zh = Zw * nu_r(3);",
    (OTTER, 238): "Kh = Kp * nu_r(4);",
    (OTTER, 239): "Mh = Mq * nu_r(5);",
    (OTTER, 240): "Nh = Nr * (1 + 10 * abs(nu_r(6))) * nu_r(6);",
    (OTTER, 245): "tau_crossflow = crossFlowDrag(L,B_pont,T,nu_r);",
    (GEN_SPH, 206): "rho = 1025;",
    (GEN_SPH, 267): "D(1,1) = D(1,1) * exp(-3*U_r);",
    (GEN_SPH, 268): "D(2,2) = D(2,2) * exp(-3*U_r);",
    (GEN_SPH, 467): "rho = 1025;",
    (GEN_CAT, 12): "mp = 25;",
    (GEN_CAT, 136): "mu = deg2rad(63.446827);",
    (GEN_CAT, 137): "g = gravity(mu);",
    (GEN_CAT, 187): "Xu = -24.4*g/Umax; Yv = -M(2,2)/T_sway; Zw = -2*0.3*w3*M(3,3);",
}


# --------------------------------------------------------------------------
# Helpers: pinned text, contract, data
# --------------------------------------------------------------------------
def _env_dir(variable):
    value = os.environ.get(variable)
    if not value:
        pytest.skip(f"{variable} is not set (agents-more rule 9); set it to run this check")
    path = Path(value).expanduser()
    if not path.is_dir():
        pytest.skip(f"{variable}={value} is not a directory")
    return path


def _live_path(ref):
    """``$VARIABLE/relative/path`` -> a file under that environment variable."""
    variable, rel = ref[1:].split("/", 1)
    path = _env_dir(variable) / rel
    if not path.exists():
        pytest.skip(f"{rel} not found under {variable}")
    return path


def _snapshot():
    return json.loads(CITED_SNAPSHOT.read_text())


def _line(rel, number):
    """A cited line, stripped, from the frozen snapshot; it must be pinned."""
    assert (rel, number) in CITED_LINES, ("cite the line first", rel, number)
    line = _snapshot()["lines"][rel][str(number)]
    assert line.startswith(CITED_LINES[(rel, number)]), (rel, number, line)
    return line


def _live_line(rel, number):
    return _live_path(rel).read_text().splitlines()[number - 1].strip()


def _value(rel, number, env=None):
    """Right-hand side of ``name = expr;`` on a pinned line, evaluated."""
    expr = _line(rel, number).split("=", 1)[1].split(";")[0].strip()
    if expr.startswith("["):
        return np.array([float(t) for t in expr.strip("[]' ").split()])
    names = {"__builtins__": {}, "deg2rad": np.deg2rad, "pi": np.pi, "sqrt": np.sqrt}
    return float(eval(expr.replace("^", "**"), names, dict(env or {})))  # pinned arithmetic


def _number_in(rel, number, pattern):
    match = re.search(pattern, _line(rel, number))
    assert match, (rel, number, pattern)
    return float(match.group(1))


TABLES = ((CYL, "CD_DATA"), (CYL, "KAPPA_SUBCRITICAL_DATA"),
          (CYL, "KAPPA_SUPERCRITICAL_DATA"), (HOERNER, "CD_DATA"))


def _parse_table(text, name):
    """A ``NAME = [ ... ];`` numeric table of an MSS file, parsed from the text."""
    block = re.split(rf"(?m)^{name} = \[", text)[1].split("];")[0]
    rows = [[float(t) for t in ln.replace("...", "").split()] for ln in block.splitlines()]
    return np.array([r for r in rows if r])


def _table(rel, name):
    """A cited MSS table from the frozen snapshot."""
    assert (rel, name) in TABLES, (rel, name)
    return np.array(_snapshot()["tables"][rel][name])


def _gravity(mu):
    """INS/functions/gravity.m lines 11-12."""
    g0 = _number_in(GRAVITY, 11, r"g = ([0-9.]+) \*")
    k1 = _number_in(GRAVITY, 11, r"1 \+ ([0-9.]+) \*")
    e2 = _number_in(GRAVITY, 12, r"1 - ([0-9.]+) \*")
    return g0 * (1 + k1 * np.sin(mu) ** 2) / np.sqrt(1 - e2 * np.sin(mu) ** 2)


def _contract(name):
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed (owner decision E-7): {exc}")
        pytest.fail(f"block not ported yet: {exc}")


def _evaluate(function, nu_r):
    out = function(nu_r=np.asarray(nu_r, float))
    return {k: np.array(v, dtype=float) for k, v in out.items()}


def _load_csv(name):
    path = DATA_DIR / name
    header = path.read_text().splitlines()[0].split(",")
    return header, np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def _columns(header, values, prefix, n):
    return values[:, [header.index(f"{prefix}_{i:02d}") for i in range(1, n + 1)]]


LEGACY_SPHEROID = "spheroid_matlab_reference.csv"
LEGACY_CATAMARAN = "catamaran_matlab_reference.csv"
# MATLAB running current MSS (job A-26, SOURCE.md): the default-path references
MSS_SPHEROID = "spheroid_matlab_reference_mss_current.csv"      # MSS 99bf0b3
MSS_CATAMARAN = "catamaran_matlab_reference_mss_current.csv"
# job U3a: the A-26 generator re-run unchanged at MSS ac77394; only the four
# non-zero tau_crossflow columns differ from MSS_SPHEROID (test below)
MSS_AC_SPHEROID = "spheroid_matlab_reference_mss_ac77394.csv"
MSS_SURGE = "surge_damping_mss_ac77394.csv"          # forceSurgeDamping.m, direct calls
MSS_XUU_ITTC = "xuu_ittc_mss_ac77394.csv"            # XuuITTC.m, direct calls
MSS_FLOATING = "floating_damping_mss_ac77394.csv"    # Dmtrx.m surface branch, direct calls


def _spheroid_reference(csv=LEGACY_SPHEROID):
    h, v = _load_csv(csv)
    return {
        "nu_r": _columns(h, v, "nu_r", 6),
        "alpha": v[:, h.index("alpha")],
        "U_r": v[:, h.index("U_r")],
        "tau_lift_drag": _columns(h, v, "tau_lift_drag", 6),
        "tau_crossflow": _columns(h, v, "tau_crossflow", 6),
        "M_RB": _columns(h, v, "M_RB", 36).reshape(-1, 6, 6),  # row-major
        "M_A": _columns(h, v, "M_A", 36).reshape(-1, 6, 6),
        "D": _columns(h, v, "D", 36).reshape(-1, 6, 6),
    }


def _catamaran_reference(csv=LEGACY_CATAMARAN):
    h, v = _load_csv(csv)
    return {
        "nu_r": _columns(h, v, "nu_r", 6),
        "tau_damp": _columns(h, v, "tau_damp", 6),
        "tau_crossflow": _columns(h, v, "tau_crossflow", 6),
        "M": _columns(h, v, "M_total", 36).reshape(-1, 6, 6),
        "G": _columns(h, v, "G", 36).reshape(-1, 6, 6),
        "D": _columns(h, v, "D", 36).reshape(-1, 6, 6),
    }


def _random_nu_r():
    """1000 seeded states: 700 uniform in +-3, 300 with slow cross-flow
    (|v|, |w| < 0.1 m/s, Re < 2e5 on the 1.6 m length: the sub-critical
    branch the 50 reference cases never reach, A-2 section 7)."""
    rng = np.random.default_rng(SEED)
    fast = rng.uniform(-3.0, 3.0, size=(700, 6))
    slow = rng.uniform(-3.0, 3.0, size=(300, 6))
    slow[:, 1:3] = rng.uniform(-0.1, 0.1, size=(300, 2))
    slow[:, 3:6] = rng.uniform(-0.05, 0.05, size=(300, 3))
    return np.vstack([fast, slow])


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _skew(v):
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def _hmtrx(r):
    return np.block([[np.eye(3), _skew(r).T], [np.zeros((3, 3)), np.eye(3)]])


# --------------------------------------------------------------------------
# Parameter sets (every number parsed from a pinned line)
# --------------------------------------------------------------------------
def _spheroid_geometry():
    l_auv, d_auv = _value(REMUS, 132), _value(REMUS, 133)
    s = _value(REMUS, 134, {"L_auv": l_auv, "D_auv": d_auv})
    b = _value(REMUS, 136, {"D_auv": d_auv})
    cd0 = _value(REMUS, 145, {"Cd": _value(REMUS, 144), "b": b, "S": s})
    return {"L": l_auv, "D": d_auv, "S": s, "CD_0": cd0}


def _spheroid_damping(fade, csv=LEGACY_SPHEROID):
    """Matrices and weight as the spheroid reference generator built them."""
    ref = _spheroid_reference(csv)
    assert np.ptp(ref["M_RB"], axis=0).max() == 0.0 and np.ptp(ref["M_A"], axis=0).max() == 0.0
    mrb, ma = ref["M_RB"][0], ref["M_A"][0]
    return {
        "rigid_body_mass_matrix": mrb,
        "added_mass_matrix": ma,
        "weight": mrb[0, 0] * _gravity(_value(REMUS, 97)),  # remus100.m 212
        "center_of_gravity": _value(REMUS, 138),
        "center_of_buoyancy": _value(REMUS, 139),
        "time_constants": [_value(REMUS, 193), _value(REMUS, 194), _value(REMUS, 197)],
        "damping_ratios": [_value(REMUS, 195), _value(REMUS, 196)],
        "sway_damping_fade": fade,
    }


def _scaled_damping(fade):
    """Second vehicle (test construction): masses x8, rotational inertia x32,
    weight x8, CG x2, time constants x sqrt(2)."""
    p = _spheroid_damping(fade)
    scale = np.sqrt(np.array([8.0, 8.0, 8.0, 32.0, 32.0, 32.0]))
    return {**p,
            "rigid_body_mass_matrix": np.outer(scale, scale) * p["rigid_body_mass_matrix"],
            "added_mass_matrix": np.outer(scale, scale) * p["added_mass_matrix"],
            "weight": 8.0 * p["weight"],
            "center_of_gravity": 2.0 * p["center_of_gravity"],
            "time_constants": list(np.sqrt(2.0) * np.asarray(p["time_constants"])),
            }


def _cross_flow_spheroid(grid):
    geo = _spheroid_geometry()
    return {"length": geo["L"], "beam": geo["D"], "draft": geo["D"],   # remus100.m 219
            "water_density": _value(CFD, 24), "drag_model": "cylinder",
            "strip_grid": grid}


def _catamaran_draft():
    m, rho, mp = _value(OTTER, 94), _value(OTTER, 91), _value(GEN_CAT, 12)
    nabla = _value(OTTER, 121, {"m": m, "mp": mp, "rho": rho})
    return _value(OTTER, 122, {"nabla": nabla, "Cb_pont": _value(OTTER, 107),
                               "B_pont": _value(OTTER, 104), "L": _value(OTTER, 92)})


def _cross_flow_catamaran(grid):
    return {"length": _value(OTTER, 92), "beam": _value(OTTER, 104),
            "draft": _catamaran_draft(), "water_density": _value(CFD, 24),
            "drag_model": "hoerner", "strip_grid": grid}                # otter.m 245


def _scaled_cross_flow(params):
    return {**params, "length": 2.0 * params["length"], "beam": 2.0 * params["beam"],
            "draft": 2.0 * params["draft"]}


def _lift_drag_spheroid(rho):
    geo = _spheroid_geometry()
    return {"span": geo["D"], "planform_area": geo["S"],                # remus100.m 218
            "parasitic_drag_coefficient": geo["CD_0"],
            "oswald_efficiency": _value(CLD, 54), "water_density": rho}


def _surface_damping_catamaran(gravity, csv=LEGACY_CATAMARAN):
    """Mass and restoring matrices of the catamaran reference; G moved from
    the CO to the CF with H(-r_f), r_f = [LCF 0 0] (otter.m 192-194)."""
    ref = _catamaran_reference(csv)
    m, g_co = ref["M"][0], ref["G"][0]
    h_inv = _hmtrx(-np.array([_value(OTTER, 192), 0.0, 0.0]))
    return {
        "mass_matrix": m,
        "restoring_matrix": h_inv.T @ g_co @ h_inv,
        "max_forward_thrust": _number_in(OTTER, 202, r"Xu = -([0-9.]+) \*") * gravity,
        "max_speed": _value(OTTER, 101),
        "time_constants": [_value(OTTER, 99), _value(OTTER, 100)],
        "damping_ratios": [_number_in(OTTER, n, r"-2 \* ([0-9.]+) \*") for n in (204, 205, 206)],
        "yaw_damping_nonlinearity": _number_in(OTTER, 240, r"\(1 \+ ([0-9.]+) \* abs"),
    }


def _matlab_gravity_catamaran():
    return _gravity(_value(GEN_CAT, 136))  # generator lines 136-137


# --------------------------------------------------------------------------
# MSS transcriptions (each line cited and pinned)
# --------------------------------------------------------------------------
def _mss_dmtrx(p, nu_r):
    """Dmtrx.m 30, 44-49 (submerged branch) + remus100.m 216 (surge fade only)."""
    m = p["rigid_body_mass_matrix"] + p["added_mass_matrix"]
    t1, t2, t6 = p["time_constants"]
    z4, z5 = p["damping_ratios"]
    t3 = t2
    dz = p["center_of_gravity"][2] - p["center_of_buoyancy"][2]
    w4 = np.sqrt(p["weight"] * dz / m[3, 3])
    w5 = np.sqrt(p["weight"] * dz / m[4, 4])
    d = np.diag([m[0, 0] / t1, m[1, 1] / t2, m[2, 2] / t3,
                 m[3, 3] * 2 * z4 * w4, m[4, 4] * 2 * z5 * w5, m[5, 5] / t6])
    d[0, 0] *= np.exp(-3 * np.linalg.norm(nu_r[:3]))  # remus100.m 128, 216
    return d


def _mss_otter_damping(p, nu_r):
    """otter.m 197-207 and 235-240; returns (derivatives, tau_damp)."""
    m, g_cf = p["mass_matrix"], p["restoring_matrix"]
    w3, w4, w5 = (np.sqrt(g_cf[i, i] / m[i, i]) for i in (2, 3, 4))
    z3, z4, z5 = p["damping_ratios"]
    t_sway, t_yaw = p["time_constants"]
    deriv = np.array([-p["max_forward_thrust"] / p["max_speed"], -m[1, 1] / t_sway,
                      -2 * z3 * w3 * m[2, 2], -2 * z4 * w4 * m[3, 3],
                      -2 * z5 * w5 * m[4, 4], -m[5, 5] / t_yaw])
    tau = deriv * nu_r
    tau[5] *= 1 + p["yaw_damping_nonlinearity"] * abs(nu_r[5])
    return deriv, tau


def _interp1_clamped(x, table):
    """cylinderDrag.m 80-86: interp1 inside the table, end values outside."""
    if x < table[0, 0]:
        return table[0, 1]
    if x > table[-1, 0]:
        return table[-1, 1]
    return float(np.interp(x, table[:, 0], table[:, 1]))


def _mss_cylinder_cd(length, beam, nu_r, reynolds_length=None):
    """cylinderDrag.m 25-110 at ac77394: ``Re = U_cf B / nu_water`` (lines
    78-80). ``reynolds_length=length`` gives the pre-fix line of ``99bf0b3``
    (``Re = U_crossflow * L * 1e6``, line 77 there) for the history tests only."""
    cd_data = _table(CYL, "CD_DATA")
    u_cf = np.sqrt(nu_r[1] ** 2 + nu_r[2] ** 2)            # line 78
    nu_water = _value(CYL, 79)                             # line 79
    re_number = u_cf * (beam if reynolds_length is None else reynolds_length) / nu_water  # line 80
    cd = _interp1_clamped(re_number, cd_data)              # lines 83-89
    kappa_table = _table(CYL, "KAPPA_SUBCRITICAL_DATA") if re_number < 2e5 \
        else _table(CYL, "KAPPA_SUPERCRITICAL_DATA")       # lines 92-108
    return cd * _interp1_clamped(length / beam, kappa_table)


def _mss_hoerner(beam, draft):
    """Hoerner.m 25-51."""
    data = _table(HOERNER, "CD_DATA")
    ratio = beam / (2 * draft)
    return float(np.interp(ratio, data[:, 0], data[:, 1])) if ratio <= data[-1, 0] \
        else _number_in(HOERNER, 50, r"= ([0-9.]+);")


def _mss_cross_flow(p, nu_r, history_reynolds_on_length=False):
    """crossFlowDrag.m 24-54 (current: 20 midpoints); ``strip_grid="endpoint"``
    gives the pre-2026-08-26 loop ``for xL = -L/2:dx:L/2`` (A-2 section 1).
    ``history_reynolds_on_length``: the cylinder Re of MSS before ``ac77394``."""
    length, beam, draft = p["length"], p["beam"], p["draft"]
    if p["drag_model"] == "cylinder":
        cd = _mss_cylinder_cd(length, beam, nu_r, length if history_reynolds_on_length else None)
    else:
        cd = _mss_hoerner(beam, draft)
    n_strips = int(_value(CFD, 25))
    dx = length / n_strips
    if p.get("strip_grid", "midpoint") == "midpoint":
        xs = -length / 2 + (np.arange(1, n_strips + 1) - 0.5) * dx   # line 41
    else:
        xs = -length / 2 + dx * np.arange(n_strips + 1)
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    u_h = np.abs(v_r + xs * r) * (v_r + xs * r)                      # line 46
    u_v = np.abs(w_r + xs * q) * (w_r + xs * q)                      # line 47
    k = -0.5 * p["water_density"] * draft * cd * dx                 # lines 48-51
    return np.array([0.0, k * u_h.sum(), k * u_v.sum(), 0.0, k * (xs * u_v).sum(), k * (xs * u_h).sum()])


def _mss_lift_drag(p, nu_r):
    """forceLiftDrag.m 26-40 + coeffLiftDrag.m 54-68 (sigma = 0), remus100.m 127-128."""
    alpha = np.arctan2(nu_r[2], nu_r[0])
    u_r = np.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2)
    b, s, e = p["span"], p["planform_area"], p["oswald_efficiency"]
    ar = b ** 2 / s
    cl = np.pi * ar / (1 + np.sqrt(1 + (ar / 2) ** 2)) * alpha
    cd = p["parasitic_drag_coefficient"] + cl ** 2 / (np.pi * e * ar)
    f_drag = 0.5 * p["water_density"] * u_r ** 2 * s * cd
    f_lift = 0.5 * p["water_density"] * u_r ** 2 * s * cl
    return np.array([np.cos(alpha) * (-f_drag) - np.sin(alpha) * (-f_lift), 0.0,
                     np.sin(alpha) * (-f_drag) + np.cos(alpha) * (-f_lift), 0.0, 0.0, 0.0])


# --------------------------------------------------------------------------
# Builders and numpy source
# --------------------------------------------------------------------------
def _build_submerged(p):
    block = _contract(LINEAR_DAMPING)
    c = block.preprocess_submerged_linear_damping(**p)
    return c, block.submerged_linear_damping_casadi(c)


def _build_surface(p):
    block = _contract(LINEAR_DAMPING)
    c = block.preprocess_surface_linear_damping(**p)
    return c, block.surface_linear_damping_casadi(c)


def _build_cross_flow(p):
    block = _contract(CROSS_FLOW)
    return block.cross_flow_drag_casadi(block.preprocess_cross_flow_drag(**p))


def _build_lift_drag(p):
    block = _contract(LIFT_DRAG)
    return block.lift_drag_casadi(block.preprocess_lift_drag(**p))


def _attribute(module_name, name):
    block = _contract(module_name)
    if not hasattr(block, name):
        pytest.fail(f"block not ported yet: {module_name}.{name}")
    return getattr(block, name)


def _build_floating(p):
    c = _attribute(LINEAR_DAMPING, "preprocess_floating_linear_damping")(**p)
    return c, _attribute(LINEAR_DAMPING, "floating_linear_damping_casadi")(c)


def _build_surge(p):
    p = dict(p)
    model = p.pop("quadratic_model")
    c = _attribute(SURGE_DAMPING, f"preprocess_surge_damping_{model}")(**p)
    return c, _attribute(SURGE_DAMPING, "surge_damping_casadi")(c)


# --------------------------------------------------------------------------
# U3a references: forceSurgeDamping.m, XuuITTC.m, Dmtrx.m surface branch
# (MATLAB direct calls at MSS ac77394, generate_u3_surge_floating_mss.m)
# --------------------------------------------------------------------------
def _named_columns(csv):
    header, values = _load_csv(csv)
    return {name: values[:, i] for i, name in enumerate(header)}


def _surge_parameters(set_id, model):
    """Block parameters of one parameter set of the surge CSV (set 1 = osv.m)."""
    ref = _named_columns(MSS_SURGE)
    row = int(np.flatnonzero(ref["set"] == set_id)[0])
    common = {"mass": ref["m"][row], "length": ref["L"][row], "water_density": ref["rho"][row],
              "time_constant": ref["T1"][row], "crossover_speed": _value(FSD, 57)}
    if model == "ittc":
        return {"quadratic_model": "ittc", **common, "wetted_surface": ref["S"][row],
                "form_factor": _value(FSD, 70)}
    return {"quadratic_model": "max_thrust", **common, "max_speed": ref["u_max"][row],
            "max_thrust": ref["thrust_max"][row]}


def _surge_rows(set_id, model, probe=0):
    ref = _named_columns(MSS_SURGE)
    branch = 0 if model == "ittc" else 1
    rows = np.flatnonzero((ref["set"] == set_id) & (ref["branch"] == branch) & (ref["probe"] == probe))
    return {k: v[rows] for k, v in ref.items()}


def _mss_added_mass_surge(m, length, rho):
    """addedMassSurge.m 33-34 (Soding 1982): A11 = 2.7 rho nabla^(5/3) / L^2,
    the exponent and the L^2 pinned by the line text."""
    factor = _number_in(AMS, 34, r"A11 = ([0-9.]+) \*")
    return factor * rho * (m / rho) ** (5 / 3) / length ** 2


def _mss_force_surge_damping(u, m, s, length, t1, rho, u_max, thrust_max=None):
    """forceSurgeDamping.m 57-82; returns (X, Xuu, Xu, A11)."""
    a11 = _mss_added_mass_surge(m, length, rho)
    xu = -(m - (-a11)) / t1                                            # lines 60-61
    if thrust_max is not None:
        xuu = -thrust_max / u_max ** 2                                 # line 66
    else:
        rn = (length / _value(FSD, 69)) * abs(u)                       # line 72
        cf = 0.075 / (np.log10(rn + _value(FSD, 71)) - 2) ** 2         # line 73
        xuu = -0.5 * rho * s * (1 + _value(FSD, 70)) * cf              # line 74
    sigma = 1 - np.tanh(u / _value(FSD, 57))                           # line 79
    return sigma * xu * u + (1 - sigma) * xuu * abs(u) * u, xuu, xu, a11   # line 82


def _mss_xuu_ittc(u, rho, length, beam, draft, c_b):
    """XuuITTC.m 32-39."""
    re_number = max(length * abs(u) / _value(XITTC, 32), _value(XITTC, 34))
    cf = 0.075 / (np.log10(re_number) - 2) ** 2
    wetted = 1.025 * length * (c_b * beam + 1.7 * draft)            # Mumford, line 38
    return -0.5 * rho * wetted * (1 + _value(XITTC, 33)) * cf


def _floating_reference():
    h, v = _load_csv(MSS_FLOATING)
    ref = {p: _columns(h, v, p, 36).reshape(-1, 6, 6) for p in ("MRB", "MA", "G", "D")}
    ref["T"] = v[:, [h.index(n) for n in ("T1", "T2", "T6")]]
    ref["zeta"] = v[:, [h.index(n) for n in ("zeta4", "zeta5")]]
    return ref


def _floating_parameters(k, ref=None):
    ref = _floating_reference() if ref is None else ref
    return {"rigid_body_mass_matrix": ref["MRB"][k], "added_mass_matrix": ref["MA"][k],
            "restoring_matrix": ref["G"][k], "time_constants": list(ref["T"][k]),
            "damping_ratios": [_value(DMTRX, 57), *ref["zeta"][k]]}   # zeta3 = 0.2, Dmtrx.m 57


def _mss_dmtrx_floating(p):
    """Dmtrx.m 30, 51-63 (surface-craft branch)."""
    m = p["rigid_body_mass_matrix"] + p["added_mass_matrix"]
    g = p["restoring_matrix"]
    t1, t2, t6 = p["time_constants"]
    z3, z4, z5 = p["damping_ratios"]
    w3, w4, w5 = (np.sqrt(g[i, i] / m[i, i]) for i in (2, 3, 4))
    return np.diag([m[0, 0] / t1, m[1, 1] / t2, m[2, 2] * 2 * z3 * w3,
                    m[3, 3] * 2 * z4 * w4, m[4, 4] * 2 * z5 * w5, m[5, 5] / t6])


def _source():
    """The numpy source, from MORE_GENERIC_MODELS_DIR only (agents-more rule 9)."""
    root = str(_env_dir("MORE_GENERIC_MODELS_DIR"))
    if root not in sys.path:
        sys.path.insert(0, root)
    try:
        from more_generic_models.dynamics.plant.asv_catamaran.asv_catamaran import ASVCatamaran
        from more_generic_models.dynamics.plant.auv_spheroid.auv_spheroid import AUVSpheroid
        from more_generic_models.dynamics.plant.matrices.hydrodynamics import HydroForces
        from more_generic_models.dynamics.plant.matrices.linear_damping import DampingMatrix
    except ImportError as exc:
        pytest.skip(f"numpy source not importable from MORE_GENERIC_MODELS_DIR={root}: {exc}")
    return {"spheroid": AUVSpheroid, "catamaran": ASVCatamaran,
            "hydro": HydroForces, "damping": DampingMatrix}


@pytest.fixture
def full_table_source(monkeypatch):
    """The numpy source with its rounded Re column replaced, in memory only, by
    the full-precision ``cylinderDrag.m`` table (A-2 section 6, edit 1)."""
    src = _source()
    monkeypatch.setattr(src["hydro"], "_CYLINDER_CD_DATA", _table(CYL, "CD_DATA"))
    return src


# --------------------------------------------------------------------------
# Reference sanity and pinned lines (no port needed)
# --------------------------------------------------------------------------
def test_reference_shapes():
    sph, cat = _spheroid_reference(), _catamaran_reference()
    for key in ("nu_r", "tau_lift_drag", "tau_crossflow"):
        assert sph[key].shape == (50, 6), key
    for key in ("M_RB", "M_A", "D"):
        assert sph[key].shape == (50, 6, 6), key
    for key in ("nu_r", "tau_damp", "tau_crossflow"):
        assert cat[key].shape == (50, 6), key
    for key in ("M", "G", "D"):
        assert cat[key].shape == (50, 6, 6), key
    # U_r and alpha of the reference are remus100.m 127-128 of its nu_r
    np.testing.assert_allclose(sph["U_r"], np.linalg.norm(sph["nu_r"][:, :3], axis=1), atol=1e-12)
    np.testing.assert_allclose(sph["alpha"], np.arctan2(sph["nu_r"][:, 2], sph["nu_r"][:, 0]), atol=1e-12)


def test_snapshot_holds_every_cited_line_and_table():
    """Runs everywhere: the frozen snapshot carries each pinned line (starting
    with its CITED_LINES text) and the four tables in their MSS shapes."""
    for (rel, number), text in CITED_LINES.items():
        assert _line(rel, number).startswith(text), (rel, number)
    assert _table(CYL, "CD_DATA").shape == (28, 2)
    assert _table(CYL, "KAPPA_SUBCRITICAL_DATA").shape == (7, 2)
    assert _table(CYL, "KAPPA_SUPERCRITICAL_DATA").shape == (7, 2)
    assert _table(HOERNER, "CD_DATA").shape == (20, 2)


@pytest.mark.parametrize("variable", ["MSS_DIR", "MORE_GENERIC_MODELS_DIR"])
def test_cited_lines_are_unchanged(variable):
    """With the variable set: every cited line and table of that checkout equals
    the snapshot (a moved or edited line fails here, as cylinderDrag.m did at
    ac77394). Without it: skips, naming the variable (agents-more rule 9)."""
    _env_dir(variable)
    snapshot = _snapshot()
    for (rel, number), text in CITED_LINES.items():
        if rel.startswith(f"${variable}/"):
            live = _live_line(rel, number)
            assert live.startswith(text), (rel, number, live)
            assert live == snapshot["lines"][rel][str(number)], (rel, number)
    if variable == "MSS_DIR":
        for rel, name in TABLES:
            live = _parse_table(_live_path(rel).read_text(), name)
            assert np.array_equal(live, _table(rel, name)), (rel, name)


def test_transcriptions_reproduce_the_stored_matlab_reference():
    """The test-side MSS transcriptions, run on the legacy settings, land on the
    stored MATLAB numbers; so they are fit to serve as default-path references."""
    sph, cat = _spheroid_reference(), _catamaran_reference()
    p_d = _spheroid_damping(True)
    p_cf = _cross_flow_spheroid("endpoint")
    p_ld = _lift_drag_spheroid(_value(GEN_SPH, 467))
    p_cat = _surface_damping_catamaran(_matlab_gravity_catamaran())
    p_cat_cf = _cross_flow_catamaran("endpoint")
    for k in range(50):
        d = _mss_dmtrx(p_d, sph["nu_r"][k])
        d[1, 1] *= np.exp(-3 * sph["U_r"][k])   # generator line 268 (sway fade)
        assert _max_diff(d, sph["D"][k]) <= G1_TOLERANCE, k
        tau_cf = _mss_cross_flow(p_cf, sph["nu_r"][k], history_reynolds_on_length=True)  # pre-ac77394 Re
        assert _max_diff(tau_cf, sph["tau_crossflow"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_lift_drag(p_ld, sph["nu_r"][k]), sph["tau_lift_drag"][k]) <= G1_TOLERANCE, k
        deriv, tau = _mss_otter_damping(p_cat, cat["nu_r"][k])
        assert _max_diff(np.diag(deriv), cat["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(tau, cat["tau_damp"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cat_cf, cat["nu_r"][k]), cat["tau_crossflow"][k]) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G1 — numpy source vs the stored MATLAB reference
# --------------------------------------------------------------------------
def test_G1_numpy_source_spheroid_damping_and_lift_drag():
    vehicle = _source()["spheroid"]()
    ref = _spheroid_reference()
    for k in range(50):
        assert _max_diff(vehicle.get_D(ref["U_r"][k]), ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        tau = vehicle.get_tau_lift_drag(ref["alpha"][k], ref["U_r"][k])
        assert _max_diff(tau, ref["tau_lift_drag"][k]) <= G1_TOLERANCE, ("lift_drag", k)


def test_G1_numpy_source_spheroid_cross_flow_with_full_table(full_table_source):
    vehicle = full_table_source["spheroid"]()
    ref = _spheroid_reference()
    for k in range(50):
        tau = vehicle.get_tau_cross_flow(ref["nu_r"][k], "cylinder")
        assert _max_diff(tau, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_G1_numpy_source_rounded_table_residual_is_pinned():
    """A-2: the unpatched source misses MATLAB by the rounded Re column only."""
    vehicle = _source()["spheroid"]()
    ref = _spheroid_reference()
    worst = max(_max_diff(vehicle.get_tau_cross_flow(ref["nu_r"][k], "cylinder"),
                          ref["tau_crossflow"][k]) for k in range(50))
    assert G1_TOLERANCE < worst <= ROUNDED_TABLE_RESIDUAL_BOUND, worst


def test_numpy_source_drag_models_copy_equals_hydro_forces_copy():
    """``drag_models.py::DragModels`` (map rows -> cross_flow.py) carries a second
    copy of both coefficient tables; it agrees with the ``HydroForces`` copy the
    vehicles call, so the G1/G2 tests on ``HydroForces`` cover both rows."""
    src = _source()
    from more_generic_models.dynamics.plant.matrices.drag_models import DragModels
    rng = np.random.default_rng(SEED)
    for b, t in rng.uniform(0.05, 2.0, size=(200, 2)):
        assert abs(DragModels.hoerner(b, t) - src["hydro"]._hoerner(b, t)) <= G2_TOLERANCE, (b, t)
    for nu_r in _random_nu_r():
        length, beam = rng.uniform(0.5, 6.0), rng.uniform(0.05, 0.5)
        assert abs(DragModels.cylinder_drag(length, beam, nu_r)
                   - src["hydro"]._cylinder_drag(length, beam, nu_r)) <= G2_TOLERANCE


def test_G1_numpy_source_catamaran_damping_and_cross_flow():
    vessel = _source()["catamaran"]()
    ref = _catamaran_reference()
    for k in range(50):
        assert _max_diff(vessel.get_D(), ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(vessel.get_tau_damping(ref["nu_r"][k]), ref["tau_damp"][k]) <= G1_TOLERANCE, k
        assert _max_diff(vessel.get_tau_cross_flow(ref["nu_r"][k]), ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G1 — block (legacy settings) vs the stored MATLAB reference
# --------------------------------------------------------------------------
def test_G1_block_spheroid_damping_with_sway_fade_flag():
    constants, function = _build_submerged(_spheroid_damping(True))
    ref = _spheroid_reference()
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], ref["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(out["tau"].reshape(-1), -ref["D"][k] @ ref["nu_r"][k]) <= G1_TOLERANCE, k


# U3a: ``test_G1_block_spheroid_cross_flow_legacy_settings`` is removed with the
# ``cross_flow_reynolds_length`` flag (owner 2026-10-06): the legacy spheroid file
# mixes the old Re-on-length line with the endpoint grid, so no flag-free block
# reproduces it. The endpoint grid keeps its G1 on the catamaran (Hoerner, no
# Reynolds number) below and its G2 on both vehicles.


def test_G1_block_spheroid_lift_drag_template_density():
    function = _build_lift_drag(_lift_drag_spheroid(_value(GEN_SPH, 467)))
    ref = _spheroid_reference()
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_lift_drag"][k]) <= G1_TOLERANCE, k


def test_G1_block_catamaran_surface_damping_and_cross_flow():
    _, damping = _build_surface(_surface_damping_catamaran(_matlab_gravity_catamaran()))
    cross_flow = _build_cross_flow(_cross_flow_catamaran("endpoint"))
    ref = _catamaran_reference()
    for k in range(50):
        out = _evaluate(damping, ref["nu_r"][k])
        assert _max_diff(out["D"], -ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), ref["tau_damp"][k]) <= G1_TOLERANCE, ("tau", k)
        tau_cf = _evaluate(cross_flow, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(tau_cf, ref["tau_crossflow"][k]) <= G1_TOLERANCE, ("cross_flow", k)


# --------------------------------------------------------------------------
# G2 — CasADi blocks vs numpy source (reference cases + 1000 seeded states)
# --------------------------------------------------------------------------
def _g2_states(reference):
    return np.vstack([reference["nu_r"], _random_nu_r()])


def test_G2_function_signatures():
    for function, outs in (
        (_build_submerged(_spheroid_damping(False))[1], ["D", "tau"]),
        (_build_surface(_surface_damping_catamaran(_matlab_gravity_catamaran()))[1], ["D", "tau"]),
        (_build_cross_flow(_cross_flow_spheroid("midpoint")), ["tau"]),
        (_build_lift_drag(_lift_drag_spheroid(_value(FLD, 26))), ["tau"]),
        (_build_floating(_floating_parameters(0))[1], ["D", "tau"]),
        (_build_surge(_surge_parameters(1, "ittc"))[1], ["tau"]),
        (_build_surge(_surge_parameters(1, "max_thrust"))[1], ["tau"]),
    ):
        assert function.name_in() == ["nu_r"]
        assert function.size_in(0) == (6, 1)
        assert function.name_out() == outs
        sizes = {"D": (6, 6), "tau": (6, 1)}
        for i, name in enumerate(outs):
            assert function.size_out(i) == sizes[name], name


@pytest.mark.parametrize("scaled", [False, True], ids=["remus_like", "scaled"])
@pytest.mark.parametrize("fade", [False, True], ids=["default", "sway_fade"])
def test_G2_submerged_damping_matches_source(fade, scaled):
    damping = _source()["damping"]
    p = _scaled_damping(fade) if scaled else _spheroid_damping(fade)
    constants, function = _build_submerged(p)
    coeff = damping.get_damping_coeff(p["weight"], p["center_of_gravity"], p["center_of_buoyancy"],
                                      p["rigid_body_mass_matrix"], p["added_mass_matrix"],
                                      p["time_constants"], p["damping_ratios"])
    assert _max_diff(constants.damping_coefficients, coeff) <= G2_TOLERANCE
    for k, nu_r in enumerate(_g2_states(_spheroid_reference())):
        expected = damping.D_linear(coeff, np.linalg.norm(nu_r[:3]))  # fades sway too
        if not fade:
            expected[1, 1] = coeff[1]   # default: surge-only fade, remus100.m 216
        out = _evaluate(function, nu_r)
        assert _max_diff(out["D"], expected) <= G2_TOLERANCE, k
        assert _max_diff(out["tau"].reshape(-1), -expected @ nu_r) <= G2_TOLERANCE, k


@pytest.mark.parametrize("scale", [1.0, 3.0], ids=["otter_like", "scaled"])
def test_G2_surface_damping_matches_source(scale):
    damping = _source()["damping"]
    g = _matlab_gravity_catamaran()
    p = _surface_damping_catamaran(g)
    p = {**p, "mass_matrix": scale ** 3 * p["mass_matrix"],
         "restoring_matrix": scale ** 2 * p["restoring_matrix"],
         "max_forward_thrust": scale ** 3 * p["max_forward_thrust"]}
    constants, function = _build_surface(p)
    # source D_linear_catamaran(G_CF, MRB, MA, Umax, g, zeta_345, T_126, coeff):
    # Xu = coeff[0] g / Umax, the rest -1 x the otter.m forms (asv_catamaran_params.py 25-30)
    coeff = np.array([-p["max_forward_thrust"] / g, -1.0, -1.0, -1.0, -1.0, -1.0])
    t_sway, t_yaw = p["time_constants"]
    d_src = damping.D_linear_catamaran(p["restoring_matrix"], p["mass_matrix"], np.zeros((6, 6)),
                                       p["max_speed"], g, p["damping_ratios"],
                                       [np.nan, t_sway, t_yaw], coeff)
    assert _max_diff(constants.damping_derivatives, np.diag(d_src)) <= G2_TOLERANCE
    for k, nu_r in enumerate(_g2_states(_catamaran_reference())):
        out = _evaluate(function, nu_r)
        assert _max_diff(out["D"], -d_src) <= G2_TOLERANCE, k
        # get_tau_damping hard-codes the yaw factor 10 (linear_damping.py 221)
        assert _max_diff(out["tau"].reshape(-1), damping.get_tau_damping(d_src, nu_r)) <= G2_TOLERANCE, k


@pytest.fixture
def diameter_reynolds_source(full_table_source, monkeypatch):
    """The numpy source with the full table and its cylinder Reynolds number on
    the diameter (MSS ``ac77394``), in memory only and by argument mapping
    alone: ``_cylinder_drag(L, B, nu_r)`` puts Re on its first argument and
    kappa on the ratio of the two, so calling the source's own function with
    ``(B, B^2/L)`` gives ``Re = U B 1e6`` and ``kappa(L/B)`` (no equation of the
    source is rewritten)."""
    hydro = full_table_source["hydro"]
    original = hydro._cylinder_drag
    monkeypatch.setattr(hydro, "_cylinder_drag",
                        staticmethod(lambda L, B, nu_r: original(B, B * B / L, nu_r)))
    return full_table_source


@pytest.mark.parametrize("vehicle", ["spheroid", "spheroid_scaled", "catamaran", "catamaran_scaled"])
def test_G2_cross_flow_legacy_settings_match_source(vehicle, diameter_reynolds_source):
    """The endpoint-grid flag vs the numpy source (21 end points); the spheroid
    runs the source with Re on the diameter (fixture), the catamaran needs no
    Reynolds number (Hoerner)."""
    hydro = diameter_reynolds_source["hydro"]
    p = (_cross_flow_spheroid("endpoint") if vehicle.startswith("spheroid")
         else _cross_flow_catamaran("endpoint"))
    if vehicle.endswith("scaled"):
        p = _scaled_cross_flow(p)
    function = _build_cross_flow(p)
    ref = _spheroid_reference() if vehicle.startswith("spheroid") else _catamaran_reference()
    for k, nu_r in enumerate(_g2_states(ref)):
        expected = hydro.cross_flow_drag(p["length"], p["beam"], p["draft"], nu_r,
                                         p["drag_model"], p["water_density"])
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, expected) <= G2_TOLERANCE, (vehicle, k)


@pytest.mark.parametrize("scale", [1.0, 2.0], ids=["remus_like", "scaled"])
def test_G2_lift_drag_matches_source(scale):
    hydro = _source()["hydro"]
    p = _lift_drag_spheroid(_value(GEN_SPH, 467))
    p = {**p, "span": scale * p["span"], "planform_area": scale ** 2 * p["planform_area"]}
    function = _build_lift_drag(p)
    for k, nu_r in enumerate(_g2_states(_spheroid_reference())):
        alpha, u_r = np.arctan2(nu_r[2], nu_r[0]), np.linalg.norm(nu_r[:3])
        expected = hydro.force_lift_drag(p["span"], p["planform_area"], p["parasitic_drag_coefficient"],
                                         alpha, u_r, p["water_density"], p["oswald_efficiency"])
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, expected) <= G2_TOLERANCE, k


# --------------------------------------------------------------------------
# G4 — a broken model must fail (control first, then one perturbation each)
# --------------------------------------------------------------------------
def _worst(function, cases, expected, key):
    return max(_max_diff(_evaluate(function, nu_r)[key].reshape(np.shape(e)), e)
               for nu_r, e in zip(cases, expected))


def test_G4_submerged_damping_roll_ratio_times_zero_is_detected():
    ref = _spheroid_reference()
    p = _spheroid_damping(True)
    assert _worst(_build_submerged(p)[1], ref["nu_r"], ref["D"], "D") <= G1_TOLERANCE
    broken = {**p, "damping_ratios": [0.0, p["damping_ratios"][1]]}
    assert _worst(_build_submerged(broken)[1], ref["nu_r"], ref["D"], "D") > G4_FACTOR * G1_TOLERANCE


def test_G4_surface_damping_heave_ratio_times_zero_is_detected():
    ref = _catamaran_reference()
    p = _surface_damping_catamaran(_matlab_gravity_catamaran())
    assert _worst(_build_surface(p)[1], ref["nu_r"], ref["tau_damp"], "tau") <= G1_TOLERANCE
    broken = {**p, "damping_ratios": [0.0] + list(p["damping_ratios"][1:])}
    assert _worst(_build_surface(broken)[1], ref["nu_r"], ref["tau_damp"], "tau") > G4_FACTOR * G1_TOLERANCE


def test_G4_cross_flow_draft_plus_1_percent_is_detected():
    """Default path against MATLAB at ac77394 (U3a; was the legacy file + flag)."""
    ref = _spheroid_reference(MSS_AC_SPHEROID)
    p = {k: v for k, v in _cross_flow_spheroid("midpoint").items() if k != "strip_grid"}
    assert _worst(_build_cross_flow(p), ref["nu_r"], ref["tau_crossflow"], "tau") <= G1_TOLERANCE
    broken = {**p, "draft": 1.01 * p["draft"]}
    assert _worst(_build_cross_flow(broken), ref["nu_r"], ref["tau_crossflow"], "tau") > G4_FACTOR * G1_TOLERANCE


def test_G4_lift_drag_area_plus_1_percent_is_detected():
    ref = _spheroid_reference()
    p = _lift_drag_spheroid(_value(GEN_SPH, 467))
    assert _worst(_build_lift_drag(p), ref["nu_r"], ref["tau_lift_drag"], "tau") <= G1_TOLERANCE
    broken = {**p, "planform_area": 1.01 * p["planform_area"]}
    assert _worst(_build_lift_drag(broken), ref["nu_r"], ref["tau_lift_drag"], "tau") > G4_FACTOR * G1_TOLERANCE


# --------------------------------------------------------------------------
# Headline (E-18, E-20 Q2 a, E-10 a, E-20 Q3 b): default path vs current MSS
# --------------------------------------------------------------------------
def _states():
    return np.vstack([_spheroid_reference()["nu_r"], _random_nu_r()])


def test_MSS_submerged_damping_default_fades_surge_only():
    """Default = Dmtrx.m + remus100.m:216; the template's sway fade is the flag."""
    p = _spheroid_damping(False)
    _, function = _build_submerged(p)
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)
        assert _max_diff(out["D"], _mss_dmtrx(p, nu_r)) <= G1_TOLERANCE, ("Dmtrx.m 30-49, remus100.m 216", k)


def test_MSS_surface_damping_equals_otter():
    """otter.m 197-207, 235-240 with otter.m's own g = 9.81 (line 90)."""
    p = _surface_damping_catamaran(_value(OTTER, 90))
    constants, function = _build_surface(p)
    for k, nu_r in enumerate(np.vstack([_catamaran_reference()["nu_r"], _random_nu_r()])):
        deriv, tau = _mss_otter_damping(p, nu_r)
        out = _evaluate(function, nu_r)
        assert _max_diff(out["tau"].reshape(-1), tau) <= G1_TOLERANCE, ("otter.m 235-240", k)
    assert _max_diff(constants.damping_derivatives, deriv) <= G1_TOLERANCE


@pytest.mark.parametrize("model", ["cylinder", "hoerner"])
def test_MSS_cross_flow_default_equals_current_mss(model):
    """Default (no flags) = crossFlowDrag.m (6a2a064) with cylinderDrag.m at
    ac77394 (Re on the diameter, lines 78-80) on 1050 states (U3a; replaces
    the length-flag test)."""
    p = _cross_flow_spheroid("midpoint") if model == "cylinder" else _cross_flow_catamaran("midpoint")
    function = _build_cross_flow({k: v for k, v in p.items() if k != "strip_grid"})
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, _mss_cross_flow(p, nu_r)) <= G1_TOLERANCE, ("crossFlowDrag.m 24-54", k)


def test_D_MSS_1_flag_is_gone():
    """Owner 2026-10-06 (E-20, "a"): with MSS on the diameter since ac77394 the
    ``cross_flow_reynolds_length`` flag is dropped; ``strip_grid`` stays with
    default "midpoint". Fails until the porter (U3b) removes the flag."""
    signature = inspect.signature(_contract(CROSS_FLOW).preprocess_cross_flow_drag)
    assert "cross_flow_reynolds_length" not in signature.parameters, list(signature.parameters)
    assert signature.parameters["strip_grid"].default == "midpoint"


def test_E10_endpoint_grid_flag_reproduces_the_template():
    """strip_grid="endpoint" equals the pre-2026-08-26 grid of the stored
    catamaran reference (Hoerner: no Reynolds number); the default grid does
    not (U3a: was the spheroid file, which needs the dropped length flag)."""
    ref = _catamaran_reference()
    legacy = _build_cross_flow(_cross_flow_catamaran("endpoint"))
    current = _build_cross_flow(_cross_flow_catamaran("midpoint"))
    worst_legacy = max(_max_diff(_evaluate(legacy, n)["tau"].reshape(-1), ref["tau_crossflow"][k])
                       for k, n in enumerate(ref["nu_r"]))
    worst_current = max(_max_diff(_evaluate(current, n)["tau"].reshape(-1), ref["tau_crossflow"][k])
                        for k, n in enumerate(ref["nu_r"]))
    assert worst_legacy <= G1_TOLERANCE
    assert worst_current > G4_FACTOR * G1_TOLERANCE


def test_MSS_lift_drag_equals_force_lift_drag():
    """forceLiftDrag.m rho = 1026 (line 26), e = 0.3 (coeffLiftDrag.m 54)."""
    p = _lift_drag_spheroid(_value(FLD, 26))
    function = _build_lift_drag(p)
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, _mss_lift_drag(p, nu_r)) <= G1_TOLERANCE, ("forceLiftDrag.m 26-40", k)


# --------------------------------------------------------------------------
# G1-MSS (job A-30): default paths vs MATLAB running current MSS (job A-26)
# --------------------------------------------------------------------------
def test_transcriptions_reproduce_the_current_mss_matlab_reference():
    """The second check: the test-side transcriptions, on current-MSS settings,
    land on the A-26 files (A-26 ledger step 4 c, now a test) and, for the
    spheroid cross-flow, on the ac77394 file (U3a). The A-26 spheroid
    cross-flow is MSS 99bf0b3 (Re on the length): history."""
    sph, cat = _spheroid_reference(MSS_SPHEROID), _catamaran_reference(MSS_CATAMARAN)
    sph_ac = _spheroid_reference(MSS_AC_SPHEROID)
    p_d = _spheroid_damping(False, MSS_SPHEROID)
    p_cf = _cross_flow_spheroid("midpoint")
    p_ld = _lift_drag_spheroid(_value(FLD, 26))
    p_cat = _surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN)
    p_cat_cf = _cross_flow_catamaran("midpoint")
    for k in range(50):
        assert _max_diff(_mss_dmtrx(p_d, sph["nu_r"][k]), sph["D"][k]) <= G1_TOLERANCE, k
        history = _mss_cross_flow(p_cf, sph["nu_r"][k], history_reynolds_on_length=True)
        assert _max_diff(history, sph["tau_crossflow"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cf, sph_ac["nu_r"][k]), sph_ac["tau_crossflow"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_lift_drag(p_ld, sph["nu_r"][k]), sph["tau_lift_drag"][k]) <= G1_TOLERANCE, k
        deriv, tau = _mss_otter_damping(p_cat, cat["nu_r"][k])
        assert _max_diff(np.diag(deriv), cat["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(tau, cat["tau_damp"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cat_cf, cat["nu_r"][k]), cat["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_G1_MSS_submerged_damping_default_equals_matlab():
    """Default (surge fade only, remus100.m 216) with M_RB, M_A of the A-26 file
    (imlay61.m rho = 1026) vs MATLAB's D and -D nu_r."""
    ref = _spheroid_reference(MSS_SPHEROID)
    p = _spheroid_damping(False, MSS_SPHEROID)
    p = {k: v for k, v in p.items() if k != "sway_damping_fade"}   # the default
    _, function = _build_submerged(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), -ref["D"][k] @ ref["nu_r"][k]) <= G1_TOLERANCE, ("tau", k)


# U3a: ``test_G1_MSS_cross_flow_length_flag_equals_matlab`` is removed with the flag
# (owner 2026-10-06); its MATLAB file (MSS 99bf0b3) stays, read by the history
# checks above and below.


def test_G1_MSS_cross_flow_default_equals_matlab_ac77394():
    """D-MSS-1 resolved: the default (Re on the diameter) equals MATLAB running
    MSS ac77394 (the A-26 generator re-run unchanged, SOURCE.md)."""
    ref = _spheroid_reference(MSS_AC_SPHEROID)
    p = {k: v for k, v in _cross_flow_spheroid("midpoint").items() if k != "strip_grid"}
    function = _build_cross_flow(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_ac77394_reference_differs_from_99bf0b3_only_in_cross_flow():
    """The re-run at ac77394 changes exactly the four non-zero tau_crossflow
    columns (the issue #81 fix); every other column is byte-equal as text."""
    header = (DATA_DIR / MSS_AC_SPHEROID).read_text().splitlines()[0].split(",")
    assert header == (DATA_DIR / MSS_SPHEROID).read_text().splitlines()[0].split(",")
    new = [row.split(",") for row in (DATA_DIR / MSS_AC_SPHEROID).read_text().splitlines()[1:]]
    old = [row.split(",") for row in (DATA_DIR / MSS_SPHEROID).read_text().splitlines()[1:]]
    changed = {header[i] for i in range(len(header)) if any(a[i] != b[i] for a, b in zip(new, old))}
    assert changed == {f"tau_crossflow_0{i}" for i in (2, 3, 5, 6)}, changed


def test_D_MSS_1_resolved_register_effect_vs_pre_fix_mss():
    """D-MSS-1 history (30_checks/README.md): the default departs from MATLAB
    at 99bf0b3 (Re on the length) on the reference cases, and the register's
    stated effect holds: REMUS Y_v|v| = -125 (diameter, the default) vs -66
    N/(m/s)^2 (length, pre-fix MSS, now a transcription only) at 0.3 m/s
    sideways (whole numbers typed from the register row)."""
    ref = _spheroid_reference(MSS_SPHEROID)
    p = _cross_flow_spheroid("midpoint")
    default = _build_cross_flow({k: v for k, v in p.items() if k != "strip_grid"})
    worst_vs_pre_fix = max(_max_diff(_evaluate(default, ref["nu_r"][k])["tau"].reshape(-1),
                                     ref["tau_crossflow"][k]) for k in range(50))
    assert worst_vs_pre_fix > G4_FACTOR * G1_TOLERANCE
    v = 0.3
    sway = np.array([0.0, v, 0.0, 0.0, 0.0, 0.0])
    y_v_diameter = _evaluate(default, sway)["tau"].reshape(-1)[1] / (v * abs(v))
    y_v_length = _mss_cross_flow(p, sway, history_reynolds_on_length=True)[1] / (v * abs(v))
    assert round(y_v_diameter) == -125, y_v_diameter
    assert round(y_v_length) == -66, y_v_length


def test_G1_MSS_lift_drag_equals_matlab():
    """forceLiftDrag.m rho = 1026 (line 26), e = 0.3 (coeffLiftDrag.m 54)."""
    ref = _spheroid_reference(MSS_SPHEROID)
    function = _build_lift_drag(_lift_drag_spheroid(_value(FLD, 26)))
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_lift_drag"][k]) <= G1_TOLERANCE, k


def test_G1_MSS_catamaran_surface_damping_equals_matlab():
    """otter.m g = 9.81 (line 90); M and G of the A-26 file (inertia about the
    combined CG)."""
    ref = _catamaran_reference(MSS_CATAMARAN)
    _, function = _build_surface(_surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN))
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], -ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), ref["tau_damp"][k]) <= G1_TOLERANCE, ("tau", k)


def test_G1_MSS_catamaran_cross_flow_default_equals_matlab():
    """Hoerner on the default grid (20 midpoints): no registered deviation."""
    ref = _catamaran_reference(MSS_CATAMARAN)
    p = {k: v for k, v in _cross_flow_catamaran("midpoint").items() if k != "strip_grid"}
    function = _build_cross_flow(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_D_MSS_2_hoerner_below_its_table_raises_naming_the_ratio():
    """D-MSS-2 (owner E-24): B/(2T) below Hoerner.m's first data point (MSS:
    interp1 -> NaN, lines 47-48) raises ValueError naming the ratio; the first
    data point itself is accepted."""
    block = _contract(CROSS_FLOW)
    first = _table(HOERNER, "CD_DATA")[0, 0]
    draft = 1.0
    beam = first * draft          # ratio = first / 2, below the table
    ratio = beam / (2 * draft)
    with pytest.raises(ValueError) as error:
        block.preprocess_cross_flow_drag(length=2.0, beam=beam, draft=draft,
                                         water_density=_value(CFD, 24), drag_model="hoerner")
    assert f"{ratio}" in str(error.value), str(error.value)
    block.preprocess_cross_flow_drag(length=2.0, beam=2 * first * draft, draft=draft,
                                     water_density=_value(CFD, 24), drag_model="hoerner")


# --------------------------------------------------------------------------
# U3a — floating linear damping: Dmtrx.m surface-craft branch
# (map row get_D_linear_surface_vessel)
# --------------------------------------------------------------------------
def test_floating_transcription_equals_matlab():
    """Second check, no block: the Dmtrx.m 51-63 transcription lands on MATLAB."""
    ref = _floating_reference()
    for k in range(len(ref["D"])):
        assert _max_diff(_mss_dmtrx_floating(_floating_parameters(k, ref)), ref["D"][k]) <= G1_TOLERANCE, k


def test_G1_numpy_source_floating_damping_equals_matlab():
    """get_D_linear_surface_vessel with coeff = 1 and zeta3 = 0.2 = Dmtrx.m."""
    damping = _source()["damping"]
    ref = _floating_reference()
    for k in range(len(ref["D"])):
        p = _floating_parameters(k, ref)
        src = damping.get_D_linear_surface_vessel(p["restoring_matrix"], p["rigid_body_mass_matrix"],
                                                  p["added_mass_matrix"], p["time_constants"],
                                                  p["damping_ratios"], np.ones(6))
        assert _max_diff(src, ref["D"][k]) <= G1_TOLERANCE, k


def test_G1_MSS_floating_damping_equals_matlab():
    """Block D and tau = -D nu_r vs MATLAB Dmtrx.m on 20 seeded craft."""
    ref = _floating_reference()
    states = _random_nu_r()[:50]
    for k in range(len(ref["D"])):
        constants, function = _build_floating(_floating_parameters(k, ref))
        assert _max_diff(constants.damping_coefficients, np.diag(ref["D"][k])) <= G1_TOLERANCE, k
        for nu_r in states:
            out = _evaluate(function, nu_r)
            assert _max_diff(out["D"], ref["D"][k]) <= G1_TOLERANCE, k
            assert _max_diff(out["tau"].reshape(-1), -ref["D"][k] @ nu_r) <= G1_TOLERANCE, k


@pytest.mark.parametrize("scaled_coeff", [False, True], ids=["coeff_ones", "coeff_mapped"])
def test_G2_floating_damping_matches_source(scaled_coeff):
    """The source's ``coeff_i`` multiplies D_ii; the contract has no coeff, so
    a coeff set is mapped onto T_i / coeff_i (DOF 1, 2, 6) and coeff_i zeta_i
    (DOF 3, 4, 5); both must give the source's D."""
    damping = _source()["damping"]
    ref = _floating_reference()
    rng = np.random.default_rng(SEED)
    for k in range(len(ref["D"])):
        p = _floating_parameters(k, ref)
        coeff = rng.uniform(0.2, 3.0, 6) if scaled_coeff else np.ones(6)
        src = damping.get_D_linear_surface_vessel(p["restoring_matrix"], p["rigid_body_mass_matrix"],
                                                  p["added_mass_matrix"], p["time_constants"],
                                                  p["damping_ratios"], coeff)
        t1, t2, t6 = p["time_constants"]
        z3, z4, z5 = p["damping_ratios"]
        mapped = {**p, "time_constants": [t1 / coeff[0], t2 / coeff[1], t6 / coeff[5]],
                  "damping_ratios": [z3 * coeff[2], z4 * coeff[3], z5 * coeff[4]]}
        _, function = _build_floating(mapped)
        for nu_r in 0.1 * _random_nu_r()[:20]:   # |tau| < 1e5 N: see G2_SCALE_NOTE
            out = _evaluate(function, nu_r)
            assert _max_diff(out["D"], src) <= G2_TOLERANCE, k
            assert _max_diff(out["tau"].reshape(-1), -src @ nu_r) <= G2_TOLERANCE, k


def test_G4_floating_roll_ratio_times_zero_is_detected():
    ref = _floating_reference()
    p = _floating_parameters(0, ref)
    assert _max_diff(_build_floating(p)[0].damping_coefficients, np.diag(ref["D"][0])) <= G1_TOLERANCE
    broken = {**p, "damping_ratios": [p["damping_ratios"][0], 0.0, p["damping_ratios"][2]]}
    assert _max_diff(_build_floating(broken)[0].damping_coefficients,
                     np.diag(ref["D"][0])) > G4_FACTOR * G1_TOLERANCE


# --------------------------------------------------------------------------
# U3a — surge damping: forceSurgeDamping.m (map row force_surge_damping)
# --------------------------------------------------------------------------
SURGE_MODELS = ("ittc", "max_thrust")
# G2_SCALE_NOTE: G2's 1e-10 is absolute; one ulp of a 7.5e5 N force is 1.2e-10,
# so two correct implementations can miss it at ship scale (U3a probe: 1.2e-10
# and 2.3e-10). The new G2 tests stay below 1e5 N; ship scale (osv.m) is held
# by G1 at 1e-9. A relative G2 for large craft needs an ADR (ledger, open item).


def test_surge_reference_set_1_is_osv():
    """Parameter set 1 of the surge CSV is osv.m's vessel (lines 62-80, 129),
    the only MSS vehicle that calls forceSurgeDamping (osv.m 180-181)."""
    ref = _named_columns(MSS_SURGE)
    row = int(np.flatnonzero(ref["set"] == 1)[0])
    length, beam, draft = _value(OSV, 62), _value(OSV, 63), _value(OSV, 64)
    rho, cb = _value(OSV, 65), _value(OSV, 66)
    k_max = _value(OSV, 70)
    assert ref["L"][row] == length and ref["rho"][row] == rho
    assert ref["m"][row] == rho * cb * length * beam * draft                  # lines 79-80
    assert ref["S"][row] == length * beam + 2 * draft * beam                  # line 67
    assert ref["T1"][row] == _value(OSV, 129) and ref["u_max"][row] == _value(OSV, 77)
    assert ref["thrust_max"][row] == k_max[2] + k_max[3]                      # line 76


def test_surge_transcription_reproduces_matlab():
    """Second check, no block: the forceSurgeDamping.m transcription lands on
    every non-probe MATLAB row (X, Xuu, Xu, A11), both branches, five sets."""
    for set_id in range(1, 6):
        for model in SURGE_MODELS:
            rows = _surge_rows(set_id, model)
            for k in range(len(rows["u_r"])):
                x, xuu, xu, a11 = _mss_force_surge_damping(
                    rows["u_r"][k], rows["m"][k], rows["S"][k], rows["L"][k], rows["T1"][k],
                    rows["rho"][k], rows["u_max"][k],
                    rows["thrust_max"][k] if model == "max_thrust" else None)
                got = np.array([x, xuu, xu, a11])
                want = np.array([rows["X"][k], rows["Xuu"][k], rows["Xu"][k], rows["A11"][k]])
                assert _max_diff(got, want) <= G1_TOLERANCE, (set_id, model, k)


def test_xuu_ittc_transcription_reproduces_matlab():
    """Second check for XuuITTC.m (MSS 97fae93, not called by forceSurgeDamping
    or any MSS vehicle): the reference for the bounded-Reynolds option the
    ledger puts to the owner (proposed D-MSS-5)."""
    ref = _named_columns(MSS_XUU_ITTC)
    for k in range(len(ref["u_r"])):
        got = _mss_xuu_ittc(ref["u_r"][k], ref["rho"][k], ref["L"][k], ref["B"][k], ref["T"][k], ref["C_B"][k])
        assert abs(got - ref["Xuu"][k]) <= G1_TOLERANCE, k


@pytest.mark.parametrize("model", SURGE_MODELS)
def test_G1_MSS_surge_damping_equals_matlab(model):
    """Block vs MATLAB forceSurgeDamping.m on five sets (osv.m + four seeded),
    40 speeds each in -8..8 m/s; constants A11, Xu (and Xuu for max_thrust)."""
    for set_id in range(1, 6):
        constants, function = _build_surge(_surge_parameters(set_id, model))
        rows = _surge_rows(set_id, model)
        assert abs(constants.added_mass - rows["A11"][0]) <= G1_TOLERANCE, set_id
        assert abs(constants.linear_coefficient - rows["Xu"][0]) <= G1_TOLERANCE, set_id
        if model == "max_thrust":
            assert abs(constants.quadratic_coefficient - rows["Xuu"][0]) <= G1_TOLERANCE, set_id
        for k, u in enumerate(rows["u_r"]):
            nu_r = np.array([u, 0.3, -0.2, 0.1, 0.05, -0.1])   # other DOFs must not enter
            tau = _evaluate(function, nu_r)["tau"].reshape(-1)
            assert _max_diff(tau, [rows["X"][k], 0, 0, 0, 0, 0]) <= G1_TOLERANCE, (set_id, k)


@pytest.fixture
def corrected_surge_source(monkeypatch):
    """The numpy source with its one sign error corrected in memory:
    ``added_mass_surge_static`` returns A11, and the source uses it as
    ``Xudot``; MSS (forceSurgeDamping.m 60) uses ``Xudot = -A11``."""
    damping = _source()["damping"]
    from more_generic_models.dynamics.plant.matrices import added_mass
    original = added_mass.AddedMass.added_mass_surge_static
    monkeypatch.setattr(added_mass.AddedMass, "added_mass_surge_static",
                        staticmethod(lambda m, L, rho: (-original(m, L, rho)[0], original(m, L, rho)[1])))
    return damping


@pytest.mark.parametrize("model", SURGE_MODELS)
def test_G2_surge_damping_matches_corrected_source(model, corrected_surge_source):
    """Sets 2-5 (|X| < 1e5 N); set 1 (osv.m, |X| ~ 1e6 N) is G1-only: see G2_SCALE_NOTE."""
    rng = np.random.default_rng(SEED)
    for set_id in range(2, 6):
        p = _surge_parameters(set_id, model)
        _, function = _build_surge(p)
        rows = _surge_rows(set_id, model)
        for u in np.concatenate([rows["u_r"], rng.uniform(-8.0, 8.0, 200)]):
            x_src = corrected_surge_source.force_surge_damping(
                u, rows["m"][0], rows["S"][0], rows["L"][0], rows["T1"][0], rows["rho"][0],
                rows["u_max"][0], rows["thrust_max"][0] if model == "max_thrust" else None)[0]
            tau = _evaluate(function, np.array([u, 0, 0, 0, 0, 0.0]))["tau"].reshape(-1)
            assert abs(tau[0] - x_src) <= G2_TOLERANCE, (set_id, u)


def test_G4_surge_time_constant_plus_1_percent_is_detected():
    rows = _surge_rows(2, "max_thrust")
    p = _surge_parameters(2, "max_thrust")
    worst = lambda q: max(abs(_evaluate(_build_surge(q)[1], np.array([u, 0, 0, 0, 0, 0.0]))["tau"][0, 0] - x)
                          for u, x in zip(rows["u_r"], rows["X"]))
    assert worst(p) <= G1_TOLERANCE
    assert worst({**p, "time_constant": 1.01 * p["time_constant"]}) > G4_FACTOR * G1_TOLERANCE


def test_numpy_source_surge_linear_coefficient_sign_is_documented():
    """Documents a source sign error (not fixed here: more_generic_models
    physics is read-only for this job): the source's Xu is -(m - A11)/T1,
    MSS's and Fossen (2011) eq. 6.76's is -(m + A11)/T1. When the source is
    fixed this test fails and is replaced by plain G2."""
    damping = _source()["damping"]
    for set_id in range(1, 6):
        rows = _surge_rows(set_id, "ittc")
        m, a11, t1 = rows["m"][0], rows["A11"][0], rows["T1"][0]
        xu_src = damping.force_surge_damping(1.0, m, rows["S"][0], rows["L"][0], t1, rows["rho"][0],
                                             rows["u_max"][0])[2]
        assert abs(xu_src - (-(m - a11) / t1)) <= G2_TOLERANCE * abs(xu_src), set_id
        assert abs(rows["Xu"][0] - (-(m + a11) / t1)) <= G2_TOLERANCE * abs(xu_src), set_id


# Physical-sign tests (owner E-23: MSS is the reference, not the truth)
def test_physical_linear_surge_coefficient_gives_the_time_constant():
    """Fossen (2011) eqs. 6.71 and 6.76 (pp. 124-125, read in full): the
    linearised surge mode (m + A11) du/dt = Xu u decays with time constant
    T1 = (m + A11) / (-Xu). Model-free on the MATLAB rows (MSS satisfies it),
    then on the block's own derivative dX/du at u = 0 (sigma(0) = 1)."""
    import casadi as ca
    for set_id in range(1, 6):
        rows = _surge_rows(set_id, "ittc")
        m, a11, t1 = rows["m"][0], rows["A11"][0], rows["T1"][0]
        assert abs((m + a11) / -rows["Xu"][0] - t1) <= 1e-12 * t1, set_id
        for model in SURGE_MODELS:
            _, function = _build_surge(_surge_parameters(set_id, model))
            x = ca.SX.sym("x", 6)
            slope = float(ca.Function("j", [x], [ca.jacobian(function(nu_r=x)["tau"][0], x)])(np.zeros(6))[0])
            assert abs((m + a11) / -slope - t1) <= 1e-12 * t1, (set_id, model, slope)


def test_physical_surge_damping_dissipation_mss_blend_is_documented():
    """Fossen (2011) Property 6.3 (p. 123): damping is dissipative, X u_r <= 0.
    MATLAB's own rows: true for every u_r >= 0; false in reverse, because
    sigma = 1 - tanh(u_r/u_cross) (forceSurgeDamping.m 79) exceeds 1 for
    u_r < 0 and flips the quadratic term. Example pinned: osv.m with its
    thrust calibration at u_r = -8 m/s gives X < 0 (pushes it further astern).
    Owner decision pending (ledger, proposed D-MSS-4)."""
    ref = _named_columns(MSS_SURGE)
    real = ref["probe"] == 0
    power = ref["X"] * ref["u_r"]
    assert np.all(power[real & (ref["u_r"] >= 0)] <= 0.0)
    assert np.any(power[real & (ref["u_r"] < 0)] > 0.0)
    osv_reverse = real & (ref["set"] == 1) & (ref["branch"] == 1) & (ref["u_r"] == -8.0)
    assert osv_reverse.sum() == 1 and ref["X"][osv_reverse][0] < 0.0


def test_physical_surge_ittc_branch_singularity_is_documented():
    """Fossen (2011) p. 125 below eq. 6.85: a minimum Rn must be used or C_F
    blows up at low speed. forceSurgeDamping.m 72-73 has none (eps only), so at
    Rn = 100 (u_r = 1e-4/L) MATLAB returns |X| > 1e9 N on osv.m; the
    max_thrust branch at the same speeds stays below 1 N. XuuITTC.m bounds
    Re at 1e5 (line 34). Owner decision pending (ledger, proposed D-MSS-5)."""
    ittc = _surge_rows(1, "ittc", probe=1)
    thrust = _surge_rows(1, "max_thrust", probe=1)
    assert np.max(np.abs(ittc["X"])) > 1e9
    assert np.max(np.abs(thrust["X"])) < 1.0


def test_physical_cross_flow_equals_first_principles_strip_integral():
    """The strip at x (from the CO, z down) moves with v + x r sideways and
    w - x q vertically (omega x r for r = [x 0 0]); its moments are
    N = +x Y_i and M = -x Z_i. Fossen (2011) eqs. 6.91-6.92 (p. 127) state the
    sway/yaw half; heave/pitch is this ledger's derivation. crossFlowDrag.m
    47 and 50 use w + x q and +x Z: two sign flips that cancel on a grid
    symmetric about the CO (A-22 finding 6), so the default block must equal
    the first-principles form; a pure pitch or yaw rate must be opposed."""
    for p in (_cross_flow_spheroid("midpoint"), _cross_flow_catamaran("midpoint")):
        function = _build_cross_flow({k: v for k, v in p.items() if k != "strip_grid"})
        n_strips = int(_value(CFD, 25))
        dx = p["length"] / n_strips
        xs = -p["length"] / 2 + (np.arange(1, n_strips + 1) - 0.5) * dx
        for k, nu_r in enumerate(_states()):
            cd = (_mss_cylinder_cd(p["length"], p["beam"], nu_r) if p["drag_model"] == "cylinder"
                  else _mss_hoerner(p["beam"], p["draft"]))
            v_pt, w_pt = nu_r[1] + xs * nu_r[5], nu_r[2] - xs * nu_r[4]
            per_strip = -0.5 * p["water_density"] * p["draft"] * cd * dx
            y_i, z_i = per_strip * np.abs(v_pt) * v_pt, per_strip * np.abs(w_pt) * w_pt
            physical = np.array([0, y_i.sum(), z_i.sum(), 0, -(xs * z_i).sum(), (xs * y_i).sum()])
            out = _evaluate(function, nu_r)["tau"].reshape(-1)
            assert _max_diff(out, physical) <= G1_TOLERANCE, (p["drag_model"], k)
        pitch = _evaluate(function, np.array([0, 0, 0, 0, 0.2, 0.0]))["tau"].reshape(-1)
        yaw = _evaluate(function, np.array([0, 0, 0, 0, 0, 0.2]))["tau"].reshape(-1)
        assert pitch[4] < 0 and abs(pitch[2]) <= G1_TOLERANCE
        assert yaw[5] < 0 and abs(yaw[1]) <= G1_TOLERANCE


# --------------------------------------------------------------------------
# A-27 finding 4: Jacobians at and near rest (documents today; owner decides)
# --------------------------------------------------------------------------
JACOBIAN_STATES = {
    "rest": [0, 0, 0, 0, 0, 0],
    "pure_sway": [0, 0.3, 0, 0, 0, 0],             # sideways current, no way on
    "pure_rotation": [0, 0, 0, 0.1, 0.1, 0.1],     # turning on the spot
}
# (block, output, state) -> non-finite Jacobian entries (row of the output,
# column of nu_r), probed 2026-10-06. Submerged damping: U_r = |nu_r[0:3]|
# (remus100.m 128) inside exp(-3 U_r) has d/dnu = nu/U_r = 0/0 when u=v=w=0.
# Lift/drag: atan2(w_r, u_r) (remus100.m 127) has 0/0 derivatives when
# u_r = w_r = 0, which includes pure sway.
JACOBIAN_NONFINITE_TODAY = {
    ("submerged", "tau", "rest"): {(0, 0), (0, 1), (0, 2)},
    ("submerged", "tau", "pure_rotation"): {(0, 0), (0, 1), (0, 2)},
    ("lift_drag", "tau", "rest"): {(0, 0), (0, 2), (2, 0), (2, 2)},
    ("lift_drag", "tau", "pure_sway"): {(0, 0), (0, 2), (2, 0), (2, 2)},
    ("lift_drag", "tau", "pure_rotation"): {(0, 0), (0, 2), (2, 0), (2, 2)},
}


def _jacobian_blocks():
    return {
        "submerged": _build_submerged({k: v for k, v in _spheroid_damping(False, MSS_SPHEROID).items()
                                       if k != "sway_damping_fade"})[1],
        "surface": _build_surface(_surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN))[1],
        "cross_flow_cylinder": _build_cross_flow({k: v for k, v in _cross_flow_spheroid("midpoint").items()
                                                  if k != "strip_grid"}),
        "cross_flow_hoerner": _build_cross_flow({k: v for k, v in _cross_flow_catamaran("midpoint").items()
                                                 if k != "strip_grid"}),
        "lift_drag": _build_lift_drag(_lift_drag_spheroid(_value(FLD, 26))),
    }


def _nonfinite_jacobian(function, output, state):
    import casadi as ca
    x = ca.SX.sym("x", 6)
    jac = np.array(ca.Function("j", [x], [ca.jacobian(function(nu_r=x)[output], x)])(np.asarray(state, float)))
    return {(int(i), int(j)) for i, j in zip(*np.nonzero(~np.isfinite(jac)))}


@pytest.mark.parametrize("state", sorted(JACOBIAN_STATES))
def test_A27_finding4_nonfinite_jacobian_entries_are_documented(state):
    """Documents today's behaviour for the owner's decision (ledger, options
    a/b/c): exactly the listed entries are non-finite; cross-flow and surface
    damping are finite; 1e-12 off the singular set every Jacobian is finite."""
    for name, function in _jacobian_blocks().items():
        expected = JACOBIAN_NONFINITE_TODAY.get((name, "tau", state), set())
        assert _nonfinite_jacobian(function, "tau", JACOBIAN_STATES[state]) == expected, (name, state)
        assert not _nonfinite_jacobian(function, "tau", np.full(6, 1e-12)), name


@pytest.mark.parametrize("model", SURGE_MODELS)
def test_surge_damping_jacobian_is_finite_at_rest(model):
    """|u_r| u_r and tanh are smooth at 0; the ITTC log10(Rn + eps) enters
    multiplied by |u_r| u_r (zero slope), so the Jacobian must be finite at
    rest and on the other two singular sets of finding 4."""
    _, function = _build_surge(_surge_parameters(2, model))
    for state in JACOBIAN_STATES.values():
        assert not _nonfinite_jacobian(function, "tau", state), (model, state)


# --------------------------------------------------------------------------
# E-25: transforms only from L0 (more_transformations)
# --------------------------------------------------------------------------
OUR_HYDRODYNAMICS_MODULES = (LINEAR_DAMPING, CROSS_FLOW, LIFT_DRAG, SURGE_DAMPING)  # not Luka's linear_surface (E-28)
ALLOWED_IMPORT_ROOTS = {"casadi", "numpy", "dataclasses", "typing", "math", "__future__",
                        "more_transformations", "more_dynamics"}
TRANSFORM_NAME = re.compile(r"(skew|smtrx|rzyx|rotation|rot_|hmtrx|euler|quat|gravity|jacobian_|tzyx)")


@pytest.mark.parametrize("module_name", OUR_HYDRODYNAMICS_MODULES)
def test_E25_no_local_transform_in_hydrodynamics(module_name):
    """No function defined here looks like a rotation, skew, H, T, J or
    gravity, and nothing is imported from outside the allowed roots; any
    transform a block needs comes from ``more_transformations`` (numpy) or
    ``more_transformations.more_casadi_transformations`` (CasADi)."""
    tree = ast.parse(Path(_contract(module_name).__file__).read_text())
    defined = [n.name for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
    assert not [n for n in defined if TRANSFORM_NAME.search(n.lower())], defined
    roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0:
            roots.add(node.module.split(".")[0])
    assert roots <= ALLOWED_IMPORT_ROOTS, roots - ALLOWED_IMPORT_ROOTS


# --------------------------------------------------------------------------
# Generic by construction (owner E-16)
# --------------------------------------------------------------------------
def _code_without_docstrings(module):
    tree = ast.parse(Path(module.__file__).read_text())
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(body, list) and body and isinstance(body[0], ast.Expr) \
                and isinstance(getattr(body[0], "value", None), ast.Constant) \
                and isinstance(body[0].value.value, str):
            body.pop(0)
    return ast.unparse(tree).lower()


GENERIC = {
    LINEAR_DAMPING: {"preprocess_submerged_linear_damping": {"sway_damping_fade"},
                     "preprocess_surface_linear_damping": set()},
    CROSS_FLOW: {"preprocess_cross_flow_drag": {"strip_grid"}},   # U3a: length flag gone
    LIFT_DRAG: {"preprocess_lift_drag": set()},
    SURGE_DAMPING: {"preprocess_surge_damping_ittc": set(),
                    "preprocess_surge_damping_max_thrust": set()},
}
GENERIC[LINEAR_DAMPING]["preprocess_floating_linear_damping"] = set()


@pytest.mark.parametrize("module_name", sorted(GENERIC))
def test_generic_block_reads_no_vehicle_name_and_has_no_vehicle_defaults(module_name):
    block = _contract(module_name)
    code = _code_without_docstrings(block)
    assert not [n for n in VEHICLE_NAMES if n in code]
    for function_name, flags in GENERIC[module_name].items():
        signature = inspect.signature(_attribute(module_name, function_name))
        defaults = {n for n, p in signature.parameters.items() if p.default is not inspect.Parameter.empty}
        assert defaults == flags, (function_name, defaults)
