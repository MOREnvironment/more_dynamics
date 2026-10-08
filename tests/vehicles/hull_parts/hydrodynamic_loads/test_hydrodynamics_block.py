"""Gate tests for the hydrodynamics blocks: linear damping, cross-flow drag, lift/drag.

Written 2026-10-05, before the blocks existed. Gates (the test names carry
them): G1 against MATLAB running MSS (frozen CSVs), 1e-9 absolute; G2 against
an independent transcription of the cited MSS lines, 1e-10; G4 a perturbed model is detected; G5 a printed number reproduced.
Section "G1-MSS" compares the default paths with MATLAB running current MSS
(2026-10-05, ``*_mss_current.csv``); the test-side transcriptions stay as the
second check and are themselves checked against those files. Hoerner below
its table raises ``ValueError`` (MSS returns NaN there; owner's decision of
2026-10-05).

Added 2026-10-06:

* MSS ``ac77394`` (MSS issue #81) puts the cross-flow Reynolds number on the
  diameter, as DNV-RP-C205 defines it, so the default path is G1 against
  MATLAB at ``ac77394`` (``spheroid_matlab_reference_mss_ac77394.csv``) and
  the ``cross_flow_reynolds_length`` flag is gone (owner, 2026-10-06;
  ``test_cross_flow_length_flag_is_gone``). The ``99bf0b3`` file stays for
  history.
* ``Dmtrx.m``'s surface-craft branch is tested (a third function in
  ``linear_damping``), and ``force_surge_damping`` is ``forceSurgeDamping.m``
  (module ``surge_damping``); both G1 against MATLAB calling those MSS
  functions.
* Physical-sign tests (MSS is the reference, not the truth): the linear surge
  coefficient reproduces its time constant; MSS's surge blend is not dissipative in reverse and its ITTC branch
  blows up at Rn = 100 (documented; our defaults are the symmetric blend and
  the Rn floor, owner's decision of 2026-10-06, MSS behind the flags, both
  tested); cross-flow pitch/yaw terms equal the first-principles strip
  integral.
* The Jacobian entries that are non-finite at rest with the MSS lines are
  pinned; ``smooth_speed_epsilon > 0`` removes them and keeps lift/drag
  dissipative (2026-10-07).
* No local transform in our hydrodynamics modules.

Rules these tests encode (the owner's decisions, 2026-10-05/06)
---------------------------------------------------------------
* The **default** path equals current MSS; each template setting stays
  behind a named flag that reproduces the template value.
* Cross-flow drag on MSS's 20 strip midpoints with the full-precision
  ``cylinderDrag.m`` table (default); the 21-end-point grid behind a flag.
* The cross-flow Reynolds number uses the **diameter**, which is newest MSS
  since ``ac77394`` (``cylinderDrag.m`` 79-80); no flag.
* Generic by construction (no vehicle name or number inside a block).
* The other developer's ``linear_surface.py`` and its tests are not touched.

MSS revisions: every MSS line below is cited at ``72656d1`` (release 2.0.2),
except ``forceSurgeDamping.m`` and the two ``osv.m`` lines that call it,
cited at ``ac77394``, the last revision that holds that file (deleted in MSS
``108ceda``, release 2.0, 2026-10-06).

Contract of the block
---------------------
All-CasADi blocks with named parameters. Each module declares its numbers
(``Parameter``: name, shape, SI unit, meaning, range) and builds one
``ca.Function`` per choice of its selectors; matrices that come from other
blocks are named coupling inputs. The tests bind the checked numbers and the
couplings into a function of ``nu_r`` alone (``_bind``, the frozen path a
plugin carries) and keep every assertion of the earlier contract.

``more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.damping.linear_damping``:

* ``linear_damping_parameters(form, *, smooth_speed=False)``,
  ``linear_damping_couplings(form)``, ``linear_damping_casadi(*, form,
  sway_damping_fade=False, smooth_speed=False)`` with ``form`` in
  ``("submerged", "floating", "surface")``, and
  ``check_linear_damping_values(values, *, form, smooth_speed=False)``.
  ``smooth_speed`` is a selector (submerged form only); with ``True`` the
  regularisation ``smooth_speed_epsilon`` (m/s, ``> 0``) is a declared
  parameter and a named input, the last one.
* submerged: couplings ``rigid_body_mass_matrix``, ``added_mass_matrix``;
  parameters ``weight`` (N), ``center_of_gravity``, ``center_of_buoyancy``
  ((3,) m from the CO), ``time_constants = (T1, T2, T6)`` s,
  ``damping_ratios = (zeta4, zeta5)``. Outputs ``D`` (6x6, positive diagonal
  as ``Dmtrx.m``) with ``D[0,0] *= exp(-3 U_r)`` (``remus100.m:218``),
  ``U_r = |nu_r[0:3]|``, and ``D[1,1]`` faded the same way only when
  ``sway_damping_fade`` is True (the template generator's setting);
  ``tau = -D @ nu_r``; ``damping_coefficients`` (the diagonal before any
  fade).
* floating (``Dmtrx.m`` 51-63, the surface-craft branch): couplings
  ``rigid_body_mass_matrix``, ``added_mass_matrix``, ``restoring_matrix``
  (only G33, G44, G55 used); ``time_constants = (T1, T2, T6)``,
  ``damping_ratios = (zeta3, zeta4, zeta5)``; MSS fixes ``zeta3 = 0.2``
  (``Dmtrx.m`` 57), so it is a parameter and the MSS tests pass 0.2. The
  source's ``coeff`` multipliers are not in the contract: ``coeff_i`` scales
  ``D_ii``, which a parameter set gets by scaling ``T_i`` or ``zeta_i``.
  Outputs ``D``, ``tau = -D @ nu_r``, ``damping_coefficients``.
* surface (``otter.m``): couplings ``mass_matrix`` (M_RB + M_A),
  ``restoring_matrix`` (G about the centre of flotation); parameters
  ``max_forward_thrust`` N, ``max_speed`` m/s, ``time_constants = (T_sway,
  T_yaw)``, ``damping_ratios = (zeta3, zeta4, zeta5)``,
  ``yaw_damping_nonlinearity`` (the 10 of ``otter.m:240``). Outputs
  ``D = -diag(damping_derivatives)`` (positive), ``tau`` = ``otter.m``
  235-240 (the force on the vehicle, quadratic yaw term included) and
  ``damping_derivatives = [Xu Yv Zw Kp Mq Nr]`` (``otter.m`` 202-207).

``more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.surge_resistance.surge_damping`` (``forceSurgeDamping.m`` at ``ac77394``):

* ``surge_damping_parameters(surge_form, *, ittc_reynolds_bound="floor")``,
  ``surge_damping_couplings(surge_form)``, ``surge_damping_casadi(*,
  surge_form, surge_blend="symmetric", ittc_reynolds_bound="floor")``,
  ``check_surge_damping_values(values, *, surge_form,
  ittc_reynolds_bound="floor")``; ``surge_form`` in ``("ittc", "max_thrust",
  "exp_ittc")``. On the ITTC forms with ``"floor"`` the Reynolds floor
  ``ittc_reynolds_floor`` (``> 100``) is a declared parameter and a named
  input, the last one; the tests give it ``XuuITTC.m`` 34's ``1e5``.
* ``"ittc"`` (``forceSurgeDamping.m`` 68-75): couplings ``mass``,
  ``wetted_surface``; parameters ``length``, ``water_density``,
  ``time_constant``, ``form_factor`` (the ``k = 0.1`` of line 70),
  ``crossover_speed`` (the ``u_cross = 2`` of line 57); ``Xuu`` from the
  ITTC-1957 line at the current speed. ``"max_thrust"`` (lines 64-66):
  coupling ``mass``; ``length``, ``water_density``, ``time_constant``,
  ``max_speed``, ``max_thrust``, ``crossover_speed``;
  ``Xuu = -max_thrust / max_speed^2``.
* The two selectors (owner's decision of 2026-10-06): ``surge_blend="symmetric"``
  (default) blends on ``|u_r|``, ``"mss_tanh"`` is MSS line 79;
  ``ittc_reynolds_bound="floor"`` (default) bounds ``Rn`` below by the
  parameter ``ittc_reynolds_floor`` as ``XuuITTC.m`` 34-36,
  ``"mss_offset"`` is MSS lines 72-73 (unbounded). The test parameter dicts
  keep the earlier keys (``ittc_reynolds_floor=None`` for the MSS line,
  ``smooth_speed_epsilon=0`` for MSS) because the committed blocks read them;
  ``_selectors_from_numbers`` turns them into the selectors and parameters.
* ``A11 = 2.7 rho (m/rho)^(5/3) / L^2`` (``addedMassSurge.m`` 32-33),
  ``Xudot = -A11`` (``forceSurgeDamping.m`` line 60), ``Xu = -(m - Xudot) / T1 =
  -(m + A11) / T1`` (line 61). Outputs ``tau = [X 0 0 0 0 0]`` with
  ``sigma = 1 - tanh(u_r / u_cross)`` and ``X = sigma Xu u_r + (1 - sigma)
  Xuu |u_r| u_r`` (lines 79-82), ``added_mass`` (A11 > 0),
  ``linear_coefficient`` (Xu < 0), ``quadratic_coefficient`` (Xuu < 0).
  The force on the vehicle (added, as ``osv.m`` 180-182 at ``ac77394``).
  The vehicle that uses it must set its linear ``D[0,0] = 0`` (``osv.m`` 188).
* ``"exp_ittc"`` (the MSS 2.0.2 ship law, ``osv.m`` 182-185 at ``cc07579``):
  couplings ``mass_matrix`` (M(1,1)), ``wetted_surface``; parameters
  ``length``, ``water_density``, ``time_constant``, ``form_factor``,
  ``exponential_decay_rate`` (``k_u``); ``D_nl,11 = exp(-k_u |u_r|) D11 -
  Xuu |u_r|``, ``D11 = M(1,1) / T1``, ``X = -D_nl,11 u_r``; outputs also
  ``surge_damping_coefficient`` (``D_nl,11``).

``more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.cross_flow``:

* ``cross_flow_drag_parameters()`` (``length``, ``beam``, ``draft``,
  ``water_density``), ``cross_flow_drag_casadi(*, drag_model,
  strip_grid="midpoint")``; ``drag_model`` in {"cylinder", "hoerner"} (no
  default), ``strip_grid`` in {"midpoint", "endpoint"};
  ``check_cross_flow_drag_values(values, *, drag_model)`` refuses a Hoerner
  ratio below the table. The Reynolds number is on the ``beam``, the
  diameter (MSS calls ``crossFlowDrag(L, D, D, ...)``; ``cylinderDrag.m``
  79-80 at ``ac77394``); ``cross_flow_reynolds_length`` is gone (owner
  2026-10-06). Output ``["tau"]``.

``more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.lift_drag.lift_drag``:

* ``lift_drag_parameters()`` (``span``, ``planform_area``,
  ``parasitic_drag_coefficient``, ``oswald_efficiency``, ``water_density``),
  ``lift_drag_casadi(*, smooth_speed=False)``; with ``smooth_speed=True``
  also ``smooth_speed_epsilon`` (m/s, ``> 0``); output ``["tau"]``;
  ``alpha = atan2(w_r, u_r)``, ``U_r = |nu_r[0:3]|`` (``remus100.m`` 126-127).

Conventions (hidden assumptions)
--------------------------------
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
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[3] / "data" / "hydrodynamics"
# Outside this repo only through environment variables (nothing relative to one machine):
# MSS_DIR = an MSS checkout. Default: not set. The gates read the frozen CSVs and
# the snapshot of every cited MSS / generator line (``cited_lines_snapshot.json``,
# made by ``snapshot_cited_lines.py``); only the live re-read of the MSS lines
# needs the variable, and it skips without it.
CITED_SNAPSHOT = DATA_DIR / "cited_lines_snapshot.json"
LINEAR_DAMPING = "more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.damping.linear_damping"
CROSS_FLOW = "more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.cross_flow"
LIFT_DRAG = "more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.lift_drag.lift_drag"
SURGE_DAMPING = "more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.surge_resistance.surge_damping"

G1_TOLERANCE = 1e-9   # G1: block vs MATLAB running MSS, absolute
G4_FACTOR = 10.0      # G4: a perturbed model must differ by more than this x G1
SEED = 20261005
N_RANDOM_STATES = 1000
# The source's rounded Re table moves the spheroid cross-flow force by 2.3e-2 N (2026-10-05)

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
# The surge reference was made with MSS ac77394, the last revision that holds
# forceSurgeDamping.m (deleted in MSS release 2.0.2, 72656d1): a cited path
# "$MSS_DIR@<revision>/..." is read from that revision of the MSS_DIR git
# checkout (git show); every other MSS line is cited at 72656d1.
FSD = "$MSS_DIR@ac77394/LIBRARY/modeling/forceSurgeDamping.m"
AMS = "$MSS_DIR/LIBRARY/modeling/addedMassSurge.m"
XITTC = "$MSS_DIR/LIBRARY/modeling/XuuITTC.m"
OSV = "$MSS_DIR/CRAFT/SHIP/models/osv.m"
OSV_SURGE_CALL = "$MSS_DIR@ac77394/CRAFT/SHIP/models/osv.m"
# osv.m 188 at 72656d1 (release 2.0.2) zeroed the linear surge coefficient after
# reading it; MSS removed the line in de9a316 (2026-10-07), so it is cited at 72656d1.
OSV_2_0_2 = "$MSS_DIR@72656d1/CRAFT/SHIP/models/osv.m"
# The template generators of the two template references, in this repository
# ("$TESTS/..." = this tests/ folder; no variable needed).
GEN_SPH = "$TESTS/data/hydrodynamics/generate_spheroid_template.m"
GEN_CAT = "$TESTS/data/hydrodynamics/generate_catamaran_template.m"

# (file relative to <more>, line) -> how the stripped line starts. Pinned by
# test_cited_lines_are_unchanged; every constant below is parsed from these.
CITED_LINES = {
    (REMUS, 96): "mu = deg2rad(63.446827);",
    (REMUS, 126): "alpha = atan2( nu_r(3), nu_r(1) );",
    (REMUS, 127): "U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );",
    (REMUS, 131): "L_auv = 1.6;",
    (REMUS, 132): "D_auv = 0.19;",
    (REMUS, 133): "S = 0.7 * L_auv * D_auv;",
    (REMUS, 135): "b = 1.0096 * D_auv/2;",
    (REMUS, 137): "r_bG = [ 0 0 0.02 ]';",
    (REMUS, 138): "r_bB = [ 0 0 0 ]';",
    (REMUS, 143): "Cd = 0.42;",
    (REMUS, 144): "CD_0 = Cd * pi * b^2 / S;",
    (REMUS, 192): "T1 = 20;",
    (REMUS, 193): "T2 = 20;",
    (REMUS, 194): "zeta4 = 0.3;",
    (REMUS, 195): "zeta5 = 0.8;",
    (REMUS, 196): "T6 = 1;",
    (REMUS, 217): "D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bG' r_bB']);",
    (REMUS, 218): "D(1,1) = D(1,1) * exp(-3 * U_r);",
    (REMUS, 220): "tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);",
    (REMUS, 221): "tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r,'cylinder');",
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
    (CFD, 36): "rho = 1025;",
    (CFD, 37): "nStrips = 20;",
    (CFD, 38): "dx = L / nStrips;",
    (CFD, 56): "xL = -L/2 + (i - 0.5) * dx;",
    (CFD, 61): "U_h = abs(v_r + xL * r) * (v_r + xL * r);",
    (CFD, 62): "U_v = abs(w_r + xL * q) * (w_r + xL * q);",
    (CFD, 63): "Yh = Yh - 0.5 * rho * T * Cd_2D * U_h * dx;",
    (CFD, 64): "Zh = Zh - 0.5 * rho * T * Cd_2D * U_v * dx;",
    (CFD, 65): "Mh = Mh - 0.5 * rho * T * Cd_2D * xL * U_v * dx;",
    (CFD, 66): "Nh = Nh - 0.5 * rho * T * Cd_2D * xL * U_h * dx;",
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
    (AMS, 32): "nabla = m / rho;",
    (AMS, 33): "A11 = 2.7 * rho * nabla^(5/3) / L^2;",
    (XITTC, 29): "C_B = 0.65;",
    (XITTC, 32): "nu_kin = 1e-6;",
    (XITTC, 33): "k = 0.1;",
    (XITTC, 34): "Re_min = 1e5;",
    (XITTC, 35): "Re = max(L * abs(u_r) / nu_kin,Re_min);",
    (XITTC, 36): "Cf = 0.075 / (log10(Re) - 2)^2;",
    (XITTC, 38): "S = 1.025 * L * (C_B*B + 1.7*T);",
    (XITTC, 39): "Xuu = -0.5 * rho * S * (1+k) * Cf;",
    (OSV, 61): "vessel.L = 83;",
    (OSV, 62): "vessel.B = 18;",
    (OSV, 63): "vessel.T = 5;",
    (OSV, 64): "vessel.rho = 1025;",
    (OSV, 65): "vessel.Cb = 0.65;",
    (OSV, 66): "vessel.S = vessel.L * vessel.B + 2 * vessel.T * vessel.B;",
    (OSV, 69): "vessel.K_max = [300e3 300e3 420e3 655e3]';",
    (OSV, 75): "vessel.thrust_max = vessel.K_max(3)+vessel.K_max(4);",
    (OSV, 76): "vessel.U_max = 7.7;",
    (OSV, 78): "vessel.nabla = vessel.Cb * vessel.L * vessel.B * vessel.T;",
    (OSV, 79): "vessel.m = vessel.rho * vessel.nabla;",
    (OSV, 128): "vessel.T1 = 100;",
    (OSV_SURGE_CALL, 180): "[X,Xuu,Xu] = forceSurgeDamping(flag,nu_r(1),vessel.m,vessel.S,vessel.L, ...",
    (OSV_SURGE_CALL, 181): "vessel.T1,vessel.rho,vessel.U_max,vessel.thrust_max);",
    (OSV_2_0_2, 188): "vessel.D(1,1) = 0;",
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
    (GEN_SPH, 231): "rho = 1025;",
    (GEN_SPH, 292): "D(1,1) = D(1,1) * exp(-3*U_r);",
    (GEN_SPH, 293): "D(2,2) = D(2,2) * exp(-3*U_r);",
    (GEN_SPH, 492): "rho = 1025;",
    (GEN_CAT, 35): "mp = 25;",
    (GEN_CAT, 165): "mu = deg2rad(63.446827);",
    (GEN_CAT, 166): "g = gravity(mu);",
    (GEN_CAT, 216): "Xu = -24.4*g/Umax; Yv = -M(2,2)/T_sway; Zw = -2*0.3*w3*M(3,3);",
}


# --------------------------------------------------------------------------
# Helpers: pinned text, contract, data
# --------------------------------------------------------------------------
def _env_dir(variable):
    value = os.environ.get(variable)
    if not value:
        pytest.skip(f"{variable} is not set (nothing relative to one machine); set it to run this check")
    path = Path(value).expanduser()
    if not path.is_dir():
        pytest.skip(f"{variable}={value} is not a directory")
    return path


def _live_text(ref):
    """``$VARIABLE/relative/path`` -> the text of a file under that environment
    variable; ``$VARIABLE@<revision>/relative/path`` -> the file at that
    revision of the git checkout the variable names (``git show``)."""
    root, rel = ref[1:].split("/", 1)
    variable, _, revision = root.partition("@")
    if variable == "TESTS":   # a file of this repository
        return (Path(__file__).resolve().parents[3] / rel).read_text()
    directory = _env_dir(variable)
    if revision:
        out = subprocess.run(["git", "-C", str(directory), "show", f"{revision}:{rel}"],
                             capture_output=True, text=True)
        if out.returncode != 0:
            pytest.skip(f"{rel} at revision {revision} not readable from {variable}: {out.stderr.strip()}")
        return out.stdout
    path = directory / rel
    if not path.exists():
        pytest.skip(f"{rel} not found under {variable}")
    return path.read_text()


def _snapshot():
    return json.loads(CITED_SNAPSHOT.read_text())


def _line(rel, number):
    """A cited line, stripped, from the frozen snapshot; it must be pinned."""
    assert (rel, number) in CITED_LINES, ("cite the line first", rel, number)
    line = _snapshot()["lines"][rel][str(number)]
    assert line.startswith(CITED_LINES[(rel, number)]), (rel, number, line)
    return line


def _live_line(rel, number):
    return _live_text(rel).splitlines()[number - 1].strip()


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
            pytest.fail(f"casadi is not installed: {exc}")
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
# MATLAB running current MSS (2026-10-05, SOURCE.md): the default-path references
MSS_SPHEROID = "spheroid_matlab_reference_mss_current.csv"      # MSS 99bf0b3
MSS_CATAMARAN = "catamaran_matlab_reference_mss_current.csv"
# MSS cc07579: otter.m KB with the waterplane area (cf349d4); G, Kp, Mq change
MSS_CATAMARAN_CC07579 = "catamaran_matlab_reference_mss_cc07579.csv"
# osv.m surge law of release 2.0.2 at cc07579 (generate_osv_surge_mss.m)
MSS_OSV_SURGE = "osv_surge_damping_mss_cc07579.csv"
# 2026-10-06: the same generator re-run unchanged at MSS ac77394; only the four
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
    branch the 50 reference cases never reach)."""
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
    l_auv, d_auv = _value(REMUS, 131), _value(REMUS, 132)
    s = _value(REMUS, 133, {"L_auv": l_auv, "D_auv": d_auv})
    b = _value(REMUS, 135, {"D_auv": d_auv})
    cd0 = _value(REMUS, 144, {"Cd": _value(REMUS, 143), "b": b, "S": s})
    return {"L": l_auv, "D": d_auv, "S": s, "CD_0": cd0}


def _spheroid_damping(fade, csv=LEGACY_SPHEROID):
    """Matrices and weight as the spheroid reference generator built them."""
    ref = _spheroid_reference(csv)
    assert np.ptp(ref["M_RB"], axis=0).max() == 0.0 and np.ptp(ref["M_A"], axis=0).max() == 0.0
    mrb, ma = ref["M_RB"][0], ref["M_A"][0]
    return {
        "rigid_body_mass_matrix": mrb,
        "added_mass_matrix": ma,
        "weight": mrb[0, 0] * _gravity(_value(REMUS, 96)),  # remus100.m 214
        "center_of_gravity": _value(REMUS, 137),
        "center_of_buoyancy": _value(REMUS, 138),
        "time_constants": [_value(REMUS, 192), _value(REMUS, 193), _value(REMUS, 196)],
        "damping_ratios": [_value(REMUS, 194), _value(REMUS, 195)],
        "sway_damping_fade": fade,
    }


def _cross_flow_spheroid(grid):
    geo = _spheroid_geometry()
    return {"length": geo["L"], "beam": geo["D"], "draft": geo["D"],   # remus100.m 221
            "water_density": _value(CFD, 36), "drag_model": "cylinder",
            "strip_grid": grid}


def _catamaran_draft():
    m, rho, mp = _value(OTTER, 94), _value(OTTER, 91), _value(GEN_CAT, 35)
    nabla = _value(OTTER, 121, {"m": m, "mp": mp, "rho": rho})
    return _value(OTTER, 122, {"nabla": nabla, "Cb_pont": _value(OTTER, 107),
                               "B_pont": _value(OTTER, 104), "L": _value(OTTER, 92)})


def _cross_flow_catamaran(grid):
    return {"length": _value(OTTER, 92), "beam": _value(OTTER, 104),
            "draft": _catamaran_draft(), "water_density": _value(CFD, 36),
            "drag_model": "hoerner", "strip_grid": grid}                # otter.m 244


def _lift_drag_spheroid(rho):
    geo = _spheroid_geometry()
    return {"span": geo["D"], "planform_area": geo["S"],                # remus100.m 220
            "parasitic_drag_coefficient": geo["CD_0"],
            "oswald_efficiency": _value(CLD, 54), "water_density": rho}


def _surface_damping_catamaran(gravity, csv=LEGACY_CATAMARAN):
    """Mass and restoring matrices of the catamaran reference; G moved from
    the CO to the CF with H(-r_f), r_f = [LCF 0 0] (otter.m 191-193)."""
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
    return _gravity(_value(GEN_CAT, 165))  # generator lines 136-137


# --------------------------------------------------------------------------
# MSS transcriptions (each line cited and pinned)
# --------------------------------------------------------------------------
def _mss_dmtrx(p, nu_r):
    """Dmtrx.m 30, 44-49 (submerged branch) + remus100.m 218 (surge fade only)."""
    m = p["rigid_body_mass_matrix"] + p["added_mass_matrix"]
    t1, t2, t6 = p["time_constants"]
    z4, z5 = p["damping_ratios"]
    t3 = t2
    dz = p["center_of_gravity"][2] - p["center_of_buoyancy"][2]
    w4 = np.sqrt(p["weight"] * dz / m[3, 3])
    w5 = np.sqrt(p["weight"] * dz / m[4, 4])
    d = np.diag([m[0, 0] / t1, m[1, 1] / t2, m[2, 2] / t3,
                 m[3, 3] * 2 * z4 * w4, m[4, 4] * 2 * z5 * w5, m[5, 5] / t6])
    d[0, 0] *= np.exp(-3 * np.linalg.norm(nu_r[:3]))  # remus100.m 127, 218
    return d


def _mss_otter_damping(p, nu_r):
    """otter.m 196-206 and 234-239; returns (derivatives, tau_damp)."""
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
    """crossFlowDrag.m 36-69 (current: 20 midpoints); ``strip_grid="endpoint"``
    gives the pre-2026-08-26 loop ``for xL = -L/2:dx:L/2`` of MSS.
    ``history_reynolds_on_length``: the cylinder Re of MSS before ``ac77394``."""
    length, beam, draft = p["length"], p["beam"], p["draft"]
    if p["drag_model"] == "cylinder":
        cd = _mss_cylinder_cd(length, beam, nu_r, length if history_reynolds_on_length else None)
    else:
        cd = _mss_hoerner(beam, draft)
    n_strips = int(_value(CFD, 37))
    dx = length / n_strips
    if p.get("strip_grid", "midpoint") == "midpoint":
        xs = -length / 2 + (np.arange(1, n_strips + 1) - 0.5) * dx   # line 56
    else:
        xs = -length / 2 + dx * np.arange(n_strips + 1)
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    u_h = np.abs(v_r + xs * r) * (v_r + xs * r)                      # line 61
    u_v = np.abs(w_r + xs * q) * (w_r + xs * q)                      # line 62
    k = -0.5 * p["water_density"] * draft * cd * dx                 # lines 63-66
    return np.array([0.0, k * u_h.sum(), k * u_v.sum(), 0.0, k * (xs * u_v).sum(), k * (xs * u_h).sum()])


def _mss_lift_drag(p, nu_r):
    """forceLiftDrag.m 26-40 + coeffLiftDrag.m 54-68 (sigma = 0), remus100.m 126-127."""
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
# Builders
# --------------------------------------------------------------------------
def _bind(block, values):
    """The block with its checked numbers and couplings bound in: a
    ``ca.Function`` of ``nu_r`` alone with every output of the block (the
    frozen path a plugin carries; couplings fixed as a vehicle class wires
    them)."""
    import casadi as ca

    nu_r = ca.SX.sym("nu_r", 6)
    out = block(nu_r=nu_r, **values)
    names = block.name_out()
    return ca.Function(f"{block.name()}_bound", [nu_r], [out[n] for n in names], ["nu_r"], names)


def _at_rest(function, *names):
    out = function(nu_r=np.zeros(6))
    return {n: np.array(out[n], dtype=float).ravel() for n in names}


def _split(p, selectors):
    p = dict(p)
    return {k: p.pop(k) for k in selectors if k in p}, p


def _smoothing(numbers):
    """The earlier selector value ``smooth_speed_epsilon`` (0 = MSS) as the
    selector ``smooth_speed`` and, when positive, the declared parameter."""
    numbers = dict(numbers)
    epsilon = numbers.pop("smooth_speed_epsilon", 0.0)
    if epsilon > 0.0:
        numbers["smooth_speed_epsilon"] = epsilon
    return {"smooth_speed": bool(epsilon > 0.0)}, numbers


def _reynolds_bound(form, numbers):
    """The earlier selector value ``ittc_reynolds_floor`` (``None`` = the MSS
    line) as the selector ``ittc_reynolds_bound`` and, with the floor, the
    declared parameter; absent = ``XuuITTC.m`` 34's floor."""
    numbers = dict(numbers)
    if form not in ("ittc", "exp_ittc"):
        numbers.pop("ittc_reynolds_floor", None)
        return {}, numbers
    floor = numbers.pop("ittc_reynolds_floor", _value(XITTC, 34))
    if floor is None:
        return {"ittc_reynolds_bound": "mss_offset"}, numbers
    numbers["ittc_reynolds_floor"] = floor
    return {"ittc_reynolds_bound": "floor"}, numbers


def _build_linear(form, p):
    block = _contract(LINEAR_DAMPING)
    selectors, numbers = _split(p, ("sway_damping_fade",))
    smoothing, numbers = _smoothing(numbers)
    function = _bind(block.linear_damping_casadi(form=form, **selectors, **smoothing),
                     block.check_linear_damping_values(numbers, form=form, **smoothing))
    name = "damping_derivatives" if form == "surface" else "damping_coefficients"
    return SimpleNamespace(**_at_rest(function, name)), function


def _build_submerged(p):
    return _build_linear("submerged", p)


def _build_surface(p):
    return _build_linear("surface", p)


def _build_cross_flow(p):
    block = _contract(CROSS_FLOW)
    selectors, numbers = _split(p, ("drag_model", "strip_grid"))
    values = block.check_cross_flow_drag_values(numbers, drag_model=selectors["drag_model"])
    return _bind(block.cross_flow_drag_casadi(**selectors), values)


def _build_lift_drag(p):
    from more_transformations.more_casadi_transformations import check_values

    block = _contract(LIFT_DRAG)
    selectors, numbers = _smoothing(p)
    return _bind(block.lift_drag_casadi(**selectors), check_values(block.lift_drag_parameters(**selectors), numbers))


def _attribute(module_name, name):
    block = _contract(module_name)
    if not hasattr(block, name):
        pytest.fail(f"block not ported yet: {module_name}.{name}")
    return getattr(block, name)


def _build_floating(p):
    return _build_linear("floating", p)


def _build_surge(p):
    block = _contract(SURGE_DAMPING)
    selectors, numbers = _split(p, ("quadratic_model", "surge_blend"))
    form = selectors.pop("quadratic_model")
    bound, numbers = _reynolds_bound(form, numbers)
    function = block.surge_damping_casadi(surge_form=form, **selectors, **bound)
    function = _bind(function, block.check_surge_damping_values(numbers, surge_form=form, **bound))
    names = [n for n in ("added_mass", "linear_coefficient", "quadratic_coefficient") if n in function.name_out()]
    rest = _at_rest(function, *names)
    return SimpleNamespace(**{n: float(v[0]) for n, v in rest.items()}), function


# --------------------------------------------------------------------------
# References of 2026-10-06: forceSurgeDamping.m, XuuITTC.m, Dmtrx.m surface branch
# (MATLAB direct calls at MSS ac77394, generate_u3_surge_floating_mss.m)
# --------------------------------------------------------------------------
def _named_columns(csv):
    header, values = _load_csv(csv)
    return {name: values[:, i] for i, name in enumerate(header)}


def _surge_parameters(set_id, model, mss=False):
    """Block parameters of one parameter set of the surge CSV (set 1 = osv.m).
    ``mss=True`` adds the flags that reproduce the MSS lines
    (``surge_blend="mss_tanh"``, ``ittc_reynolds_floor=None``).
    ``kinematic_viscosity`` (the ITTC forms) and ``surge_added_mass_factor``
    (both forms) are declared parameters now (rule 16): given here at the
    values the cited MSS lines fix (``forceSurgeDamping.m`` 69,
    ``addedMassSurge.m`` 33)."""
    ref = _named_columns(MSS_SURGE)
    row = int(np.flatnonzero(ref["set"] == set_id)[0])
    common = {"mass": ref["m"][row], "length": ref["L"][row], "water_density": ref["rho"][row],
              "time_constant": ref["T1"][row], "crossover_speed": _value(FSD, 57),
              "surge_added_mass_factor": _number_in(AMS, 33, r"A11 = ([0-9.]+) \*")}
    if mss:
        common["surge_blend"] = "mss_tanh"
    if model == "ittc":
        flags = {"ittc_reynolds_floor": None} if mss else {}
        return {"quadratic_model": "ittc", **common, "wetted_surface": ref["S"][row],
                "form_factor": _value(FSD, 70), "kinematic_viscosity": _value(FSD, 69), **flags}
    return {"quadratic_model": "max_thrust", **common, "max_speed": ref["u_max"][row],
            "max_thrust": ref["thrust_max"][row]}


def _surge_rows(set_id, model, probe=0):
    ref = _named_columns(MSS_SURGE)
    branch = 0 if model == "ittc" else 1
    rows = np.flatnonzero((ref["set"] == set_id) & (ref["branch"] == branch) & (ref["probe"] == probe))
    return {k: v[rows] for k, v in ref.items()}


def _mss_added_mass_surge(m, length, rho):
    """addedMassSurge.m 32-33 (Soding 1982): A11 = 2.7 rho nabla^(5/3) / L^2,
    the exponent and the L^2 pinned by the line text."""
    factor = _number_in(AMS, 33, r"A11 = ([0-9.]+) \*")
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
    # U_r and alpha of the reference are remus100.m 126-127 of its nu_r
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


@pytest.mark.parametrize("variable", ["MSS_DIR", "TESTS"])
def test_cited_lines_are_unchanged(variable):
    """With the variable set: every cited line and table of that checkout equals
    the snapshot (a moved or edited line fails here, as cylinderDrag.m did at
    ac77394). Without it: skips, naming the variable (nothing relative to one machine).
    TESTS = the template generators in this repository: always checked."""
    if variable != "TESTS":
        _env_dir(variable)
    snapshot = _snapshot()
    for (rel, number), text in CITED_LINES.items():
        if rel[1:].split("/", 1)[0].partition("@")[0] == variable:
            live = _live_line(rel, number)
            assert live.startswith(text), (rel, number, live)
            assert live == snapshot["lines"][rel][str(number)], (rel, number)
    if variable == "MSS_DIR":
        for rel, name in TABLES:
            live = _parse_table(_live_text(rel), name)
            assert np.array_equal(live, _table(rel, name)), (rel, name)


def test_transcriptions_reproduce_the_stored_matlab_reference():
    """The test-side MSS transcriptions, run on the legacy settings, land on the
    stored MATLAB numbers; so they are fit to serve as default-path references."""
    sph, cat = _spheroid_reference(), _catamaran_reference()
    p_d = _spheroid_damping(True)
    p_cf = _cross_flow_spheroid("endpoint")
    p_ld = _lift_drag_spheroid(_value(GEN_SPH, 492))
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
# G1 — block (legacy settings) vs the stored MATLAB reference
# --------------------------------------------------------------------------
def test_G1_block_spheroid_damping_with_sway_fade_flag():
    constants, function = _build_submerged(_spheroid_damping(True))
    ref = _spheroid_reference()
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], ref["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(out["tau"].reshape(-1), -ref["D"][k] @ ref["nu_r"][k]) <= G1_TOLERANCE, k


# ``test_G1_block_spheroid_cross_flow_legacy_settings`` was removed with the
# ``cross_flow_reynolds_length`` flag (owner 2026-10-06): the legacy spheroid file
# mixes the old Re-on-length line with the endpoint grid, so no flag-free block
# reproduces it. The endpoint grid keeps its G1 on the catamaran (Hoerner, no
# Reynolds number) below and its G2 on both vehicles.


def test_G1_block_spheroid_lift_drag_template_density():
    function = _build_lift_drag(_lift_drag_spheroid(_value(GEN_SPH, 492)))
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
# G2 — CasADi block signatures
# --------------------------------------------------------------------------


def test_G2_function_signatures():
    """The unfrozen blocks name ``nu_r``, their couplings and their parameters
    as inputs and name their outputs; the bound functions take ``nu_r`` only."""
    linear, surge = _contract(LINEAR_DAMPING), _contract(SURGE_DAMPING)
    for form in linear.LINEAR_DAMPING_FORMS:
        block = linear.linear_damping_casadi(form=form)
        names = [d.name for d in (*linear.linear_damping_couplings(form), *linear.linear_damping_parameters(form))]
        assert block.name_in() == ["nu_r", *names], form
        assert block.name_out() == list(linear.LINEAR_DAMPING_OUTPUTS[form]), form
    for form in surge.SURGE_FORMS:
        block = surge.surge_damping_casadi(surge_form=form)
        names = [d.name for d in (*surge.surge_damping_couplings(form), *surge.surge_damping_parameters(form))]
        assert block.name_in() == ["nu_r", *names], form
        assert block.name_out() == list(surge.SURGE_DAMPING_OUTPUTS[form]), form
    smooth = linear.linear_damping_casadi(form="submerged", smooth_speed=True)
    assert smooth.name_in() == [*linear.linear_damping_casadi(form="submerged").name_in(), "smooth_speed_epsilon"]
    lift = _contract(LIFT_DRAG)
    assert lift.lift_drag_casadi(smooth_speed=True).name_in() == [
        *lift.lift_drag_casadi().name_in(), "smooth_speed_epsilon"]
    for form in surge.SURGE_FORMS:
        floor = "ittc_reynolds_floor" in surge.surge_damping_casadi(surge_form=form).name_in()
        assert floor == (form in ("ittc", "exp_ittc")), form
    for form in ("ittc", "exp_ittc"):
        assert surge.surge_damping_casadi(surge_form=form).name_in()[-1] == "ittc_reynolds_floor", form
        assert "ittc_reynolds_floor" not in surge.surge_damping_casadi(
            surge_form=form, ittc_reynolds_bound="mss_offset").name_in(), form
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
        assert function.name_out()[:len(outs)] == outs
        sizes = {"D": (6, 6), "tau": (6, 1)}
        for i, name in enumerate(outs):
            assert function.size_out(i) == sizes[name], name


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
    """Default path against MATLAB at ac77394 (before 2026-10-06: the legacy file + flag)."""
    ref = _spheroid_reference(MSS_AC_SPHEROID)
    p = {k: v for k, v in _cross_flow_spheroid("midpoint").items() if k != "strip_grid"}
    assert _worst(_build_cross_flow(p), ref["nu_r"], ref["tau_crossflow"], "tau") <= G1_TOLERANCE
    broken = {**p, "draft": 1.01 * p["draft"]}
    assert _worst(_build_cross_flow(broken), ref["nu_r"], ref["tau_crossflow"], "tau") > G4_FACTOR * G1_TOLERANCE


def test_G4_lift_drag_area_plus_1_percent_is_detected():
    ref = _spheroid_reference()
    p = _lift_drag_spheroid(_value(GEN_SPH, 492))
    assert _worst(_build_lift_drag(p), ref["nu_r"], ref["tau_lift_drag"], "tau") <= G1_TOLERANCE
    broken = {**p, "planform_area": 1.01 * p["planform_area"]}
    assert _worst(_build_lift_drag(broken), ref["nu_r"], ref["tau_lift_drag"], "tau") > G4_FACTOR * G1_TOLERANCE


# --------------------------------------------------------------------------
# Headline (owner, 2026-10-05/06): default path vs current MSS
# --------------------------------------------------------------------------
def _states():
    return np.vstack([_spheroid_reference()["nu_r"], _random_nu_r()])


def test_MSS_submerged_damping_default_fades_surge_only():
    """Default = Dmtrx.m + remus100.m:218; the template's sway fade is the flag."""
    p = _spheroid_damping(False)
    _, function = _build_submerged(p)
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)
        assert _max_diff(out["D"], _mss_dmtrx(p, nu_r)) <= G1_TOLERANCE, ("Dmtrx.m 30-49, remus100.m 218", k)


def test_MSS_surface_damping_equals_otter():
    """otter.m 196-206, 234-239 with otter.m's own g = 9.81 (line 89)."""
    p = _surface_damping_catamaran(_value(OTTER, 90))
    constants, function = _build_surface(p)
    for k, nu_r in enumerate(np.vstack([_catamaran_reference()["nu_r"], _random_nu_r()])):
        deriv, tau = _mss_otter_damping(p, nu_r)
        out = _evaluate(function, nu_r)
        assert _max_diff(out["tau"].reshape(-1), tau) <= G1_TOLERANCE, ("otter.m 234-239", k)
    assert _max_diff(constants.damping_derivatives, deriv) <= G1_TOLERANCE


@pytest.mark.parametrize("model", ["cylinder", "hoerner"])
def test_MSS_cross_flow_default_equals_current_mss(model):
    """Default (no flags) = crossFlowDrag.m with cylinderDrag.m at ac77394 or
    later (Re on the diameter, lines 78-80) on 1050 states."""
    p = _cross_flow_spheroid("midpoint") if model == "cylinder" else _cross_flow_catamaran("midpoint")
    function = _build_cross_flow({k: v for k, v in p.items() if k != "strip_grid"})
    for k, nu_r in enumerate(_states()):
        out = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert _max_diff(out, _mss_cross_flow(p, nu_r)) <= G1_TOLERANCE, ("crossFlowDrag.m 36-69", k)


def test_cross_flow_length_flag_is_gone():
    """Owner 2026-10-06: with MSS on the diameter since ac77394 the
    ``cross_flow_reynolds_length`` flag is dropped; ``strip_grid`` stays with
    default "midpoint"."""
    block = _contract(CROSS_FLOW)
    signature = inspect.signature(block.cross_flow_drag_casadi)
    assert "cross_flow_reynolds_length" not in signature.parameters, list(signature.parameters)
    assert "cross_flow_reynolds_length" not in [d.name for d in block.cross_flow_drag_parameters()]
    assert signature.parameters["strip_grid"].default == "midpoint"


def test_endpoint_grid_flag_reproduces_the_template():
    """strip_grid="endpoint" equals the pre-2026-08-26 grid of the stored
    catamaran reference (Hoerner: no Reynolds number); the default grid does
    not (the spheroid file would need the dropped length flag)."""
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
# G1-MSS: default paths vs MATLAB running current MSS (2026-10-05)
# --------------------------------------------------------------------------
def test_transcriptions_reproduce_the_current_mss_matlab_reference():
    """The second check: the test-side transcriptions, on current-MSS settings,
    land on the 2026-10-05 files and, for the spheroid cross-flow, on the
    ac77394 file (2026-10-06). The 2026-10-05 spheroid
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


def test_transcriptions_reproduce_the_cc07579_catamaran_reference():
    """The otter.m transcriptions land on the cc07579 catamaran file (the KB
    fix of cf349d4 moved G44, G55, Kp, Mq; the damping law did not change)."""
    cat = _catamaran_reference(MSS_CATAMARAN_CC07579)
    p_cat = _surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN_CC07579)
    p_cat_cf = _cross_flow_catamaran("midpoint")
    for k in range(50):
        deriv, tau = _mss_otter_damping(p_cat, cat["nu_r"][k])
        assert _max_diff(np.diag(deriv), cat["D"][k]) <= G1_TOLERANCE, k
        assert _max_diff(tau, cat["tau_damp"][k]) <= G1_TOLERANCE, k
        assert _max_diff(_mss_cross_flow(p_cat_cf, cat["nu_r"][k]), cat["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_G1_MSS_submerged_damping_default_equals_matlab():
    """Default (surge fade only, remus100.m 218) with M_RB, M_A of the 2026-10-05 file
    (imlay61.m rho = 1026) vs MATLAB's D and -D nu_r."""
    ref = _spheroid_reference(MSS_SPHEROID)
    p = _spheroid_damping(False, MSS_SPHEROID)
    p = {k: v for k, v in p.items() if k != "sway_damping_fade"}   # the default
    _, function = _build_submerged(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), -ref["D"][k] @ ref["nu_r"][k]) <= G1_TOLERANCE, ("tau", k)


# ``test_G1_MSS_cross_flow_length_flag_equals_matlab`` was removed with the flag
# (owner 2026-10-06); its MATLAB file (MSS 99bf0b3) stays, read by the history
# checks above and below.


def test_G1_MSS_cross_flow_default_equals_matlab_ac77394():
    """Reynolds number on the diameter: the default equals MATLAB running
    MSS ac77394 (the 2026-10-05 generator re-run unchanged, SOURCE.md)."""
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


def test_cross_flow_reynolds_on_diameter_effect_vs_pre_fix_mss():
    """History of the Reynolds-number fix: the default departs from MATLAB
    at 99bf0b3 (Re on the length) on the reference cases, and the stated
    effect holds: REMUS Y_v|v| = -125 (diameter, the default) vs -66
    N/(m/s)^2 (length, pre-fix MSS, now a transcription only) at 0.3 m/s
    sideways (whole numbers recorded when the fix was proposed to MSS,
    2026-10-05)."""
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


@pytest.mark.parametrize("csv", [MSS_CATAMARAN, MSS_CATAMARAN_CC07579])
def test_G1_MSS_catamaran_surface_damping_equals_matlab(csv):
    """otter.m g = 9.81 (line 90); M and G of the 2026-10-05 file (inertia about the
    combined CG) and of the cc07579 file (the KB fix changes G44, G55 and with
    them Kp, Mq: the block takes G as input, so it must follow)."""
    ref = _catamaran_reference(csv)
    _, function = _build_surface(_surface_damping_catamaran(_value(OTTER, 90), csv))
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])
        assert _max_diff(out["D"], -ref["D"][k]) <= G1_TOLERANCE, ("D", k)
        assert _max_diff(out["tau"].reshape(-1), ref["tau_damp"][k]) <= G1_TOLERANCE, ("tau", k)


@pytest.mark.parametrize("csv", [MSS_CATAMARAN, MSS_CATAMARAN_CC07579])
def test_G1_MSS_catamaran_cross_flow_default_equals_matlab(csv):
    """Hoerner on the default grid (20 midpoints): no deviation from MSS (at
    cc07579 too)."""
    ref = _catamaran_reference(csv)
    p = {k: v for k, v in _cross_flow_catamaran("midpoint").items() if k != "strip_grid"}
    function = _build_cross_flow(p)
    for k in range(50):
        out = _evaluate(function, ref["nu_r"][k])["tau"].reshape(-1)
        assert _max_diff(out, ref["tau_crossflow"][k]) <= G1_TOLERANCE, k


def test_hoerner_below_its_table_raises_naming_the_ratio():
    """Deviation from MSS (owner, 2026-10-05): B/(2T) below Hoerner.m's first data point (MSS:
    interp1 -> NaN, lines 47-48) raises ValueError naming the ratio; the first
    data point itself is accepted."""
    block = _contract(CROSS_FLOW)
    first = _table(HOERNER, "CD_DATA")[0, 0]
    draft = 1.0
    beam = first * draft          # ratio = first / 2, below the table
    ratio = beam / (2 * draft)
    with pytest.raises(ValueError) as error:
        block.check_cross_flow_drag_values(
            {"length": 2.0, "beam": beam, "draft": draft, "water_density": _value(CFD, 36)},
            drag_model="hoerner")
    assert f"{ratio}" in str(error.value), str(error.value)
    block.check_cross_flow_drag_values(
        {"length": 2.0, "beam": 2 * first * draft, "draft": draft, "water_density": _value(CFD, 36)},
        drag_model="hoerner")


# --------------------------------------------------------------------------
# Floating linear damping: Dmtrx.m surface-craft branch
# --------------------------------------------------------------------------
def test_floating_transcription_equals_matlab():
    """Second check, no block: the Dmtrx.m 51-63 transcription lands on MATLAB."""
    ref = _floating_reference()
    for k in range(len(ref["D"])):
        assert _max_diff(_mss_dmtrx_floating(_floating_parameters(k, ref)), ref["D"][k]) <= G1_TOLERANCE, k


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


def test_G4_floating_roll_ratio_times_zero_is_detected():
    ref = _floating_reference()
    p = _floating_parameters(0, ref)
    assert _max_diff(_build_floating(p)[0].damping_coefficients, np.diag(ref["D"][0])) <= G1_TOLERANCE
    broken = {**p, "damping_ratios": [p["damping_ratios"][0], 0.0, p["damping_ratios"][2]]}
    assert _max_diff(_build_floating(broken)[0].damping_coefficients,
                     np.diag(ref["D"][0])) > G4_FACTOR * G1_TOLERANCE


# --------------------------------------------------------------------------
# Surge damping: forceSurgeDamping.m at ac77394 (source force_surge_damping)
# --------------------------------------------------------------------------
SURGE_MODELS = ("ittc", "max_thrust")
# G2_SCALE_NOTE: G2's 1e-10 is absolute; one ulp of a 7.5e5 N force is 1.2e-10,
# so two correct implementations can miss it at ship scale (probe of 2026-10-06:
# 1.2e-10 and 2.3e-10). The G2 tests stay below 1e5 N; ship scale (osv.m) is held
# by G1 at 1e-9. A relative G2 for large craft is an open question.


def test_surge_reference_set_1_is_osv():
    """Parameter set 1 of the surge CSV is osv.m's vessel (lines 61-79, 128),
    the only MSS vehicle that called forceSurgeDamping (osv.m 180-181 at
    ac77394, the revision of the surge reference)."""
    ref = _named_columns(MSS_SURGE)
    row = int(np.flatnonzero(ref["set"] == 1)[0])
    length, beam, draft = _value(OSV, 61), _value(OSV, 62), _value(OSV, 63)
    rho, cb = _value(OSV, 64), _value(OSV, 65)
    k_max = _value(OSV, 69)
    assert ref["L"][row] == length and ref["rho"][row] == rho
    assert ref["m"][row] == rho * cb * length * beam * draft                  # lines 78-79
    assert ref["S"][row] == length * beam + 2 * draft * beam                  # line 66
    assert ref["T1"][row] == _value(OSV, 128) and ref["u_max"][row] == _value(OSV, 76)
    assert ref["thrust_max"][row] == k_max[2] + k_max[3]                      # line 75


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
    or any MSS vehicle): the reference for the bounded-Reynolds default
    (owner's decision of 2026-10-06)."""
    ref = _named_columns(MSS_XUU_ITTC)
    for k in range(len(ref["u_r"])):
        got = _mss_xuu_ittc(ref["u_r"][k], ref["rho"][k], ref["L"][k], ref["B"][k], ref["T"][k], ref["C_B"][k])
        assert abs(got - ref["Xuu"][k]) <= G1_TOLERANCE, k


@pytest.mark.parametrize("model", SURGE_MODELS)
def test_G1_MSS_surge_damping_equals_matlab(model):
    """Block vs MATLAB forceSurgeDamping.m on five sets (osv.m + four seeded),
    40 speeds each in -8..8 m/s; constants A11, Xu (and Xuu for max_thrust)."""
    for set_id in range(1, 6):
        constants, function = _build_surge(_surge_parameters(set_id, model, mss=True))
        rows = _surge_rows(set_id, model)
        assert abs(constants.added_mass - rows["A11"][0]) <= G1_TOLERANCE, set_id
        assert abs(constants.linear_coefficient - rows["Xu"][0]) <= G1_TOLERANCE, set_id
        if model == "max_thrust":
            assert abs(constants.quadratic_coefficient - rows["Xuu"][0]) <= G1_TOLERANCE, set_id
        for k, u in enumerate(rows["u_r"]):
            nu_r = np.array([u, 0.3, -0.2, 0.1, 0.05, -0.1])   # other DOFs must not enter
            tau = _evaluate(function, nu_r)["tau"].reshape(-1)
            assert _max_diff(tau, [rows["X"][k], 0, 0, 0, 0, 0]) <= G1_TOLERANCE, (set_id, k)


def test_G4_surge_time_constant_plus_1_percent_is_detected():
    rows = _surge_rows(2, "max_thrust")
    p = _surge_parameters(2, "max_thrust", mss=True)
    worst = lambda q: max(abs(_evaluate(_build_surge(q)[1], np.array([u, 0, 0, 0, 0, 0.0]))["tau"][0, 0] - x)
                          for u, x in zip(rows["u_r"], rows["X"]))
    assert worst(p) <= G1_TOLERANCE
    assert worst({**p, "time_constant": 1.01 * p["time_constant"]}) > G4_FACTOR * G1_TOLERANCE


# Physical-sign tests (MSS is the reference, not the truth)
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
    Owner's decision of 2026-10-06: our default is ``surge_blend="symmetric"``,
    tested below; MSS stays behind ``surge_blend="mss_tanh"``."""
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
    Re at 1e5 (line 34). Owner's decision of 2026-10-06: our default is
    ``ittc_reynolds_floor=1e5``, tested below; MSS stays behind
    ``ittc_reynolds_floor=None``."""
    ittc = _surge_rows(1, "ittc", probe=1)
    thrust = _surge_rows(1, "max_thrust", probe=1)
    assert np.max(np.abs(ittc["X"])) > 1e9
    assert np.max(np.abs(thrust["X"])) < 1.0


# Owner defaults (2026-10-06): symmetric blend, Rn floor
def _surge_x(function, u):
    return _evaluate(function, np.array([u, 0, 0, 0, 0, 0.0]))["tau"][0, 0]


@pytest.mark.parametrize("model", SURGE_MODELS)
def test_surge_defaults_equal_matlab_ahead_above_the_floor(model):
    """Where MSS is dissipative and above the Rn floor (u_r >= 0 and, for
    ITTC, L |u_r| / nu >= 1e5) the default equals MATLAB at G1 (1e-9)."""
    for set_id in range(1, 6):
        _, function = _build_surge(_surge_parameters(set_id, model))
        rows = _surge_rows(set_id, model)
        keep = rows["u_r"] >= 0.0
        if model == "ittc":
            keep &= rows["L"] * np.abs(rows["u_r"]) / _value(FSD, 69) >= _value(XITTC, 34)
        assert keep.sum() >= 10, (set_id, model)
        for u, x in zip(rows["u_r"][keep], rows["X"][keep]):
            assert abs(_surge_x(function, u) - x) <= G1_TOLERANCE, (set_id, model, u)


@pytest.mark.parametrize("model", SURGE_MODELS)
def test_surge_symmetric_blend_is_dissipative_on_every_row_and_probe(model):
    """Fossen (2011) Property 6.3 (p. 123): X u_r <= 0 on all 440 speeds of
    the surge CSV (ordinary rows and the singular probes), five sets."""
    ref = _named_columns(MSS_SURGE)
    for set_id in range(1, 6):
        _, function = _build_surge(_surge_parameters(set_id, model))
        for u in ref["u_r"][ref["set"] == set_id]:
            assert _surge_x(function, u) * u <= 0.0, (set_id, model, u)


def test_surge_ittc_floor_is_finite_on_the_singular_probes():
    """The 40 probe rows (Rn = 100 .. 100.1, MSS's pole, both branches, five
    sets): the default is finite on every one (MSS: |X| > 1e9 N on osv.m),
    and below 1 N on set 1 (osv.m) for the ITTC branch."""
    for set_id in range(1, 6):
        for model in SURGE_MODELS:
            _, function = _build_surge(_surge_parameters(set_id, model))
            rows = _surge_rows(set_id, model, probe=1)
            assert len(rows["u_r"]) == 4, (set_id, model)
            x = np.array([_surge_x(function, u) for u in rows["u_r"]])
            assert np.all(np.isfinite(x)), (set_id, model)
            if set_id == 1 and model == "ittc":
                assert np.max(np.abs(x)) < 1.0, x


@pytest.mark.parametrize("floor", [100.0, 50.0, 0.0, -1.0, np.inf, np.nan])
def test_surge_ittc_floor_at_or_below_the_pole_raises(floor):
    p = _surge_parameters(2, "ittc")
    with pytest.raises(ValueError, match="ittc_reynolds_floor"):
        _build_surge({**p, "ittc_reynolds_floor": floor})


def test_surge_reynolds_bound_selector_is_checked():
    """``ittc_reynolds_bound`` is a selector of the ITTC forms only; an
    unknown value, or the MSS line on ``"max_thrust"``, raises naming it."""
    surge = _contract(SURGE_DAMPING)
    with pytest.raises(ValueError, match="ittc_reynolds_bound"):
        surge.surge_damping_casadi(surge_form="ittc", ittc_reynolds_bound="none")
    with pytest.raises(ValueError, match="ittc_reynolds_bound"):
        surge.surge_damping_casadi(surge_form="max_thrust", ittc_reynolds_bound="mss_offset")


# --------------------------------------------------------------------------
# kinematic_viscosity and surge_added_mass_factor are declared parameters
# (rule 16, E-65): a water property and an empirical coefficient, not fixed
# module constants.
# --------------------------------------------------------------------------
def test_viscosity_and_added_mass_factor_are_declared_with_unit_and_range():
    """``kinematic_viscosity`` (m^2/s, > 0, the ITTC forms) and
    ``surge_added_mass_factor`` (1, > 0, ``"ittc"`` and ``"max_thrust"``) are
    ``Parameter`` declarations of ``surge_damping_parameters``."""
    surge = _contract(SURGE_DAMPING)
    for form in ("ittc", "exp_ittc"):
        declared = {d.name: d for d in surge.surge_damping_parameters(form)}
        viscosity = declared["kinematic_viscosity"]
        assert (viscosity.shape, viscosity.unit) == ((1, 1), "m^2/s")
        assert (viscosity.minimum, viscosity.minimum_exclusive, viscosity.maximum) == (0.0, True, None)
    for form in ("ittc", "max_thrust"):
        declared = {d.name: d for d in surge.surge_damping_parameters(form)}
        factor = declared["surge_added_mass_factor"]
        assert (factor.shape, factor.unit) == ((1, 1), "1")
        assert (factor.minimum, factor.minimum_exclusive, factor.maximum) == (0.0, True, None)
    assert "kinematic_viscosity" not in [d.name for d in surge.surge_damping_parameters("max_thrust")]
    assert "surge_added_mass_factor" not in [d.name for d in surge.surge_damping_parameters("exp_ittc")]


def test_committed_equality_of_surge_at_the_declared_defaults():
    """The block called with ``kinematic_viscosity`` and
    ``surge_added_mass_factor`` at the values the cited MSS lines fix
    (``forceSurgeDamping.m`` 69, ``addedMassSurge.m`` 33) equals the MATLAB
    rows at G1 (1e-9): turning the fixed constants into named inputs changes
    nothing at those values (five sets, both forms)."""
    for set_id in range(1, 6):
        for model in SURGE_MODELS:
            _, function = _build_surge(_surge_parameters(set_id, model))
            rows = _surge_rows(set_id, model)
            keep = rows["u_r"] >= 0.0
            if model == "ittc":
                keep &= rows["L"] * np.abs(rows["u_r"]) / _value(FSD, 69) >= _value(XITTC, 34)
            for u, x in zip(rows["u_r"][keep], rows["X"][keep]):
                assert abs(_surge_x(function, u) - x) <= G1_TOLERANCE, (set_id, model, u)


def test_G4_kinematic_viscosity_moves_the_ittc_line():
    """``nu`` is a live input of the ITTC forms: at the declared default
    (``forceSurgeDamping.m`` 69 / ``XuuITTC.m`` 32) ``X`` still equals MATLAB
    at G1 (1e-9, the committed-equality check at the default); +1 % moves
    ``Rn`` and so ``X`` by more than G4's factor times G1, on ``"ittc"`` (5.0
    m/s, well above the Reynolds floor) and on the OSV ``"exp_ittc"`` rows."""
    p = _surge_parameters(2, "ittc")
    u = 5.0
    base = _build_surge(p)[1]
    moved = _build_surge({**p, "kinematic_viscosity": 1.01 * p["kinematic_viscosity"]})[1]
    assert abs(_surge_x(moved, u) - _surge_x(base, u)) > G4_FACTOR * G1_TOLERANCE

    ref = _osv_surge_reference()
    rows = np.flatnonzero((ref["set"] == 2) & (ref["probe"] == 0))
    worst = lambda scale: max(
        abs(_evaluate(_build_surge({**_exp_ittc_parameters(ref, k),
                                    "kinematic_viscosity": scale * _exp_ittc_parameters(ref, k)["kinematic_viscosity"]})[1],
                      np.array([ref["u_r"][k], 0, 0, 0, 0, 0.0]))["tau"][0, 0] - ref["X"][k])
        for k in rows)
    assert worst(1.0) <= G1_TOLERANCE
    assert worst(1.01) > G4_FACTOR * G1_TOLERANCE


@pytest.mark.parametrize("value", [0.0, -1e-6, np.inf, np.nan])
def test_kinematic_viscosity_at_or_below_zero_is_refused(value):
    """``check_values`` refuses ``nu <= 0`` and a non-finite ``nu``, naming
    the parameter and its unit, on both ITTC forms."""
    p = _surge_parameters(2, "ittc")
    with pytest.raises(ValueError, match=r"'kinematic_viscosity' \[m\^2/s\]"):
        _build_surge({**p, "kinematic_viscosity": value})
    ref = _osv_surge_reference()
    exp_p = _exp_ittc_parameters(ref, 0)
    with pytest.raises(ValueError, match=r"'kinematic_viscosity' \[m\^2/s\]"):
        _build_surge({**exp_p, "kinematic_viscosity": value})


@pytest.mark.parametrize("value", [0.0, -2.7, np.inf, np.nan])
def test_surge_added_mass_factor_at_or_below_zero_is_refused(value):
    """``check_values`` refuses ``c <= 0`` and a non-finite ``c``, naming the
    parameter and its unit, on both ``"ittc"`` and ``"max_thrust"``."""
    for model in SURGE_MODELS:
        p = _surge_parameters(2, model)
        with pytest.raises(ValueError, match=r"'surge_added_mass_factor' \[1\]"):
            _build_surge({**p, "surge_added_mass_factor": value})


# --------------------------------------------------------------------------
# Regularisations are declared parameters (named inputs); their form is a selector
# --------------------------------------------------------------------------
def test_regularisation_parameters_are_declared_with_unit_and_range():
    """``smooth_speed_epsilon`` (m/s, > 0) and ``ittc_reynolds_floor`` (1,
    > 100, the pole of the ITTC line, Fossen 2011, eq. 6.83, p. 125) are
    ``Parameter`` declarations, present only with the selector that uses them."""
    linear, lift, surge = (_contract(m) for m in (LINEAR_DAMPING, LIFT_DRAG, SURGE_DAMPING))
    for declared in (linear.linear_damping_parameters("submerged", smooth_speed=True),
                     lift.lift_drag_parameters(smooth_speed=True)):
        epsilon = declared[-1]
        assert (epsilon.name, epsilon.shape, epsilon.unit) == ("smooth_speed_epsilon", (1, 1), "m/s")
        assert (epsilon.minimum, epsilon.minimum_exclusive, epsilon.maximum) == (0.0, True, None)
    for declared in (linear.linear_damping_parameters("submerged"), lift.lift_drag_parameters()):
        assert "smooth_speed_epsilon" not in [d.name for d in declared]
    for form in ("ittc", "exp_ittc"):
        floor = surge.surge_damping_parameters(form)[-1]
        assert (floor.name, floor.shape, floor.unit) == ("ittc_reynolds_floor", (1, 1), "1")
        assert (floor.minimum, floor.minimum_exclusive, floor.maximum) == (100.0, True, None)
        assert "ittc_reynolds_floor" not in [
            d.name for d in surge.surge_damping_parameters(form, ittc_reynolds_bound="mss_offset")]
    assert "ittc_reynolds_floor" not in [d.name for d in surge.surge_damping_parameters("max_thrust")]
    for form in ("floating", "surface"):
        with pytest.raises(ValueError, match="smooth_speed"):
            linear.linear_damping_casadi(form=form, smooth_speed=True)
    for build in (lambda: linear.linear_damping_casadi(form="submerged", smooth_speed=1e-3),
                  lambda: lift.lift_drag_casadi(smooth_speed=1e-3)):
        with pytest.raises(ValueError, match="smooth_speed must be True or False"):
            build()


@pytest.mark.parametrize("epsilon", [0.0, -1e-3, np.inf, np.nan])
def test_smooth_speed_epsilon_at_or_below_zero_is_refused(epsilon):
    """``check_values`` (and the linear-damping check) refuse ``eps <= 0`` and
    a non-finite ``eps``, naming the parameter and its unit."""
    from more_transformations.more_casadi_transformations import check_values

    linear, lift = _contract(LINEAR_DAMPING), _contract(LIFT_DRAG)
    submerged = {k: v for k, v in _spheroid_damping(False, MSS_SPHEROID).items() if k != "sway_damping_fade"}
    with pytest.raises(ValueError, match=r"'smooth_speed_epsilon' \[m/s\]"):
        linear.check_linear_damping_values({**submerged, "smooth_speed_epsilon": epsilon},
                                           form="submerged", smooth_speed=True)
    with pytest.raises(ValueError, match=r"'smooth_speed_epsilon' \[m/s\]"):
        check_values(lift.lift_drag_parameters(smooth_speed=True),
                     {**_lift_drag_spheroid(_value(FLD, 26)), "smooth_speed_epsilon": epsilon})


REGULARISATION_PROBE = 0.1   # m/s; a probe value of this test, large enough to move the force near rest


def test_G4_smooth_speed_epsilon_moves_the_output():
    """The regularisation is a live input: near rest, eps +1 % moves the
    submerged damping force and the lift/drag force by more than G4's
    factor times G1; the same blocks called with eps and with the bound
    numbers are equal (the input is the one used)."""
    eps = REGULARISATION_PROBE
    nu_r = np.array([eps, 0.0, eps, 0.0, 0.0, 0.0])
    submerged = {k: v for k, v in _spheroid_damping(False, MSS_SPHEROID).items()}
    lift = _lift_drag_spheroid(_value(FLD, 26))
    for label, build, p in (("submerged", lambda q: _build_submerged(q)[1], submerged),
                            ("lift_drag", _build_lift_drag, lift)):
        base = _evaluate(build({**p, "smooth_speed_epsilon": eps}), nu_r)["tau"].reshape(-1)
        moved = _evaluate(build({**p, "smooth_speed_epsilon": 1.01 * eps}), nu_r)["tau"].reshape(-1)
        assert _max_diff(moved, base) > G4_FACTOR * G1_TOLERANCE, (label, base, moved)


def test_G4_ittc_reynolds_floor_moves_the_output_below_the_floor():
    """Below the floor (``L |u_r| / nu < Rn_min``) the ITTC coefficient is
    evaluated at ``Rn_min``: the floor +1 % moves ``X`` by more than G4's
    factor times G1; above the floor it does not enter."""
    p = _surge_parameters(2, "ittc")
    floor = _value(XITTC, 34)
    low = 0.5 * floor * _value(FSD, 69) / p["length"]     # Rn = floor / 2
    high = 10.0 * floor * _value(FSD, 69) / p["length"]   # Rn = 10 floor
    base = _build_surge({**p, "ittc_reynolds_floor": floor})[1]
    moved = _build_surge({**p, "ittc_reynolds_floor": 1.01 * floor})[1]
    assert abs(_surge_x(moved, low) - _surge_x(base, low)) > G4_FACTOR * G1_TOLERANCE
    assert _surge_x(moved, high) == _surge_x(base, high)


# Inputs outside their domain raise ValueError naming the input (MSS returns NaN/inf/complex)
def _raises_naming(build, name):
    with pytest.raises(ValueError, match=name):
        build()


def test_inputs_outside_their_domain_raise_naming_the_input():
    """Floating G33 < 0; surge mass <= 0, max_speed <= 0, an unknown blend.
    Negative damping ratios (added 2026-10-07): our extension of the input
    checks, MSS does not reject them; zero is allowed (an undamped mode)."""
    floating = _floating_parameters(0)
    g = np.array(floating["restoring_matrix"], float)
    g[2, 2] = -abs(g[2, 2]) - 1.0
    _raises_naming(lambda: _build_floating({**floating, "restoring_matrix": g}), "G33")
    for model in SURGE_MODELS:
        p = _surge_parameters(2, model)
        _raises_naming(lambda: _build_surge({**p, "mass": 0.0}), "mass")
        _raises_naming(lambda: _build_surge({**p, "mass": -1.0}), "mass")
        _raises_naming(lambda: _build_surge({**p, "surge_blend": "x"}), "surge_blend")
    thrust = _surge_parameters(2, "max_thrust")
    _raises_naming(lambda: _build_surge({**thrust, "max_speed": 0.0}), "max_speed")
    _raises_naming(lambda: _build_surge({**thrust, "max_speed": -1.0}), "max_speed")

    submerged = {k: v for k, v in _spheroid_damping(False, MSS_SPHEROID).items()}
    surface = _surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN)
    for build, p, n in ((_build_submerged, submerged, 2), (_build_floating, floating, 3),
                        (_build_surface, surface, 3)):
        for k in range(n):
            ratios = list(p["damping_ratios"])
            ratios[k] = -0.1
            _raises_naming(lambda: build({**p, "damping_ratios": ratios}), "damping_ratios")
            ratios[k] = 0.0
            constants, _ = build({**p, "damping_ratios": ratios})   # zero is allowed


def test_physical_cross_flow_equals_first_principles_strip_integral():
    """The strip at x (from the CO, z down) moves with v + x r sideways and
    w - x q vertically (omega x r for r = [x 0 0]); its moments are
    N = +x Y_i and M = -x Z_i. Fossen (2011) eqs. 6.91-6.92 (p. 127) state the
    sway/yaw half; heave/pitch is derived here the same way. crossFlowDrag.m
    62 and 65 use w + x q and +x Z: two sign flips that cancel on a grid
    symmetric about the CO (found 2026-10-05), so the default block must equal
    the first-principles form; a pure pitch or yaw rate must be opposed."""
    for p in (_cross_flow_spheroid("midpoint"), _cross_flow_catamaran("midpoint")):
        function = _build_cross_flow({k: v for k, v in p.items() if k != "strip_grid"})
        n_strips = int(_value(CFD, 37))
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
# Jacobians at and near rest with the MSS lines (smooth_speed_epsilon = 0)
# --------------------------------------------------------------------------
JACOBIAN_STATES = {
    "rest": [0, 0, 0, 0, 0, 0],
    "pure_sway": [0, 0.3, 0, 0, 0, 0],             # sideways current, no way on
    "pure_rotation": [0, 0, 0, 0.1, 0.1, 0.1],     # turning on the spot
}
# (block, output, state) -> non-finite Jacobian entries (row of the output,
# column of nu_r), probed 2026-10-06. Submerged damping: U_r = |nu_r[0:3]|
# (remus100.m 127) inside exp(-3 U_r) has d/dnu = nu/U_r = 0/0 when u=v=w=0.
# Lift/drag: atan2(w_r, u_r) (remus100.m 126) has 0/0 derivatives when
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
def test_nonfinite_jacobian_entries_of_the_mss_lines_are_documented(state):
    """Documents the MSS lines' behaviour (the reason for
    ``smooth_speed_epsilon``): exactly the listed entries are non-finite; cross-flow and surface
    damping are finite; 1e-12 off the singular set every Jacobian is finite."""
    for name, function in _jacobian_blocks().items():
        expected = JACOBIAN_NONFINITE_TODAY.get((name, "tau", state), set())
        assert _nonfinite_jacobian(function, "tau", JACOBIAN_STATES[state]) == expected, (name, state)
        assert not _nonfinite_jacobian(function, "tau", np.full(6, 1e-12)), name


SMOOTH_EPSILONS = (1e-6, 1e-3)   # m/s; probe values of 2026-10-06


def _smooth_blocks(epsilon):
    submerged = {k: v for k, v in _spheroid_damping(False, MSS_SPHEROID).items()}
    return {
        "submerged": _build_submerged({**submerged, "smooth_speed_epsilon": epsilon})[1],
        "lift_drag": _build_lift_drag({**_lift_drag_spheroid(_value(FLD, 26)),
                                       "smooth_speed_epsilon": epsilon}),
    }


# States on both sides of the negative-surge cut of atan2 (u_r < 0, w_r -> 0+
# and 0-), inside and outside the band |u_r w_r| < K eps^2 of the lift_drag
# docstring, with sway and rotation on and off.
ASTERN_CUT_STATES = [
    [-1, 0, 0, 0, 0, 0],                          # pure astern probe
    *[[u, v, side * w, p, q, r]
      for u in (-0.01, -1.0, -5.0)
      for w in (1e-12, 1e-9, 1e-6, 1e-3, 0.1, 1.0)
      for side in (1, -1)
      for v, p, q, r in ((0, 0, 0, 0), (0.3, 0.1, -0.05, 0.2))],
]


@pytest.mark.parametrize("epsilon", SMOOTH_EPSILONS)
def test_smooth_epsilon_gives_finite_jacobians(epsilon):
    """With ``smooth_speed_epsilon > 0`` the submerged
    damping and lift/drag Jacobians are finite on every JACOBIAN_STATES entry
    and at the astern-cut states; eps = 0 reproduces the pinned set
    (``test_nonfinite_jacobian_entries_of_the_mss_lines_are_documented`` above)."""
    for name, function in _smooth_blocks(epsilon).items():
        for state in [*JACOBIAN_STATES.values(), *ASTERN_CUT_STATES]:
            assert not _nonfinite_jacobian(function, "tau", state), (name, epsilon, state)


@pytest.mark.parametrize("epsilon", SMOOTH_EPSILONS)
def test_smooth_lift_drag_opposes_the_motion(epsilon):
    """Dissipation (Fossen 2011, Property 6.3, p. 123): with eps > 0 the
    lift/drag force does no positive work, tau . nu_r <= 0, at the pure
    astern probe, on both sides of the negative-surge cut, near rest and on
    the 1000 seeded states of G2; at the probe X > 0 (it opposes the astern motion).
    The default eps = 0 (MSS) is dissipative off the cut too (lift does no
    work there); the first smooth form pushed the craft astern (X < 0)."""
    function = _smooth_blocks(epsilon)["lift_drag"]
    near_rest = [[s * a, 0.2, t * b, 0, 0, 0] for a in (0, 1e-9, epsilon, 10 * epsilon)
                 for b in (0, 1e-9, epsilon) for s in (1, -1) for t in (1, -1)]
    states = [*ASTERN_CUT_STATES, *near_rest, *_random_nu_r()]
    for state in states:
        nu_r = np.asarray(state, float)
        tau = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert np.all(np.isfinite(tau)), state
        assert tau @ nu_r <= 0.0, (epsilon, state, tau)
    probe = _evaluate(function, np.array([-1.0, 0, 0, 0, 0, 0]))["tau"].reshape(-1)
    assert probe[0] > 0.0, probe
    mss = _build_lift_drag(_lift_drag_spheroid(_value(FLD, 26)))
    for state in [[-1, 0, 0.3, 0, 0, 0], [-1, 0, -0.3, 0, 0, 0], [2, 0.1, 0.5, 0, 0, 0]]:
        nu_r = np.asarray(state, float)
        assert _evaluate(mss, nu_r)["tau"].reshape(-1) @ nu_r <= 0.0, state


@pytest.mark.parametrize("epsilon", SMOOTH_EPSILONS)
def test_smooth_lift_drag_equals_mss_outside_the_band(epsilon):
    """The band stated in the lift_drag docstring, K = 1000: outside
    {u^2 + w^2 < (K eps)^2} U {u < 0, |u w| < K eps^2} the smooth force is
    within 1e-3 of the MSS force's magnitude (forceLiftDrag.m transcription),
    on 1000 seeded log-spaced states of both signs outside the band (|u|, |w|
    up to 30 m/s)."""
    k_band = 1000.0
    p = _lift_drag_spheroid(_value(FLD, 26))
    function = _build_lift_drag({**p, "smooth_speed_epsilon": epsilon})
    rng = np.random.default_rng(SEED)
    checked = 0
    while checked < 1000:
        u = rng.choice([-1.0, 1.0]) * 10 ** rng.uniform(-9, 1.5)
        w = rng.choice([-1.0, 1.0]) * 10 ** rng.uniform(-12, 1.5)
        if u * u + w * w < (k_band * epsilon) ** 2 or (u < 0 and abs(u * w) < k_band * epsilon ** 2):
            continue
        nu_r = np.array([u, rng.uniform(-1, 1), w, 0, 0, 0])
        want = _mss_lift_drag(p, nu_r)
        got = _evaluate(function, nu_r)["tau"].reshape(-1)
        assert np.linalg.norm(got - want) <= 1e-3 * np.linalg.norm(want), (epsilon, nu_r)
        checked += 1


@pytest.mark.parametrize("model", SURGE_MODELS)
def test_surge_damping_jacobian_is_finite_at_rest(model):
    """|u_r| u_r and tanh are smooth at 0; the ITTC log10(Rn + eps) enters
    multiplied by |u_r| u_r (zero slope), so the Jacobian must be finite at
    rest and on the other two singular sets of ``JACOBIAN_STATES`` (pure
    sway, pure rotation)."""
    _, function = _build_surge(_surge_parameters(2, model))
    for state in JACOBIAN_STATES.values():
        assert not _nonfinite_jacobian(function, "tau", state), (model, state)


# --------------------------------------------------------------------------
# Transforms only from more_transformations
# --------------------------------------------------------------------------
OUR_HYDRODYNAMICS_MODULES = (LINEAR_DAMPING, CROSS_FLOW, LIFT_DRAG, SURGE_DAMPING)  # not the other developer's linear_surface
ALLOWED_IMPORT_ROOTS = {"casadi", "numpy", "dataclasses", "typing", "math", "__future__",
                        "more_transformations", "more_dynamics"}
TRANSFORM_NAME = re.compile(r"(skew|smtrx|rzyx|rotation|rot_|hmtrx|euler|quat|gravity|jacobian_|tzyx)")


@pytest.mark.parametrize("module_name", OUR_HYDRODYNAMICS_MODULES)
def test_transforms_no_local_transform_in_hydrodynamics(module_name):
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
# Generic by construction (no vehicle name or number)
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


# builder -> its selectors with a default (vehicle numbers never have one)
GENERIC = {
    LINEAR_DAMPING: {"linear_damping_casadi": {"sway_damping_fade", "smooth_speed"},
                     "linear_damping_parameters": {"smooth_speed"}, "linear_damping_couplings": set()},
    CROSS_FLOW: {"cross_flow_drag_casadi": {"strip_grid"},   # the length flag is gone (owner, 2026-10-06)
                 "cross_flow_drag_parameters": set()},
    LIFT_DRAG: {"lift_drag_casadi": {"smooth_speed"},           # off by default (= MSS)
                "lift_drag_parameters": {"smooth_speed"}},
    SURGE_DAMPING: {"surge_damping_casadi": {"surge_blend", "ittc_reynolds_bound"},  # owner defaults, MSS behind them
                    "surge_damping_parameters": {"ittc_reynolds_bound"}, "surge_damping_couplings": set()},
}


@pytest.mark.parametrize("module_name", sorted(GENERIC))
def test_generic_block_reads_no_vehicle_name_and_has_no_vehicle_defaults(module_name):
    block = _contract(module_name)
    code = _code_without_docstrings(block)
    assert not [n for n in VEHICLE_NAMES if n in code]
    for function_name, flags in GENERIC[module_name].items():
        signature = inspect.signature(_attribute(module_name, function_name))
        defaults = {n for n, p in signature.parameters.items() if p.default is not inspect.Parameter.empty}
        assert defaults == flags, (function_name, defaults)


# --------------------------------------------------------------------------
# Third surge form: the MSS 2.0.2 ship law (osv.m 182-185 at cc07579)
# --------------------------------------------------------------------------
def _osv_surge_reference():
    return _named_columns(MSS_OSV_SURGE)


def _mumford_wetted_surface(length, beam, draft, c_b):
    """XuuITTC.m 38 (Mumford), its two numbers parsed from the pinned line."""
    factor = _number_in(XITTC, 38, r"S = ([0-9.]+) \*")
    draft_factor = _number_in(XITTC, 38, r"\+ ([0-9.]+)\*T")
    return factor * length * (c_b * beam + draft_factor * draft)


def _exp_ittc_parameters(ref, k):
    """Block inputs of one row of the OSV surge file: M(1,1) in an identity
    mass matrix (only M(1,1) is read), the Mumford wetted surface of
    XuuITTC.m 38, the form factor of XuuITTC.m 33, k_u of the row, the
    kinematic viscosity of XuuITTC.m 32 (declared parameter, rule 16)."""
    mass_matrix = np.eye(6)
    mass_matrix[0, 0] = ref["M11"][k]
    return {"quadratic_model": "exp_ittc", "mass_matrix": mass_matrix,
            "wetted_surface": _mumford_wetted_surface(ref["L"][k], ref["B"][k], ref["T"][k], ref["C_B"][k]),
            "length": ref["L"][k], "water_density": ref["rho"][k], "time_constant": ref["T1"][k],
            "form_factor": _value(XITTC, 33), "exponential_decay_rate": ref["k_u"][k],
            "kinematic_viscosity": _value(XITTC, 32)}


def test_G1_MSS_exp_ittc_surge_law_equals_matlab():
    """``surge_form="exp_ittc"`` vs MATLAB on the 220 rows of the OSV surge file
    (MSS cc07579; set 1 = osv.m itself, sets 2-5 seeded, 4 probe rows each):
    X, Xuu, D_nl,11 and D11 = M11/T1 within 1e-9 (default Rn floor = XuuITTC.m
    34); on set 1 also the unmodified osv.m call (``*_first``, = ``*_second``
    since de9a316)."""
    ref = _osv_surge_reference()
    for k in range(len(ref["u_r"])):
        constants, function = _build_surge(_exp_ittc_parameters(ref, k))
        nu_r = np.array([ref["u_r"][k], 0.3, -0.2, 0.1, 0.05, -0.1])   # other DOFs must not enter
        out = _evaluate(function, nu_r)
        assert _max_diff(out["tau"].reshape(-1), [ref["X"][k], 0, 0, 0, 0, 0]) <= G1_TOLERANCE, k
        assert abs(out["quadratic_coefficient"].item() - ref["Xuu"][k]) <= G1_TOLERANCE, k
        assert abs(out["surge_damping_coefficient"].item() - ref["D_nl_11"][k]) <= G1_TOLERANCE, k
        assert abs(-constants.linear_coefficient - ref["D11"][k]) <= G1_TOLERANCE, k
        if ref["set"][k] == 1:
            assert abs(out["surge_damping_coefficient"].item() - ref["D_nl_11_osv_first"][k]) <= G1_TOLERANCE, k
            assert ref["D_nl_11_osv_first"][k] == ref["D_nl_11_osv_second"][k], k


def test_G4_exp_ittc_decay_rate_plus_1_percent_is_detected():
    ref = _osv_surge_reference()
    rows = np.flatnonzero((ref["set"] == 2) & (ref["probe"] == 0))
    worst = lambda scale: max(
        abs(_evaluate(_build_surge({**_exp_ittc_parameters(ref, k),
                                    "exponential_decay_rate": scale * ref["k_u"][k]})[1],
                      np.array([ref["u_r"][k], 0, 0, 0, 0, 0.0]))["tau"][0, 0] - ref["X"][k]) for k in rows)
    assert worst(1.0) <= G1_TOLERANCE
    assert worst(1.01) > G4_FACTOR * G1_TOLERANCE


def test_physical_exp_ittc_surge_law_is_dissipative():
    """Fossen (2011) Property 6.3 (p. 123): X u_r <= 0 on every row and probe
    of the OSV surge file (D_nl,11 >= 0: both terms are non-negative)."""
    ref = _osv_surge_reference()
    for k in range(len(ref["u_r"])):
        _, function = _build_surge(_exp_ittc_parameters(ref, k))
        x = _evaluate(function, np.array([ref["u_r"][k], 0, 0, 0, 0, 0.0]))["tau"][0, 0]
        assert x * ref["u_r"][k] <= 0.0, k


def test_exp_ittc_refuses_the_tanh_blend():
    with pytest.raises(ValueError, match="surge_blend"):
        _contract(SURGE_DAMPING).surge_damping_casadi(surge_form="exp_ittc", surge_blend="mss_tanh")


# --------------------------------------------------------------------------
# The all-CasADi blocks: committed blocks, frozen path, gradients, no numpy
# --------------------------------------------------------------------------
COMMITTED_REVISION = "5d67caa"  # last revision of the blocks with numpy pre-processing
COMMITTED_FILES = ("__init__.py", "linear_surface.py", "cross_flow.py", "lift_drag.py",
                   "linear_damping.py", "surge_damping.py")
COMMITTED_SEED = 20261007
N_COMMITTED_STATES = 300
G2_TOLERANCE = 1e-10
GRADIENT_TOLERANCE = 1e-8


def _committed_package(tmp_path):
    """The hydrodynamics package as committed at ``COMMITTED_REVISION``, read
    with ``git show`` into a temporary package (the old modules are imported,
    no values are frozen in a file). Skips when git or the revision is not
    available (a source archive without history)."""
    root = Path(__file__).resolve().parents[4]
    package = tmp_path / "hydrodynamics_committed"
    package.mkdir()
    for name in COMMITTED_FILES:
        try:
            shown = subprocess.run(
                ["git", "-C", str(root), "show",
                 f"{COMMITTED_REVISION}:more_dynamics/models/hydrodynamics/{name}"],
                capture_output=True, text=True, check=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            pytest.skip(f"revision {COMMITTED_REVISION} not readable with git here: {exc}")
        (package / name).write_text(shown.stdout)
    sys.path.insert(0, str(tmp_path))
    try:
        return SimpleNamespace(**{m: importlib.import_module(f"hydrodynamics_committed.{m}")
                                  for m in ("cross_flow", "lift_drag", "linear_damping", "surge_damping")})
    finally:
        sys.path.remove(str(tmp_path))


def _committed_cases(old):
    """(label, new function of nu_r, old function of nu_r, output names) for
    every selector choice of the four blocks, on small-craft sets (G2's 1e-10
    is absolute; ship scale is held by G1, G2_SCALE_NOTE)."""
    cases = []
    for fade in (False, True):
        for eps in (0.0, 1e-3):
            p = {**_spheroid_damping(fade, MSS_SPHEROID), "smooth_speed_epsilon": eps}
            new = _build_submerged(p)[1]
            oldf = old.linear_damping.submerged_linear_damping_casadi(
                old.linear_damping.preprocess_submerged_linear_damping(**p))
            cases.append((f"submerged fade={fade} eps={eps}", new, oldf, ("D", "tau")))
    for k in (0, 1, 2):
        p = _floating_parameters(k)
        cases.append((f"floating {k}", _build_floating(p)[1],
                      old.linear_damping.floating_linear_damping_casadi(
                          old.linear_damping.preprocess_floating_linear_damping(**p)), ("D", "tau")))
    p = _surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN)
    cases.append(("surface", _build_surface(p)[1],
                  old.linear_damping.surface_linear_damping_casadi(
                      old.linear_damping.preprocess_surface_linear_damping(**p)), ("D", "tau")))
    for make in (_cross_flow_spheroid, _cross_flow_catamaran):
        for grid in ("midpoint", "endpoint"):
            p = make(grid)
            cases.append((f"cross_flow {p['drag_model']} {grid}", _build_cross_flow(p),
                          old.cross_flow.cross_flow_drag_casadi(old.cross_flow.preprocess_cross_flow_drag(**p)),
                          ("tau",)))
    for eps in (0.0, 1e-6, 1e-3):
        p = {**_lift_drag_spheroid(_value(FLD, 26)), "smooth_speed_epsilon": eps}
        cases.append((f"lift_drag eps={eps}", _build_lift_drag(p),
                      old.lift_drag.lift_drag_casadi(old.lift_drag.preprocess_lift_drag(**p)), ("tau",)))
    for set_id in (2, 4):
        for model in SURGE_MODELS:
            for mss in (False, True):
                p = _surge_parameters(set_id, model, mss=mss)
                # kinematic_viscosity and surge_added_mass_factor are declared parameters
                # only from this round (rule 16); COMMITTED_REVISION still fixes them as
                # module constants at the same values, so the old preprocessing function
                # does not take them as keywords.
                q = {k: v for k, v in p.items()
                    if k not in ("quadratic_model", "kinematic_viscosity", "surge_added_mass_factor")}
                oldf = old.surge_damping.surge_damping_casadi(
                    getattr(old.surge_damping, f"preprocess_surge_damping_{model}")(**q))
                cases.append((f"surge {set_id} {model} mss={mss}", _build_surge(p)[1], oldf, ("tau",)))
    return cases


def test_G2_blocks_equal_the_committed_blocks(tmp_path):
    """Every selector choice of linear damping (3 forms, sway fade, smoothing),
    cross-flow (2 models x 2 grids), lift/drag (3 smoothing values) and surge
    damping (2 forms x MSS flags): the outputs equal the blocks committed at
    ``COMMITTED_REVISION`` (numpy pre-processing) on 300 seeded states, 1e-10."""
    old = _committed_package(tmp_path)
    rng = np.random.default_rng(COMMITTED_SEED)
    states = rng.uniform(-3.0, 3.0, size=(N_COMMITTED_STATES, 6))
    cases = _committed_cases(old)
    assert len(cases) == 4 + 3 + 1 + 4 + 3 + 8
    for label, new, oldf, names in cases:
        for k, nu_r in enumerate(states):
            a, b = new(nu_r=nu_r), oldf(nu_r=nu_r)
            for name in names:
                assert _max_diff(np.array(a[name], dtype=float), np.array(b[name], dtype=float)) <= G2_TOLERANCE, \
                    (label, name, k)


def _frozen_cases():
    """(block, declared parameters, couplings + numbers) per module."""
    linear, cross, lift, surge = (_contract(m) for m in (LINEAR_DAMPING, CROSS_FLOW, LIFT_DRAG, SURGE_DAMPING))
    sub = {k: v for k, v in _spheroid_damping(False, MSS_SPHEROID).items() if k != "sway_damping_fade"}
    flo = _floating_parameters(0)
    sur = _surface_damping_catamaran(_value(OTTER, 90), MSS_CATAMARAN)
    cf = {k: v for k, v in _cross_flow_spheroid("midpoint").items() if k not in ("drag_model", "strip_grid")}
    ld = _lift_drag_spheroid(_value(FLD, 26))
    sg = {**{k: v for k, v in _surge_parameters(2, "ittc").items() if k != "quadratic_model"},
          "ittc_reynolds_floor": _value(XITTC, 34)}
    sub_smooth = {**sub, "smooth_speed_epsilon": SMOOTH_EPSILONS[1]}
    ld_smooth = {**ld, "smooth_speed_epsilon": SMOOTH_EPSILONS[1]}
    return {
        "submerged": (linear.linear_damping_casadi(form="submerged"), linear.linear_damping_parameters("submerged"), sub),
        "floating": (linear.linear_damping_casadi(form="floating"), linear.linear_damping_parameters("floating"), flo),
        "surface": (linear.linear_damping_casadi(form="surface"), linear.linear_damping_parameters("surface"), sur),
        "cross_flow": (cross.cross_flow_drag_casadi(drag_model="cylinder"), cross.cross_flow_drag_parameters(), cf),
        "lift_drag": (lift.lift_drag_casadi(), lift.lift_drag_parameters(), ld),
        "surge": (surge.surge_damping_casadi(surge_form="ittc"), surge.surge_damping_parameters("ittc"), sg),
        "submerged_smooth": (linear.linear_damping_casadi(form="submerged", smooth_speed=True),
                             linear.linear_damping_parameters("submerged", smooth_speed=True), sub_smooth),
        "lift_drag_smooth": (lift.lift_drag_casadi(smooth_speed=True),
                             lift.lift_drag_parameters(smooth_speed=True), ld_smooth),
    }


@pytest.mark.parametrize("name", ["submerged", "floating", "surface", "cross_flow", "lift_drag", "surge",
                                  "submerged_smooth", "lift_drag_smooth"])
def test_frozen_block_equals_the_block_called_with_numbers(name):
    """``freeze`` of the declared parameters (what a plugin carries; the
    couplings stay inputs, wired by the vehicle class) against the block
    called with the numbers, 300 seeded states; entry by entry within 1e-14
    relative to max(1, |value|)."""
    from more_transformations.more_casadi_transformations import freeze

    block, declared, values = _frozen_cases()[name]
    parameters = {d.name: values[d.name] for d in declared}
    couplings = {k: v for k, v in values.items() if k not in parameters}
    frozen = freeze(block, declared, parameters)
    assert frozen.name_in() == ["nu_r", *couplings]
    rng = np.random.default_rng(COMMITTED_SEED)
    for k in range(N_COMMITTED_STATES):
        nu_r = rng.uniform(-3.0, 3.0, 6)
        a, b = frozen(nu_r=nu_r, **couplings), block(nu_r=nu_r, **values)
        for output in frozen.name_out():
            unfrozen = np.array(b[output], dtype=float)
            error = np.abs(np.array(a[output], dtype=float) - unfrozen) / np.maximum(1.0, np.abs(unfrozen))
            assert error.max() <= 1e-14, (name, output, k)


# block -> (parameter, central-difference step); tau is linear (floating
# damping ratio, cross-flow and lift/drag density) or quadratic (surge form
# factor enters linearly) in each, so the difference is exact up to rounding
GRADIENT_CASES = {"floating": ("damping_ratios", 1e-2), "cross_flow": ("water_density", 1.0),
                  "lift_drag": ("water_density", 1.0), "surge": ("form_factor", 1e-2)}


@pytest.mark.parametrize("name", sorted(GRADIENT_CASES))
def test_gradient_of_tau_with_respect_to_a_parameter(name):
    """Identification path: the block called with one declared parameter left
    as a symbol; d tau / d parameter from CasADi equals a central difference
    of the block called with numbers, 20 seeded states, 1e-8."""
    import casadi as ca

    block, _, values = _frozen_cases()[name]
    parameter, step = GRADIENT_CASES[name]
    base = np.array(values[parameter], dtype=float).reshape(-1)
    symbol = ca.SX.sym(parameter, base.size)
    rng = np.random.default_rng(COMMITTED_SEED)

    def tau_of(x, nu_r):
        return block(nu_r=nu_r, **{**values, parameter: x})["tau"]

    found = 0.0
    for k in range(20):
        nu_r = rng.uniform(-3.0, 3.0, 6)
        gradient = ca.Function("gradient", [symbol], [ca.jacobian(tau_of(symbol, nu_r), symbol)])
        exact = np.array(gradient(base), dtype=float)
        for i in range(base.size):
            shift = np.zeros(base.size)
            shift[i] = step
            upper = np.array(tau_of(base + shift, nu_r), dtype=float).ravel()
            lower = np.array(tau_of(base - shift, nu_r), dtype=float).ravel()
            assert _max_diff(exact[:, i], (upper - lower) / (2.0 * step)) <= GRADIENT_TOLERANCE, (name, k, i)
        found = max(found, np.abs(exact).max())
    assert found > 1e-2, name  # the gradient is not trivially zero


@pytest.mark.parametrize("module_name", OUR_HYDRODYNAMICS_MODULES)
def test_block_imports_no_numpy(module_name):
    """Our hydrodynamics modules import neither numpy nor scipy, nor the numpy
    ``more_transformations`` modules (AST scan)."""
    for node in ast.walk(ast.parse(Path(_contract(module_name).__file__).read_text())):
        modules = []
        if isinstance(node, ast.Import):
            modules = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            modules = [node.module]
        for module in modules:
            assert module.split(".")[0] not in ("numpy", "scipy"), (module_name, module)
            if module.startswith("more_transformations"):
                assert module.startswith("more_transformations.more_casadi_transformations"), (module_name, module)
