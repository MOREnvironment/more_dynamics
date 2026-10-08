%% ================================================================
%  REMUS 100 (MSS remus100.m) references for the torpedo-AUV vehicle class
%
%  Reference data for more_dynamics/models/vehicles (torpedo-AUV class) on the
%  REMUS 100 parameter set of MSS CRAFT/AUV/models/remus100.m (T. I. Fossen,
%  MIT). MSS is run as written; three marked copies beside this file change
%  only the lines they mark:
%    remus100_hull.m       actuator inputs and forces replaced by an external
%                          generalized force tau_ext (6x1, N and N m, BODY, CO)
%    remus100_munk.m       the zeroed added-mass Coriolis couplings of
%                          remus100.m 207-210 kept (C_A = m2c(M_A, nu_r) in full)
%    remus100_hull_munk.m  both changes
%  Before anything runs, each copy is rebuilt into the original (drop the
%  '%<added>' lines, un-comment the '%<removed> ' lines) and compared with
%  MSS's remus100.m byte for byte; a difference stops the run.
%
%  Outputs (written under <OUT_DIR>/vehicles/remus100/, never this folder):
%    remus100_parameters.json         the parameter set, every value read from
%                                     the MSS workspace or text (no number typed)
%    remus100_derivative_mss.csv      (a) remus100.m xdot, its own actuator
%                                     wrench tau, and (c) remus100_munk.m xdot,
%                                     on structured + seeded states and inputs
%    remus100_hull_derivative_mss.csv stage 1: remus100_hull.m and
%                                     remus100_hull_munk.m xdot on seeded
%                                     states and external wrenches
%    remus100_mass_matrix_mss.csv     [~,~,M] = remus100()
%    remus100_trajectory_<name>.csv   open-loop time histories, RK4 at the
%                                     SIMremus100.m step (h = 0.05 s, rk4.m)
%    remus100_trajectory_tolerances.csv  per scenario and state: the measured
%                                     difference of two integrations (below)
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS checkout,
%  OUT_DIR (optional) = output folder, default <tempdir>/more_mss_references;
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%  Byte-compare each output with the frozen file of the same name here.
%
%  Author:    Enio Krizman
%  Date:      2026-10-07
% ================================================================
clear functions;
format long g
format compact

scriptDir = fileparts(mfilename('fullpath'));
mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout root.');
end
addpath(genpath(mssDir));
addpath(scriptDir);          % the three marked copies
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outDir, 'vehicles', 'remus100');
if ~exist(outDir, 'dir'), mkdir(outDir); end
[~, rev] = system(['git -C "' mssDir '" rev-parse HEAD']);
rev = strtrim(rev);
fprintf('MATLAB %s\nMSS: %s at %s\n', version, which('remus100'), rev);

% ------------------------------------------------
% 0. The copies rebuild the original
% ------------------------------------------------
original = fileread(which('remus100'));
for name = {'remus100_hull', 'remus100_munk', 'remus100_hull_munk'}
    rebuilt = rebuild_original(fullfile(scriptDir, [name{1} '.m']));
    if ~strcmp(rebuilt, original)
        error('%s.m does not rebuild MSS remus100.m at %s', name{1}, rev);
    end
    fprintf('%s.m rebuilds remus100.m byte for byte\n', name{1});
end

remus100_ws = mss_instrumented('remus100');

% ------------------------------------------------
% 1. The parameter set, read from the MSS workspace and text
% ------------------------------------------------
[~, ~, ~, wsR] = remus100_ws(zeros(12, 1), zeros(3, 1), 0, 0, 0);
[~, ~, wsS]  = call_ws('spheroid', 3, wsR.a, wsR.b, zeros(3, 1), wsR.r_bG);
[~, ~, wsI]  = call_ws('imlay61', 3, wsR.a, wsR.b, zeros(6, 1), wsR.r44);
[~, wsLD]    = call_ws('forceLiftDrag', 2, wsR.D_auv, wsR.S, wsR.CD_0, 0, 0);
[~, ~, wsCL] = call_ws('coeffLiftDrag', 3, wsR.D_auv, wsR.S, wsR.CD_0, 0, 0);
[~, wsCF]    = call_ws('crossFlowDrag', 2, wsR.L_auv, wsR.D_auv, wsR.D_auv, zeros(6, 1), 'cylinder');
if ~(wsI.rho == wsLD.rho && wsLD.rho == wsR.rho)
    error('imlay61.m, forceLiftDrag.m and remus100.m water densities differ');
end
% Numbers that remus100.m writes inside an expression, read from its text
txt = strsplit(original, newline, 'CollapseDelimiters', false);
wake = regexp(txt{150}, '^Va = ([0-9.]+) \* U_r;', 'tokens', 'once');
roll = regexp(txt{252}, '^tau\(4\) = K_prop / ([0-9.]+);', 'tokens', 'once');
wag  = regexp(txt{156}, 'wageningen\(0,([0-9.]+),([0-9.]+),([0-9]+)\)', 'tokens', 'once');
if isempty(wake) || isempty(roll) || isempty(wag)
    error('remus100.m lines 150, 156 or 252 moved: re-read the line numbers');
end
oneMinusW = str2double(wake{1});
rollScale = 1 / str2double(roll{1});

P = {};   % name, value, unit, block, declared, MSS place, stage
P(end+1, :) = {'semi_major_axis', wsR.a, 'm', 'rigid_body', 'new', 'remus100.m:134 a = 1.0096 * L_auv/2', 'hull'};
P(end+1, :) = {'semi_minor_axis', wsR.b, 'm', 'rigid_body', 'new', 'remus100.m:135 b = 1.0096 * D_auv/2', 'hull'};
P(end+1, :) = {'body_density', wsS.rho, 'kg/m^3', 'rigid_body', 'new', 'spheroid.m:35 rho = 1025 (called at remus100.m:199)', 'hull'};
P(end+1, :) = {'water_density', wsI.rho, 'kg/m^3', 'rigid_body, lift_drag, propeller, fins', 'new, old, old, old', 'imlay61.m:31 rho = 1026 (remus100.m:200); equal to forceLiftDrag.m:26 and remus100.m:98', 'hull'};
P(end+1, :) = {'roll_added_inertia_ratio', wsR.r44, '1', 'rigid_body', 'new', 'remus100.m:136 r44', 'hull'};
P(end+1, :) = {'body_center_of_gravity', wsR.r_bG, 'm', 'rigid_body', 'new', 'remus100.m:137 r_bG', 'hull'};
P(end+1, :) = {'center_of_buoyancy', wsR.r_bB, 'm', 'submerged_hydrostatics, submerged_linear_damping', 'old, old', 'remus100.m:138 r_bB', 'hull'};
P(end+1, :) = {'weight', wsR.W, 'N', 'submerged_hydrostatics, submerged_linear_damping', 'old, old', 'remus100.m:214 W = m * g_mu (m = MRB(1,1); g_mu = gravity(mu), remus100.m:96-97)', 'hull'};
P(end+1, :) = {'buoyancy', wsR.B, 'N', 'submerged_hydrostatics', 'old', 'remus100.m:214 B = W', 'hull'};
P(end+1, :) = {'time_constants', [wsR.T1; wsR.T2; wsR.T6], 's', 'submerged_linear_damping', 'old', 'remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217)', 'hull'};
P(end+1, :) = {'damping_ratios', [wsR.zeta4; wsR.zeta5], '1', 'submerged_linear_damping', 'old', 'remus100.m:194, 195 [zeta4 zeta5]', 'hull'};
P(end+1, :) = {'span', wsR.D_auv, 'm', 'lift_drag', 'old', 'remus100.m:132 D_auv, passed as b at :220', 'hull'};
P(end+1, :) = {'planform_area', wsR.S, 'm^2', 'lift_drag', 'old', 'remus100.m:133 S = 0.7 * L_auv * D_auv', 'hull'};
P(end+1, :) = {'parasitic_drag_coefficient', wsR.CD_0, '1', 'lift_drag', 'old', 'remus100.m:143-144 CD_0 = Cd * pi * b^2 / S', 'hull'};
P(end+1, :) = {'oswald_efficiency', wsCL.e, '1', 'lift_drag', 'old', 'coeffLiftDrag.m:54 e = 0.3', 'hull'};
P(end+1, :) = {'length', wsR.L_auv, 'm', 'cross_flow', 'old', 'remus100.m:131 L_auv, passed as L at :221', 'hull'};
P(end+1, :) = {'beam', wsR.D_auv, 'm', 'cross_flow', 'old', 'remus100.m:132 D_auv, passed as B at :221', 'hull'};
P(end+1, :) = {'draft', wsR.D_auv, 'm', 'cross_flow', 'old', 'remus100.m:132 D_auv, passed as T at :221', 'hull'};
P(end+1, :) = {'cross_flow_water_density', wsCF.rho, 'kg/m^3', 'cross_flow (its water_density)', 'old', 'crossFlowDrag.m:36 rho = 1025 (remus100.m:221)', 'hull'};
P(end+1, :) = {'diameter', wsR.D_prop, 'm', 'propeller', 'old', 'remus100.m:148 D_prop', 'actuators'};
P(end+1, :) = {'max_speed', wsR.n_max, 'rpm', 'propeller', 'old', 'remus100.m:110 n_max', 'actuators'};
P(end+1, :) = {'thrust_deduction', wsR.t_prop, '1', 'propeller', 'old', 'remus100.m:149 t_prop', 'actuators'};
P(end+1, :) = {'wake_fraction', 1 - oneMinusW, '1', 'propeller', 'old', sprintf('remus100.m:150 Va = %s * U_r, w = 1 - %s', wake{1}, wake{1}), 'actuators'};
P(end+1, :) = {'pitch_diameter_ratio', str2double(wag{1}), '1', 'propeller', 'old', 'remus100.m:156 wageningen(0,1,0.718,3), second argument', 'actuators'};
P(end+1, :) = {'blade_area_ratio', str2double(wag{2}), '1', 'propeller', 'old', 'remus100.m:155-156 blade-area ratio 0.718', 'actuators'};
P(end+1, :) = {'blade_count', str2double(wag{3}), '1', 'propeller', 'old', 'remus100.m:155-156 3 blades', 'actuators'};
P(end+1, :) = {'max_advance_number', wsR.Ja_max, '1', 'propeller', 'old', 'remus100.m:153 Ja_max', 'actuators'};
P(end+1, :) = {'roll_moment_scale', rollScale, '1', 'propeller', 'old', sprintf('remus100.m:252 tau(4) = K_prop / %s', roll{1}), 'actuators'};
P(end+1, :) = {'position', zeros(3, 1), 'm', 'propeller', 'old', 'remus100.m:249-252 thrust on x_b through the CO (no moment arm)', 'actuators'};
P(end+1, :) = {'orientation', zeros(3, 1), 'rad', 'propeller', 'old', 'remus100.m:249, 252 shaft along x_b (thrust in tau(1), torque in tau(4))', 'actuators'};
P(end+1, :) = {'thrust_torque_coefficients', [wsR.KT_0; wsR.KQ_0; wsR.KT_max; wsR.KQ_max], '1', 'propeller', 'old', 'remus100.m:157, 158, 160, 161 [KT_0 KQ_0 KT_max KQ_max]', 'actuators'};
P(end+1, :) = {'rudder_area', wsR.A_r, 'm^2', 'fins', 'old', 'remus100.m:183 A_r = 2 * S_fin (S_fin :179)', 'actuators'};
P(end+1, :) = {'stern_plane_area', wsR.A_s, 'm^2', 'fins', 'old', 'remus100.m:188 A_s = 2 * S_fin (S_fin :179)', 'actuators'};
P(end+1, :) = {'rudder_lift_coefficient', wsR.CL_delta_r, '1/rad', 'fins', 'old', 'remus100.m:182 CL_delta_r', 'actuators'};
P(end+1, :) = {'stern_plane_lift_coefficient', wsR.CL_delta_s, '1/rad', 'fins', 'old', 'remus100.m:187 CL_delta_s', 'actuators'};
P(end+1, :) = {'rudder_position', wsR.x_r, 'm', 'fins', 'old', 'remus100.m:184 x_r = -a', 'actuators'};
P(end+1, :) = {'stern_plane_position', wsR.x_s, 'm', 'fins', 'old', 'remus100.m:189 x_s = -a', 'actuators'};
P(end+1, :) = {'max_deflection', wsR.delta_max, 'rad', 'fins', 'old', 'remus100.m:109 delta_max = deg2rad(20)', 'actuators'};
write_parameters_json(fullfile(outDir, 'remus100_parameters.json'), P, rev, wsR);

% ------------------------------------------------
% 2. Mass matrix (empty call)
% ------------------------------------------------
[~, ~, M0] = remus100();
write_csv(fullfile(outDir, 'remus100_mass_matrix_mss.csv'), ...
    arrayfun(@(j) sprintf('M_%d', j), 1:6, 'UniformOutput', false), M0);

% ------------------------------------------------
% 3. (a) + (c) derivative on structured and seeded cases
% ------------------------------------------------
SEED = 20261007;
x0 = [1.5 0.1 -0.05 0.02 -0.03 0.05 0 0 5 0.05 -0.1 0.3]';
n_list = [-2000 -1525 -800 -1e-3 0 1e-3 800 1525 2000];
d_list = [0 0; deg2rad(20) -deg2rad(20); 0.6 -0.6; -0.2 0.1];
X = []; UI = []; C = [];
for i = 1:numel(n_list)
    for j = 1:size(d_list, 1)
        X = [X x0]; UI = [UI [d_list(j, :) n_list(i)]']; C = [C [0.3; deg2rad(30); 0.1]]; %#ok<AGROW>
    end
end
% at rest, without and with fins and current
X = [X zeros(12, 4)];
UI = [UI [0; 0; 0] [0.2; -0.2; 0] [0; 0; 1000] [0.2; 0.2; 1000]];
C = [C zeros(3, 2) [0.5; deg2rad(30); 0.1] zeros(3, 1)];
rng(SEED, 'twister');
N_RANDOM = 600;
lo = [-1 -1 -1 -0.5 -0.5 -0.5 -10 -10 0 -pi -1.2 -pi]';
hi = [ 3  1  1  0.5  0.5  0.5  10  10 50  pi  1.2  pi]';
Xr = lo + (hi - lo) .* rand(12, N_RANDOM);
UIr = [-0.6; -0.6; -2000] + [1.2; 1.2; 4000] .* rand(3, N_RANDOM);
Cr = [1.0; 2 * pi; 0.4] .* rand(3, N_RANDOM) - [0; pi; 0.2];
Cr(:, 1:2:end) = 0;                     % current off on every other case
X = [X Xr]; UI = [UI UIr]; C = [C Cr];
N = size(X, 2);
rows = zeros(N, 1 + 12 + 3 + 3 + 6 + 12 + 12);
maxSelf = 0; maxSplit = 0;
for k = 1:N
    x = X(:, k); ui = UI(:, k); c = C(:, k);
    xdot = remus100(x, ui, c(1), c(2), c(3));
    [xdotW, ~, ~, ws] = remus100_ws(x, ui, c(1), c(2), c(3));
    maxSelf = max([maxSelf; abs(xdotW - xdot)]);
    xdotSplit = remus100_hull(x, ws.tau, c(1), c(2), c(3));
    maxSplit = max([maxSplit; abs(xdotSplit - xdot)]);
    xdotMunk = remus100_munk(x, ui, c(1), c(2), c(3));
    rows(k, :) = [k x' ui' c' ws.tau' xdot' xdotMunk'];
end
fprintf('(a) %d cases; instrumented vs unmodified max |diff| = %.3g\n', N, maxSelf);
fprintf('split: remus100_hull(x, tau of remus100.m) vs remus100.m, max |diff| = %.3g\n', maxSplit);
if maxSplit > 1e-12, error('the hull split is not exact (%.3g)', maxSplit); end
names = [{'case_id'}, vec_names('x', 12), vec_names('ui', 3), {'Vc', 'betaVc', 'w_c'}, ...
    vec_names('tau', 6), vec_names('xdot', 12), vec_names('xdot_munk', 12)];
write_csv(fullfile(outDir, 'remus100_derivative_mss.csv'), names, rows);

% ------------------------------------------------
% 4. (i) stage 1: hull derivative on seeded states and external wrenches
% ------------------------------------------------
rng(SEED + 1, 'twister');
NH = 600;
Xh = lo + (hi - lo) .* rand(12, NH);
Th = [50; 50; 50; 10; 10; 10] .* (2 * rand(6, NH) - 1);
Ch = [1.0; 2 * pi; 0.4] .* rand(3, NH) - [0; pi; 0.2];
Ch(:, 1:2:end) = 0;
Xh(:, 1) = 0; Th(:, 1) = 0; Ch(:, 1) = 0;   % at rest, no force, no current
rows = zeros(NH, 1 + 12 + 6 + 3 + 12 + 12);
for k = 1:NH
    x = Xh(:, k); t6 = Th(:, k); c = Ch(:, k);
    rows(k, :) = [k x' t6' c' remus100_hull(x, t6, c(1), c(2), c(3))' ...
        remus100_hull_munk(x, t6, c(1), c(2), c(3))'];
end
names = [{'case_id'}, vec_names('x', 12), vec_names('tau_ext', 6), {'Vc', 'betaVc', 'w_c'}, ...
    vec_names('xdot', 12), vec_names('xdot_munk', 12)];
write_csv(fullfile(outDir, 'remus100_hull_derivative_mss.csv'), names, rows);

% ------------------------------------------------
% 5. Trajectories: rk4.m at h = 0.05 s (SIMremus100.m 60, 348)
% ------------------------------------------------
h = 0.05;               % SIMremus100.m 60
HOLD = 5;               % inputs decided and held over 5 steps (0.25 s); a row every 5 steps
OFFSET = 1e-9;          % the G1 tolerance, added to every derivative entry
deg = pi / 180;
S = struct([]);
S = add(S, 'straight_from_rest', 'full', zeros(12, 1), [0 0 0], 30, @(t, x, s) [0 0 1300]);
S = add(S, 'rudder_step_20deg', 'full', [1.5 zeros(1, 11)]', [0 0 0], 40, ...
    @(t, x, s) [20 * deg * (t >= 5) 0 1300]);
% Prestero (2001) p. 53: rudder roughly 4 deg; Table 8.1 p. 60: u = 1.54 m/s.
% The propeller speed that holds 1.54 m/s straight (surge acceleration zero
% at that state) is found with fzero on remus100.m itself.
U_PRESTERO = 1.54;
surge = @(n) subsref(remus100([U_PRESTERO zeros(1, 11)]', [0 0 n]', 0, 0, 0), substruct('()', {1}));
nTrim = fzero(surge, [100 1525], optimset('TolX', 1e-12));
fprintf('propeller speed holding %.2f m/s straight: %.12g rpm (surge acceleration %.3g)\n', ...
    U_PRESTERO, nTrim, surge(nTrim));
S = add(S, 'rudder_step_4deg_at_1p54ms', 'full', [U_PRESTERO zeros(1, 11)]', [0 0 0], 60, ...
    @(t, x, s) [4 * deg * (t >= 10) 0 nTrim]);
S = add(S, 'stern_plane_step_10deg', 'full', [1.5 zeros(1, 11)]', [0 0 0], 30, ...
    @(t, x, s) [0 10 * deg * (t >= 5) 1300]);
S = add(S, 'zigzag_20_20', 'full', [1.5 zeros(1, 11)]', [0 0 0], 60, @zigzag);
S = add(S, 'dive_and_level_off', 'full', [1.5 zeros(1, 11)]', [0 0 0], 40, @dive_level);
S = add(S, 'rudder_step_10deg_in_current', 'full', [1.5 zeros(1, 11)]', [0.5 30 * deg 0.1], 30, ...
    @(t, x, s) [10 * deg * (t >= 5) 0 1300]);
S = add(S, 'munk_straight_perturbed', 'full_munk', [1.5 0.01 0 0 0 0.01 zeros(1, 6)]', [0 0 0], 30, ...
    @(t, x, s) [0 0 1300]);
S = add(S, 'munk_rudder_step_20deg', 'full_munk', [1.5 zeros(1, 11)]', [0 0 0], 40, ...
    @(t, x, s) [20 * deg * (t >= 5) 0 1300]);
S = add(S, 'hull_glide', 'hull', [2 zeros(1, 11)]', [0 0 0], 30, @(t, x, s) zeros(1, 6));
S = add(S, 'hull_sway_yaw_release', 'hull', [1 0.3 0 0 0 0.2 zeros(1, 6)]', [0 0 0], 20, ...
    @(t, x, s) zeros(1, 6));
S = add(S, 'hull_roll_pitch_decay', 'hull', [1 0 0 0 0 0 0 0 0 20 * deg 15 * deg 0]', [0 0 0], 20, ...
    @(t, x, s) zeros(1, 6));
S = add(S, 'hull_munk_glide', 'hull_munk', [2 0.01 0 0 0 0.01 zeros(1, 6)]', [0 0 0], 20, ...
    @(t, x, s) zeros(1, 6));
S = add(S, 'hull_zoh_rudder_step_20deg', 'hull_zoh', [1.5 zeros(1, 11)]', [0 0 0], 40, ...
    @(t, x, s) [20 * deg * (t >= 5) 0 1300]);
S = add(S, 'hull_zoh_stern_plane_step_10deg', 'hull_zoh', [1.5 zeros(1, 11)]', [0 0 0], 30, ...
    @(t, x, s) [0 10 * deg * (t >= 5) 1300]);
S = add(S, 'hull_zoh_zigzag_20_20', 'hull_zoh', [1.5 zeros(1, 11)]', [0 0 0], 60, @zigzag);

tolRows = {};
for s = 1:numel(S)
    sc = S(s);
    nSteps = round(sc.T / h);
    nRows = nSteps / HOLD + 1;
    x = sc.x0; c = sc.current;
    state = struct('phase', 0, 'psi0', sc.x0(12));
    zohAll = {};
    nIn = 3 + 3 * strcmp(sc.kind, 'hull_zoh'); if any(strcmp(sc.kind, {'hull', 'hull_munk'})), nIn = 6; end
    inputs = zeros(nRows, nIn); states = zeros(nRows, 12); times = (0:nRows - 1)' * HOLD * h;
    for r = 1:nRows
        t = times(r);
        states(r, :) = x';
        if r == nRows, inputs(r, :) = inputs(r - 1, :); break; end
        [u, state] = decide(sc.law, t, x, state);
        switch sc.kind
            case 'full'
                inputs(r, :) = u;
                for k = 1:HOLD, x = rk4(@remus100, h, x, u', c(1), c(2), c(3)); end
            case 'full_munk'
                inputs(r, :) = u;
                for k = 1:HOLD, x = rk4(@remus100_munk, h, x, u', c(1), c(2), c(3)); end
            case 'hull'
                inputs(r, :) = u;
                for k = 1:HOLD, x = rk4(@remus100_hull, h, x, u', c(1), c(2), c(3)); end
            case 'hull_munk'
                inputs(r, :) = u;
                for k = 1:HOLD, x = rk4(@remus100_hull_munk, h, x, u', c(1), c(2), c(3)); end
            case 'hull_zoh'
                % the actuator wrench of remus100.m at the start of each step,
                % held over that step
                zoh = zeros(HOLD, 6);
                for k = 1:HOLD
                    [~, ~, ~, ws] = remus100_ws(x, u', c(1), c(2), c(3));
                    zoh(k, :) = ws.tau';
                    x = rk4(@remus100_hull, h, x, ws.tau, c(1), c(2), c(3));
                end
                inputs(r, :) = zoh(1, :);
                zohAll{r} = zoh; %#ok<AGROW>
        end
    end
    if strcmp(sc.kind, 'hull_zoh')
        % one row per step so the wrench of each step is recorded
        [inputs, states, times] = expand_zoh(zohAll, states, times, HOLD, h, sc, c);
    end
    % the measured difference of two integrations: (1) the derivative with
    % +-OFFSET on every entry, (2) the same RK4 with its sum reordered
    dev = zeros(3, 12);
    variants = {+OFFSET, -OFFSET, 'reordered'};
    for v = 1:3
        xa = sc.x0; dmax = zeros(1, 12);
        stride = size(states, 1) - 1; stepsPerRow = nSteps / stride;
        for r = 1:size(states, 1)
            dmax = max(dmax, abs(xa' - states(r, :)));
            if r == size(states, 1), break; end
            u = inputs(r, :);
            for k = 1:stepsPerRow
                f = model_of(sc.kind);
                if ischar(variants{v})
                    xa = rk4_reordered(f, h, xa, u', c);
                else
                    g = @(xx, uu, a1, a2, a3) f(xx, uu, a1, a2, a3) + variants{v};
                    xa = rk4(g, h, xa, u', c(1), c(2), c(3));
                end
            end
        end
        dev(v, :) = dmax;
    end
    for j = 1:12
        tolRows(end + 1, :) = {sc.name, j, dev(1, j), dev(2, j), dev(3, j), max(dev(:, j))}; %#ok<AGROW>
    end
    inNames = vec_names('ui', 3);
    if any(strcmp(sc.kind, {'hull', 'hull_munk', 'hull_zoh'})), inNames = vec_names('tau_ext', 6); end
    write_csv(fullfile(outDir, ['remus100_trajectory_' sc.name '.csv']), ...
        [{'t'}, inNames, {'Vc', 'betaVc', 'w_c'}, vec_names('x', 12)], ...
        [times inputs repmat(c, size(states, 1), 1) states]);
    if strcmp(sc.name, 'rudder_step_4deg_at_1p54ms')
        last = times >= times(end) - 10;
        fprintf('steady yaw rate, last 10 s: %.6g deg/s at mean speed %.6g m/s\n', ...
            mean(states(last, 6)) * 180 / pi, mean(sqrt(sum(states(last, 1:3).^2, 2))));
    end
    fprintf('%-34s %-9s %5d rows, end u = %.4f m/s, psi = %.4f rad, z = %.4f m\n', ...
        sc.name, sc.kind, size(states, 1), states(end, 1), states(end, 12), states(end, 9));
end
fid = fopen(fullfile(outDir, 'remus100_trajectory_tolerances.csv'), 'w');
fprintf(fid, 'scenario,state,offset_plus,offset_minus,reordered,tolerance\n');
for i = 1:size(tolRows, 1)
    fprintf(fid, '%s,%d,%.17g,%.17g,%.17g,%.17g\n', tolRows{i, :});
end
fclose(fid);
fprintf('Saved to %s\n', outDir);

%% ================================================================
% Scenario helpers
% ================================================================
function S = add(S, name, kind, x0, current, T, law)
    sc = struct('name', name, 'kind', kind, 'x0', x0, 'current', current, 'T', T, 'law', law);
    if isempty(S), S = sc; else, S(end + 1) = sc; end
end

function [u, state] = decide(law, t, x, state)
    if nargin(law) == 3 && ~isempty(regexp(func2str(law), '^(zigzag|dive_level)$', 'once'))
        [u, state] = law(t, x, state);
    else
        u = law(t, x, state);
    end
end

function [u, state] = zigzag(t, x, state)
% 20/20 zig-zag: 1300 rpm; rudder +20 deg from t = 5 s; the rudder reverses
% when the heading has moved 20 deg past the initial heading (decided at the
% 0.25 s input ticks)
    deg = pi / 180;
    if t < 5
        u = [0 0 1300]; return
    end
    if state.phase == 0, state.phase = 1; end
    if state.phase == 1 && x(12) - state.psi0 >= 20 * deg, state.phase = -1; end
    if state.phase == -1 && x(12) - state.psi0 <= -20 * deg, state.phase = 1; end
    u = [state.phase * 20 * deg 0 1300];
end

function [u, state] = dive_level(t, x, state)
% 1300 rpm; stern plane +15 deg (nose down) from t = 5 s until the depth is
% 3 m, then -15 deg until the pitch angle is back to zero, then 0 (decided at
% the 0.25 s input ticks)
    deg = pi / 180;
    if t < 5, u = [0 0 1300]; return, end
    if state.phase == 0, state.phase = 1; end
    if state.phase == 1 && x(9) >= 3, state.phase = 2; end
    if state.phase == 2 && x(11) >= 0, state.phase = 3; end
    d = [15 -15 0];
    u = [0 d(state.phase) * deg 1300];
end

function f = model_of(kind)
    switch kind
        case 'full',       f = @remus100;
        case 'full_munk',  f = @remus100_munk;
        case 'hull_munk',  f = @remus100_hull_munk;
        otherwise,         f = @remus100_hull;
    end
end

function [inputs, states, times] = expand_zoh(zohAll, statesCoarse, timesCoarse, HOLD, h, sc, c)
% Re-run the held-wrench integration and keep one row per step (the wrench
% of that step in the row), so a reader can replay it step by step.
    nCoarse = size(statesCoarse, 1);
    inputs = zeros((nCoarse - 1) * HOLD + 1, 6);
    states = zeros((nCoarse - 1) * HOLD + 1, 12);
    x = sc.x0; r = 0;
    for i = 1:nCoarse - 1
        for k = 1:HOLD
            r = r + 1;
            states(r, :) = x'; inputs(r, :) = zohAll{i}(k, :);
            x = rk4(@remus100_hull, h, x, zohAll{i}(k, :)', c(1), c(2), c(3));
        end
    end
    states(end, :) = x'; inputs(end, :) = inputs(end - 1, :);
    if max(abs(x - statesCoarse(end, :)')) ~= 0
        error('%s: step-wise replay differs from the run', sc.name);
    end
    times = (0:size(states, 1) - 1)' * h;
end

function x = rk4_reordered(f, h, x, u, c)
% The same Runge-Kutta 4 as rk4.m, the weighted sum written in another order
    k1 = f(x, u, c(1), c(2), c(3));
    k2 = f(x + 0.5 * h * k1, u, c(1), c(2), c(3));
    k3 = f(x + 0.5 * h * k2, u, c(1), c(2), c(3));
    k4 = f(x + h * k3, u, c(1), c(2), c(3));
    x = x + h * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6);
end

%% ================================================================
% I/O helpers
% ================================================================
function names = vec_names(prefix, n)
    names = arrayfun(@(i) sprintf('%s_%02d', prefix, i), 1:n, 'UniformOutput', false);
end

function write_csv(file, names, data)
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', strjoin(names, ','));
    fmt = [strjoin(repmat({'%.17g'}, 1, size(data, 2)), ','), '\n'];
    fprintf(fid, fmt, data');
    fclose(fid);
end

function write_parameters_json(file, P, rev, wsR)
    fid = fopen(file, 'w');
    fprintf(fid, '{\n');
    fprintf(fid, '  "vehicle": "REMUS 100, the MSS 1.6 m REMUS-like model",\n');
    fprintf(fid, '  "source": {"repository": "MSS (T. I. Fossen), MIT licence", "revision": "%s", "file": "CRAFT/AUV/models/remus100.m"},\n', rev);
    fprintf(fid, '  "gravity": {"latitude_rad": %.17g, "g": %.17g, "place": "remus100.m:96-97 mu = deg2rad(63.446827), g_mu = gravity(mu); weight = mass * g"},\n', wsR.mu, wsR.g_mu);
    fprintf(fid, '  "selectors": {"mass_properties": "spheroid", "coriolis": "co", "drag_model": "cylinder", "strip_grid": "midpoint", "sway_damping_fade": false, "smooth_speed_epsilon": 0.0, "open_water_model": "linearized", "clip_advance_ratio": false, "convention": "starboard_down_positive"},\n');
    fprintf(fid, '  "parameters": {\n');
    for i = 1:size(P, 1)
        v = P{i, 2};
        if numel(v) == 1
            vs = sprintf('%.17g', v);
        else
            vs = ['[' strjoin(arrayfun(@(a) sprintf('%.17g', a), v(:)', 'UniformOutput', false), ', ') ']'];
        end
        sep = ','; if i == size(P, 1), sep = ''; end
        % the name the block itself declares (differs only where the class
        % must keep two values MSS uses for one block quantity apart)
        inner = P{i, 1};
        if strcmp(inner, 'cross_flow_water_density'), inner = 'water_density'; end
        fprintf(fid, '    "%s": {"value": %s, "unit": "%s", "block": "%s", "block_parameter": "%s", "declared": "%s", "place": "%s", "stage": "%s", "kind": "published"}%s\n', ...
            P{i, 1}, vs, P{i, 3}, P{i, 4}, inner, P{i, 5}, P{i, 6}, P{i, 7}, sep);
    end
    fprintf(fid, '  }\n}\n');
    fclose(fid);
end

function text = rebuild_original(copyFile)
    lines = strsplit(fileread(copyFile), newline, 'CollapseDelimiters', false);
    keep = {};
    for i = 1:numel(lines)
        t = lines{i};
        if endsWith(t, '%<added>'), continue, end
        if startsWith(t, '%<removed> '), t = t(numel('%<removed> ') + 1:end); end
        keep{end + 1} = t; %#ok<AGROW>
    end
    text = strjoin(keep, newline);
end

function varargout = call_ws(name, nOut, varargin)
% Call an MSS function through an instrumented copy and return its outputs
% and, last, its workspace.
    varargout = cell(1, nOut);
    try
        fh = mss_instrumented(name, 'before_last_end');
        [varargout{:}] = fh(varargin{:});
    catch
        % the file's last 'end' closes a block, not the function
        fh = mss_instrumented(name, 'append');
        [varargout{:}] = fh(varargin{:});
    end
end

function fh = mss_instrumented(name, mode)
% Copy of the MSS file <name>.m, written to tempdir at run time, that also
% returns its workspace as a last output. Only the function line changes and
% one capture line is added before the final 'end' (or at the end of the
% file when that 'end' closes a block); the body is MSS's text.
    src = fileread(which(name));
    hdr = ['function\s*(\[([^\]]*)\]|(\w+))\s*=\s*' name '\s*\('];
    tok = regexp(src, hdr, 'tokens', 'once');
    if isempty(tok)
        error('%s.m: function line not found', name);
    end
    outs = regexprep(tok{1}, '[\[\]]', '');
    src = regexprep(src, hdr, ['function [' outs ',ws__] = ' name '_ws('], 'once');
    if numel(regexp(src, '(^|\n)\s*function\s')) ~= 1
        error('%s.m: expected exactly one function', name);
    end
    capture = [newline 'ws__ = struct(); vl__ = who; ' ...
        'for i__ = 1:numel(vl__), ws__.(vl__{i__}) = eval(vl__{i__}); end' newline];
    if nargin < 2, mode = 'before_last_end'; end
    k = regexp(src, '\n[ \t]*end[ \t\r\n]*$');
    if isempty(k) || strcmp(mode, 'append')
        src = [src capture];
    else
        src = [src(1:k-1) capture src(k:end)];
    end
    d = fullfile(tempdir, 'mss_current_ws');
    if ~exist(d, 'dir'), mkdir(d); end
    fid = fopen(fullfile(d, [name '_ws.m']), 'w');
    fwrite(fid, src);
    fclose(fid);
    addpath(d);
    rehash;
    fh = str2func([name '_ws']);
end
