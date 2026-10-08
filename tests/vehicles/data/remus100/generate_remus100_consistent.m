%% ================================================================
%  REMUS 100 references with one value per physical quantity
%
%  Reference data for more_dynamics/models/vehicles (torpedo-AUV class) on the
%  REMUS 100 parameter set of MSS CRAFT/AUV/models/remus100.m (T. I. Fossen,
%  MIT), with every quantity holding one value in the whole model:
%    - one water density, 1026 kg/m^3 (remus100.m 98), in the added mass, the
%      hull lift and drag, the cross-flow, the propeller and the fins
%      (remus100.m runs spheroid.m at 1025, imlay61.m at 1026,
%      forceLiftDrag.m at 1026 and crossFlowDrag.m at 1025);
%    - one geometry, the hull's L = 1.6 m and D = 0.19 m (remus100.m 131-132),
%      also for the mass, the added mass, CD_0 and the fin positions
%      (remus100.m 134-135 scale both by 1.0096 for those only);
%    - the mass given, 31.9 kg (remus100.m 3), W = m g_mu and B = W
%      (remus100.m 214).
%  The model is the marked copy remus100_consistent.m beside this file, with
%  marked copies spheroid_consistent.m, imlay61_consistent.m,
%  forceLiftDrag_consistent.m and crossFlowDrag_consistent.m (each takes the
%  mass or the density as an input). Before anything runs, every copy is
%  rebuilt into its original (drop the '%<added>' lines, un-comment the
%  '%<removed> ' lines) and compared with MSS byte for byte.
%
%  Outputs (written under <OUT_DIR>/vehicles/remus100/, never this folder):
%    remus100_parameters_consistent.json      the parameter set
%    remus100_mass_matrix_consistent.csv      [~,~,M] = remus100_consistent()
%    remus100_derivative_consistent.csv       the cases of
%        remus100_derivative_mss.csv: xdot, its actuator wrench tau,
%        xdot_munk (couplings of remus100.m 207-210 kept) and
%        xdot_one_density (MSS geometry and mass rule, one density: history)
%    remus100_hull_derivative_consistent.csv  the cases of
%        remus100_hull_derivative_mss.csv: xdot, xdot_munk
%    remus100_trajectory_<name>_consistent.csv  the scenarios of
%        generate_remus100_mss.m plus rudder_step_15deg_at_925rpm
%    remus100_trajectory_tolerances_consistent.csv
%    remus100_trajectory_rudder_step_15deg_at_925rpm_{mss,one_density}.csv
%        the same turn by unmodified remus100.m and by the one-density
%        history model (no gate; they show the size of the change)
%  Printed: the density invariance of the copy, the hull split, the
%  differences to MSS.
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
addpath(scriptDir);          % the marked copies
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outDir, 'vehicles', 'remus100');
if ~exist(outDir, 'dir'), mkdir(outDir); end
[~, rev] = system(['git -C "' mssDir '" rev-parse HEAD']);
rev = strtrim(rev);
fprintf('MATLAB %s\nMSS: %s at %s\n', version, which('remus100'), rev);

% ------------------------------------------------
% 0. The copies rebuild their originals; the given mass is line 3's
% ------------------------------------------------
copies = {'remus100_consistent', 'remus100'; 'spheroid_consistent', 'spheroid'; ...
          'imlay61_consistent', 'imlay61'; 'forceLiftDrag_consistent', 'forceLiftDrag'; ...
          'crossFlowDrag_consistent', 'crossFlowDrag'};
for i = 1:size(copies, 1)
    rebuilt = rebuild_original(fullfile(scriptDir, [copies{i, 1} '.m']));
    if ~strcmp(rebuilt, fileread(which(copies{i, 2})))
        error('%s.m does not rebuild MSS %s.m at %s', copies{i, 1}, copies{i, 2}, rev);
    end
    fprintf('%s.m rebuilds %s.m byte for byte\n', copies{i, 1}, copies{i, 2});
end
original = fileread(which('remus100'));
txt = strsplit(original, newline, 'CollapseDelimiters', false);
massText = regexp(txt{3}, 'mass of the vehicle is ([0-9.]+) kg', 'tokens', 'once');
copyText = regexp(fileread(fullfile(scriptDir, 'remus100_consistent.m')), ...
    '\nm_given = ([0-9.]+);', 'tokens', 'once');
if isempty(massText) || isempty(copyText) || ~strcmp(massText{1}, copyText{1})
    error('the given mass of remus100_consistent.m is not the mass of remus100.m line 3');
end
fprintf('given mass %s kg = remus100.m line 3\n', massText{1});

remusC_ws = mss_instrumented('remus100_consistent');
full = @(x, u, a1, a2, a3) remus100_consistent(x, u, a1, a2, a3);
full_munk = @(x, u, a1, a2, a3) remus100_consistent(x, u, a1, a2, a3, struct('munk', true));
hull = @(x, t, a1, a2, a3) remus100_consistent(x, zeros(3, 1), a1, a2, a3, struct('tau_ext', t));
hull_munk = @(x, t, a1, a2, a3) remus100_consistent(x, zeros(3, 1), a1, a2, a3, struct('tau_ext', t, 'munk', true));
one_density = @(x, u, a1, a2, a3) remus100_consistent(x, u, a1, a2, a3, struct('geometry', 'mss'));

% ------------------------------------------------
% 1. The parameter set, read from the workspace of the copy and MSS text
% ------------------------------------------------
[~, ~, ~, ~, wsR] = remusC_ws(zeros(12, 1), zeros(3, 1), 0, 0, 0);
[~, ~, wsCL] = call_ws('coeffLiftDrag', 3, wsR.D_auv, wsR.S, wsR.CD_0, 0, 0);
% Numbers that remus100.m writes inside an expression, read from its text
wake = regexp(txt{150}, '^Va = ([0-9.]+) \* U_r;', 'tokens', 'once');
roll = regexp(txt{252}, '^tau\(4\) = K_prop / ([0-9.]+);', 'tokens', 'once');
wag  = regexp(txt{156}, 'wageningen\(0,([0-9.]+),([0-9.]+),([0-9]+)\)', 'tokens', 'once');
if isempty(wake) || isempty(roll) || isempty(wag)
    error('remus100.m lines 150, 156 or 252 moved: re-read the line numbers');
end
oneMinusW = str2double(wake{1});
rollScale = 1 / str2double(roll{1});
volume = 4/3 * pi * wsR.a * wsR.b^2;
bodyDensity = wsR.m_given / volume;
if abs(wsR.MRB(1, 1) - wsR.m_given) > 0
    error('the copy''s rigid-body mass is not the given mass');
end
fprintf('one geometry: a = %.17g m, b = %.17g m; volume %.17g m^3; body density = m/volume = %.17g kg/m^3\n', ...
    wsR.a, wsR.b, volume, bodyDensity);
fprintf('mass %.17g kg, displaced water at rho = %g: %.17g kg; W = B = %.17g N\n', ...
    wsR.m_given, wsR.rho, wsR.rho * volume, wsR.W);

P = {};   % name, value, unit, block, MSS place, stage, kind
P(end+1, :) = {'semi_major_axis', wsR.a, 'm', 'rigid_body', 'remus100.m:131 L_auv / 2 (one geometry; MSS scales by 1.0096 at :134)', 'hull', 'published'};
P(end+1, :) = {'semi_minor_axis', wsR.b, 'm', 'rigid_body', 'remus100.m:132 D_auv / 2 (one geometry; MSS scales by 1.0096 at :135)', 'hull', 'published'};
P(end+1, :) = {'body_density', bodyDensity, 'kg/m^3', 'rigid_body', 'mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: 1025 at spheroid.m:35', 'hull', 'derived'};
P(end+1, :) = {'water_density', wsR.rho, 'kg/m^3', 'rigid_body, lift_drag, cross_flow, propeller, fins', 'remus100.m:98 rho = 1026, the one value (MSS: 1026 at imlay61.m:31 and forceLiftDrag.m:26, 1025 at crossFlowDrag.m:36)', 'hull', 'published'};
P(end+1, :) = {'roll_added_inertia_ratio', wsR.r44, '1', 'rigid_body', 'remus100.m:136 r44', 'hull', 'published'};
P(end+1, :) = {'body_center_of_gravity', wsR.r_bG, 'm', 'rigid_body', 'remus100.m:137 r_bG', 'hull', 'published'};
P(end+1, :) = {'center_of_buoyancy', wsR.r_bB, 'm', 'submerged_hydrostatics, submerged_linear_damping', 'remus100.m:138 r_bB', 'hull', 'published'};
P(end+1, :) = {'weight', wsR.W, 'N', 'submerged_hydrostatics, submerged_linear_damping', 'remus100.m:214 W = m * g_mu, m = 31.9 kg given (remus100.m:3), g_mu = gravity(mu) (remus100.m:96-97)', 'hull', 'derived'};
P(end+1, :) = {'buoyancy', wsR.B, 'N', 'submerged_hydrostatics', 'remus100.m:214 B = W (neutral, as REMUS states)', 'hull', 'derived'};
P(end+1, :) = {'time_constants', [wsR.T1; wsR.T2; wsR.T6], 's', 'submerged_linear_damping', 'remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217)', 'hull', 'published'};
P(end+1, :) = {'damping_ratios', [wsR.zeta4; wsR.zeta5], '1', 'submerged_linear_damping', 'remus100.m:194, 195 [zeta4 zeta5]', 'hull', 'published'};
P(end+1, :) = {'span', wsR.D_auv, 'm', 'lift_drag', 'remus100.m:132 D_auv, passed as b at :220', 'hull', 'published'};
P(end+1, :) = {'planform_area', wsR.S, 'm^2', 'lift_drag', 'remus100.m:133 S = 0.7 * L_auv * D_auv', 'hull', 'derived'};
P(end+1, :) = {'parasitic_drag_coefficient', wsR.CD_0, '1', 'lift_drag', 'remus100.m:143-144 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry)', 'hull', 'derived'};
P(end+1, :) = {'oswald_efficiency', wsCL.e, '1', 'lift_drag', 'coeffLiftDrag.m:54 e = 0.3', 'hull', 'published'};
P(end+1, :) = {'length', wsR.L_auv, 'm', 'cross_flow', 'remus100.m:131 L_auv, passed as L at :221', 'hull', 'published'};
P(end+1, :) = {'beam', wsR.D_auv, 'm', 'cross_flow', 'remus100.m:132 D_auv, passed as B at :221', 'hull', 'published'};
P(end+1, :) = {'draft', wsR.D_auv, 'm', 'cross_flow', 'remus100.m:132 D_auv, passed as T at :221', 'hull', 'published'};
P(end+1, :) = {'diameter', wsR.D_prop, 'm', 'propeller', 'remus100.m:148 D_prop', 'actuators', 'published'};
P(end+1, :) = {'max_speed', wsR.n_max, 'rpm', 'propeller', 'remus100.m:110 n_max', 'actuators', 'published'};
P(end+1, :) = {'thrust_deduction', wsR.t_prop, '1', 'propeller', 'remus100.m:149 t_prop', 'actuators', 'published'};
P(end+1, :) = {'wake_fraction', 1 - oneMinusW, '1', 'propeller', sprintf('remus100.m:150 Va = %s * U_r, w = 1 - %s', wake{1}, wake{1}), 'actuators', 'published'};
P(end+1, :) = {'pitch_diameter_ratio', str2double(wag{1}), '1', 'propeller', 'remus100.m:156 wageningen(0,1,0.718,3), second argument', 'actuators', 'published'};
P(end+1, :) = {'blade_area_ratio', str2double(wag{2}), '1', 'propeller', 'remus100.m:155-156 blade-area ratio 0.718', 'actuators', 'published'};
P(end+1, :) = {'blade_count', str2double(wag{3}), '1', 'propeller', 'remus100.m:155-156 3 blades', 'actuators', 'published'};
P(end+1, :) = {'max_advance_number', wsR.Ja_max, '1', 'propeller', 'remus100.m:153 Ja_max', 'actuators', 'published'};
P(end+1, :) = {'roll_moment_scale', rollScale, '1', 'propeller', sprintf('remus100.m:252 tau(4) = K_prop / %s', roll{1}), 'actuators', 'published'};
P(end+1, :) = {'position', zeros(3, 1), 'm', 'propeller', 'remus100.m:249-252 thrust on x_b through the CO (no moment arm)', 'actuators', 'published'};
P(end+1, :) = {'orientation', zeros(3, 1), 'rad', 'propeller', 'remus100.m:249, 252 shaft along x_b (thrust in tau(1), torque in tau(4))', 'actuators', 'published'};
P(end+1, :) = {'thrust_torque_coefficients', [wsR.KT_0; wsR.KQ_0; wsR.KT_max; wsR.KQ_max], '1', 'propeller', 'remus100.m:157, 158, 160, 161 [KT_0 KQ_0 KT_max KQ_max]', 'actuators', 'published'};
P(end+1, :) = {'rudder_area', wsR.A_r, 'm^2', 'fins', 'remus100.m:183 A_r = 2 * S_fin (S_fin :179)', 'actuators', 'published'};
P(end+1, :) = {'stern_plane_area', wsR.A_s, 'm^2', 'fins', 'remus100.m:188 A_s = 2 * S_fin (S_fin :179)', 'actuators', 'published'};
P(end+1, :) = {'rudder_lift_coefficient', wsR.CL_delta_r, '1/rad', 'fins', 'remus100.m:182 CL_delta_r', 'actuators', 'published'};
P(end+1, :) = {'stern_plane_lift_coefficient', wsR.CL_delta_s, '1/rad', 'fins', 'remus100.m:187 CL_delta_s', 'actuators', 'published'};
P(end+1, :) = {'rudder_position', wsR.x_r, 'm', 'fins', 'remus100.m:184 x_r = -a, a = L_auv / 2 (one geometry)', 'actuators', 'derived'};
P(end+1, :) = {'stern_plane_position', wsR.x_s, 'm', 'fins', 'remus100.m:189 x_s = -a, a = L_auv / 2 (one geometry)', 'actuators', 'derived'};
P(end+1, :) = {'max_deflection', wsR.delta_max, 'rad', 'fins', 'remus100.m:109 delta_max = deg2rad(20)', 'actuators', 'published'};
write_parameters_json(fullfile(outDir, 'remus100_parameters_consistent.json'), P, rev, wsR, str2double(massText{1}));

% ------------------------------------------------
% 2. Mass matrix (empty call)
% ------------------------------------------------
[~, ~, M0] = remus100_consistent();
write_csv(fullfile(outDir, 'remus100_mass_matrix_consistent.csv'), ...
    arrayfun(@(j) sprintf('M_%d', j), 1:6, 'UniformOutput', false), M0);

% ------------------------------------------------
% 3. Derivative on the cases of remus100_derivative_mss.csv (same seed)
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
Cr(:, 1:2:end) = 0;
X = [X Xr]; UI = [UI UIr]; C = [C Cr];
N = size(X, 2);
mssTable = readmatrix(fullfile(scriptDir, 'remus100_derivative_mss.csv'));
if max(max(abs(mssTable(:, 2:13) - X'))) ~= 0 || max(max(abs(mssTable(:, 14:16) - UI'))) ~= 0 ...
        || max(max(abs(mssTable(:, 17:19) - C'))) ~= 0
    error('the cases differ from remus100_derivative_mss.csv');
end
rows = zeros(N, 1 + 12 + 3 + 3 + 6 + 12 + 12 + 12);
maxSplit = 0; invC = zeros(1, 3); invMss = 0;
dens = [0 0]; geom = [0 0]; total = [0 0];
for k = 1:N
    x = X(:, k); ui = UI(:, k); c = C(:, k);
    [xdot, ~, ~, tau] = remus100_consistent(x, ui, c(1), c(2), c(3));
    xdotSplit = hull(x, tau, c(1), c(2), c(3));
    maxSplit = max([maxSplit; abs(xdotSplit - xdot)]);
    xdotMunk = full_munk(x, ui, c(1), c(2), c(3));
    xdotOne = one_density(x, ui, c(1), c(2), c(3));
    % density invariance: one density rho and the mass scaled together
    s = max(abs(xdot(1:6)));
    for kk = 1:3
        f = [1000 / 1026, 0.5, 2.0];
        xk = remus100_consistent(x, ui, c(1), c(2), c(3), ...
            struct('rho', f(kk) * wsR.rho, 'mass', f(kk) * wsR.m_given));
        invC(kk) = max(invC(kk), max(abs(xk(1:6) - xdot(1:6))) / max(s, 1e-300));
    end
    % MSS geometry and mass rule with one density: rho alone
    xm = remus100_consistent(x, ui, c(1), c(2), c(3), struct('geometry', 'mss', 'rho', 1000));
    invMss = max(invMss, max(abs(xm(1:6) - xdotOne(1:6))) / max(max(abs(xdotOne(1:6))), 1e-300));
    % differences, max |d nudot| / max |nudot| (as A-47 invariance.m)
    xMss = mssTable(k, 26:37)';
    if s > 0
        dens = max(dens, [max(abs(xdotOne(1:6) - xMss(1:6))) / max(abs(xMss(1:6))), ...
                          norm(xdotOne(1:6) - xMss(1:6)) / norm(xMss(1:6))]);
        geom = max(geom, [max(abs(xdot(1:6) - xdotOne(1:6))) / max(abs(xdotOne(1:6))), ...
                          norm(xdot(1:6) - xdotOne(1:6)) / norm(xdotOne(1:6))]);
        total = max(total, [max(abs(xdot(1:6) - xMss(1:6))) / max(abs(xMss(1:6))), ...
                            norm(xdot(1:6) - xMss(1:6)) / norm(xMss(1:6))]);
    end
    rows(k, :) = [k x' ui' c' tau' xdot' xdotMunk' xdotOne'];
end
fprintf('split: hull copy with the actuator wrench of the full copy vs the full copy, max |diff| = %.3g\n', maxSplit);
if maxSplit > 1e-12, error('the hull split is not exact (%.3g)', maxSplit); end
fprintf('density invariance, (rho, m) scaled by 1000/1026, 0.5, 2: max |d nudot| / max |nudot| = %.3g, %.3g, %.3g (%d cases)\n', invC, N);
fprintf('density invariance, MSS geometry and mass = rho * volume, rho 1000 vs 1026: %.3g\n', invMss);
% 0.5 and 2 scale every number exactly (powers of two): the derivative must
% not move at all; 1000/1026 rounds rho and m, so rounding is all that moves
if any(invC(2:3) ~= 0) || invC(1) > 1e-13 || invMss > 1e-13
    error('the copy is not density-invariant');
end
fprintf('difference, max over cases of max|d nudot|/max|nudot| and |d nudot|/|nudot|:\n');
fprintf('  density  (MSS vs one density, MSS geometry):          %.4g, %.4g\n', dens);
fprintf('  geometry (one density vs one geometry and given mass): %.4g, %.4g\n', geom);
fprintf('  total    (MSS vs consistent):                          %.4g, %.4g\n', total);
names = [{'case_id'}, vec_names('x', 12), vec_names('ui', 3), {'Vc', 'betaVc', 'w_c'}, ...
    vec_names('tau', 6), vec_names('xdot', 12), vec_names('xdot_munk', 12), vec_names('xdot_one_density', 12)];
write_csv(fullfile(outDir, 'remus100_derivative_consistent.csv'), names, rows);

% ------------------------------------------------
% 4. Stage 1: hull derivative on the cases of remus100_hull_derivative_mss.csv
% ------------------------------------------------
rng(SEED + 1, 'twister');
NH = 600;
Xh = lo + (hi - lo) .* rand(12, NH);
Th = [50; 50; 50; 10; 10; 10] .* (2 * rand(6, NH) - 1);
Ch = [1.0; 2 * pi; 0.4] .* rand(3, NH) - [0; pi; 0.2];
Ch(:, 1:2:end) = 0;
Xh(:, 1) = 0; Th(:, 1) = 0; Ch(:, 1) = 0;
hullTable = readmatrix(fullfile(scriptDir, 'remus100_hull_derivative_mss.csv'));
if max(max(abs(hullTable(:, 2:13) - Xh'))) ~= 0 || max(max(abs(hullTable(:, 14:19) - Th'))) ~= 0
    error('the hull cases differ from remus100_hull_derivative_mss.csv');
end
rows = zeros(NH, 1 + 12 + 6 + 3 + 12 + 12);
invH = 0;
for k = 1:NH
    x = Xh(:, k); t6 = Th(:, k); c = Ch(:, k);
    xd = hull(x, t6, c(1), c(2), c(3));
    f = 1000 / 1026;
    xk = remus100_consistent(x, zeros(3, 1), c(1), c(2), c(3), ...
        struct('tau_ext', f * t6, 'rho', f * wsR.rho, 'mass', f * wsR.m_given));
    if max(abs(xd(1:6))) > 0
        invH = max(invH, max(abs(xk(1:6) - xd(1:6))) / max(abs(xd(1:6))));
    end
    rows(k, :) = [k x' t6' c' xd' hull_munk(x, t6, c(1), c(2), c(3))'];
end
fprintf('density invariance, hull, (rho, m, tau_ext) scaled by 1000/1026: %.3g (%d cases)\n', invH, NH);
if invH > 1e-13, error('the hull copy is not density-invariant'); end
names = [{'case_id'}, vec_names('x', 12), vec_names('tau_ext', 6), {'Vc', 'betaVc', 'w_c'}, ...
    vec_names('xdot', 12), vec_names('xdot_munk', 12)];
write_csv(fullfile(outDir, 'remus100_hull_derivative_consistent.csv'), names, rows);

% ------------------------------------------------
% 5. Trajectories: rk4.m at h = 0.05 s (SIMremus100.m 60, 348)
% ------------------------------------------------
h = 0.05;
HOLD = 5;
OFFSET = 1e-9;
deg = pi / 180;
models = struct('full', full, 'full_munk', full_munk, 'hull', hull, 'hull_munk', hull_munk, ...
    'mss', @remus100, 'one_density', one_density);
S = struct([]);
S = add(S, 'straight_from_rest', 'full', zeros(12, 1), [0 0 0], 30, @(t, x, s) [0 0 1300]);
S = add(S, 'rudder_step_20deg', 'full', [1.5 zeros(1, 11)]', [0 0 0], 40, ...
    @(t, x, s) [20 * deg * (t >= 5) 0 1300]);
U_PRESTERO = 1.54;
surge = @(n) subsref(remus100_consistent([U_PRESTERO zeros(1, 11)]', [0 0 n]', 0, 0, 0), substruct('()', {1}));
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
% A-47's turn: rudder 15 deg at 925 rpm from 1.54 m/s, 60 s
S = add(S, 'rudder_step_15deg_at_925rpm', 'full', [1.54 zeros(1, 11)]', [0 0 0], 60, ...
    @(t, x, s) [15 * deg 0 925]);
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
% history: the same turn by unmodified MSS and by the one-density model
S = add(S, 'rudder_step_15deg_at_925rpm', 'mss', [1.54 zeros(1, 11)]', [0 0 0], 60, ...
    @(t, x, s) [15 * deg 0 925]);
S = add(S, 'rudder_step_15deg_at_925rpm', 'one_density', [1.54 zeros(1, 11)]', [0 0 0], 60, ...
    @(t, x, s) [15 * deg 0 925]);

tolRows = {};
turnEnd = struct();
for s = 1:numel(S)
    sc = S(s);
    nSteps = round(sc.T / h);
    nRows = nSteps / HOLD + 1;
    x = sc.x0; c = sc.current;
    state = struct('phase', 0, 'psi0', sc.x0(12));
    zohAll = {};
    nIn = 3; if any(strcmp(sc.kind, {'hull', 'hull_munk', 'hull_zoh'})), nIn = 6; end
    inputs = zeros(nRows, nIn); states = zeros(nRows, 12); times = (0:nRows - 1)' * HOLD * h;
    for r = 1:nRows
        t = times(r);
        states(r, :) = x';
        if r == nRows, inputs(r, :) = inputs(r - 1, :); break; end
        [u, state] = decide(sc.law, t, x, state);
        if strcmp(sc.kind, 'hull_zoh')
            % the actuator wrench of the full copy at the start of each step,
            % held over that step
            zoh = zeros(HOLD, 6);
            for k = 1:HOLD
                [~, ~, ~, tau] = remus100_consistent(x, u', c(1), c(2), c(3));
                zoh(k, :) = tau';
                x = rk4(hull, h, x, tau, c(1), c(2), c(3));
            end
            inputs(r, :) = zoh(1, :);
            zohAll{r} = zoh; %#ok<AGROW>
        else
            inputs(r, :) = u;
            for k = 1:HOLD, x = rk4(models.(sc.kind), h, x, u', c(1), c(2), c(3)); end
        end
    end
    if strcmp(sc.kind, 'hull_zoh')
        [inputs, states, times] = expand_zoh(zohAll, states, times, HOLD, h, sc, c, hull);
    end
    inNames = vec_names('ui', 3);
    if nIn == 6, inNames = vec_names('tau_ext', 6); end
    suffix = '_consistent';
    if any(strcmp(sc.kind, {'mss', 'one_density'})), suffix = ['_' sc.kind]; end
    write_csv(fullfile(outDir, ['remus100_trajectory_' sc.name suffix '.csv']), ...
        [{'t'}, inNames, {'Vc', 'betaVc', 'w_c'}, vec_names('x', 12)], ...
        [times inputs repmat(c, size(states, 1), 1) states]);
    if strcmp(sc.name, 'rudder_step_15deg_at_925rpm')
        turnEnd.(strrep(suffix(2:end), '-', '_')) = states(end, :)';
    end
    if strcmp(sc.name, 'rudder_step_4deg_at_1p54ms')
        last = times >= times(end) - 10;
        fprintf('steady yaw rate, last 10 s: %.6g deg/s at mean speed %.6g m/s\n', ...
            mean(states(last, 6)) * 180 / pi, mean(sqrt(sum(states(last, 1:3).^2, 2))));
    end
    fprintf('%-34s %-11s %5d rows, end u = %.4f m/s, r = %.5f rad/s, psi = %.4f rad, z = %.4f m\n', ...
        sc.name, sc.kind, size(states, 1), states(end, 1), states(end, 6), states(end, 12), states(end, 9));
    if any(strcmp(sc.kind, {'mss', 'one_density'})), continue, end
    % the measured difference of two integrations: (1) the derivative with
    % +-OFFSET on every entry, (2) the same RK4 with its sum reordered
    dev = zeros(3, 12);
    variants = {+OFFSET, -OFFSET, 'reordered'};
    f = models.(strrep(sc.kind, 'hull_zoh', 'hull'));
    for v = 1:3
        xa = sc.x0; dmax = zeros(1, 12);
        stride = size(states, 1) - 1; stepsPerRow = nSteps / stride;
        for r = 1:size(states, 1)
            dmax = max(dmax, abs(xa' - states(r, :)));
            if r == size(states, 1), break; end
            u = inputs(r, :);
            for k = 1:stepsPerRow
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
end
fid = fopen(fullfile(outDir, 'remus100_trajectory_tolerances_consistent.csv'), 'w');
fprintf(fid, 'scenario,state,offset_plus,offset_minus,reordered,tolerance\n');
for i = 1:size(tolRows, 1)
    fprintf(fid, '%s,%d,%.17g,%.17g,%.17g,%.17g\n', tolRows{i, :});
end
fclose(fid);
pairs = {'mss', 'one_density', 'density'; 'one_density', 'consistent', 'geometry'; 'mss', 'consistent', 'total'};
for i = 1:3
    A = turnEnd.(pairs{i, 1}); B = turnEnd.(pairs{i, 2});
    fprintf('turn 15 deg, 925 rpm, 60 s, %-8s (%s -> %s): heading %+.4f deg, turn rate %.6f vs %.6f rad/s, position %.4f m, speed %.6f vs %.6f m/s\n', ...
        pairs{i, 3}, pairs{i, 1}, pairs{i, 2}, rad2deg(B(12) - A(12)), A(6), B(6), norm(B(7:9) - A(7:9)), norm(A(1:3)), norm(B(1:3)));
end
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

function [inputs, states, times] = expand_zoh(zohAll, statesCoarse, timesCoarse, HOLD, h, sc, c, hull) %#ok<INUSL>
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
            x = rk4(hull, h, x, zohAll{i}(k, :)', c(1), c(2), c(3));
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

function write_parameters_json(file, P, rev, wsR, massGiven)
    fid = fopen(file, 'w');
    fprintf(fid, '{\n');
    fprintf(fid, '  "vehicle": "REMUS 100, the MSS 1.6 m REMUS-like model with one value per quantity",\n');
    fprintf(fid, '  "source": {"repository": "MSS (T. I. Fossen), MIT licence", "revision": "%s", "file": "CRAFT/AUV/models/remus100.m", "copy": "remus100_consistent.m"},\n', rev);
    fprintf(fid, '  "one_value": {"water_density": "%.17g kg/m^3 everywhere (remus100.m:98)", "geometry": "L = %.17g m, D = %.17g m everywhere (remus100.m:131-132)", "mass": %.17g, "mass_place": "remus100.m:3 the mass of the vehicle is 31.9 kg (given; body_density = mass / (4/3 pi a b^2))", "buoyancy": "B = W (remus100.m:214)"},\n', ...
        wsR.rho, wsR.L_auv, wsR.D_auv, massGiven);
    fprintf(fid, '  "gravity": {"latitude_rad": %.17g, "g": %.17g, "place": "remus100.m:96-97 mu = deg2rad(63.446827), g_mu = gravity(mu); weight = mass * g"},\n', wsR.mu, wsR.g_mu);
    fprintf(fid, '  "selectors": {"mass_properties": "spheroid", "coriolis": "co", "form": "submerged", "drag_model": "cylinder", "strip_grid": "midpoint", "sway_damping_fade": false, "smooth_speed_epsilon": 0.0, "open_water_model": "linearized", "open_water_coefficients": "given", "clip_advance_ratio": false, "convention": "starboard_down_positive"},\n');
    fprintf(fid, '  "parameters": {\n');
    for i = 1:size(P, 1)
        v = P{i, 2};
        if numel(v) == 1
            vs = sprintf('%.17g', v);
        else
            vs = ['[' strjoin(arrayfun(@(a) sprintf('%.17g', a), v(:)', 'UniformOutput', false), ', ') ']'];
        end
        sep = ','; if i == size(P, 1), sep = ''; end
        fprintf(fid, '    "%s": {"value": %s, "unit": "%s", "block": "%s", "block_parameter": "%s", "declared": "new", "place": "%s", "stage": "%s", "kind": "%s"}%s\n', ...
            P{i, 1}, vs, P{i, 3}, P{i, 4}, P{i, 1}, P{i, 5}, P{i, 6}, P{i, 7}, sep);
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
        fh = mss_instrumented(name, 'append');
        [varargout{:}] = fh(varargin{:});
    end
end

function fh = mss_instrumented(name, mode)
% Copy of <name>.m, written to tempdir at run time, that also returns its
% workspace as a last output. Only the function line changes and one capture
% line is added before the final 'end' (or at the end of the file when that
% 'end' closes a block); the body is the file's text.
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
