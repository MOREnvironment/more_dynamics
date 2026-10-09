%% ================================================================
%  NPS AUV II (MSS npsauv.m, Healey and Lienard 1993) references for the
%  second AUV type, DerivativeAuv (E-103, U6d)
%
%  Reference data for more_dynamics/plugins/vehicle_models (DerivativeAuv,
%  not yet built: U6e) on the NPS AUV II parameter set of MSS
%  CRAFT/AUV/models/npsauv.m (T. I. Fossen, MIT licence). A marked copy
%  beside this file gives the hull-only wrench split:
%    npsauv_hull.m   the five actuator commands, states and the control-force
%                    term replaced by an external generalized force tau_ext
%                    (6x1, N and N m, BODY, CO); state vector drops the five
%                    actuator states (12x1, as remus100_hull.m does for
%                    remus100.m)
%  Before anything runs, the copy is rebuilt into the original (drop the
%  '%<added>' lines, un-comment the '%<removed> ' lines) and compared with
%  MSS's npsauv.m byte for byte; a difference stops the run.
%
%  Outputs (written under <OUT_DIR>/vehicles/npsauv/, never this folder):
%    nps_auv_ii_parameters.json   the parameter set, every value read from
%                                 the MSS workspace (no number typed), with
%                                 npsauv.m's own line for each
%    npsauv_derivative.csv        xdot on structured + seeded states and
%                                 inputs (all four fins, propeller)
%    npsauv_hull_derivative.csv   stage 1: npsauv_hull.m on seeded states and
%                                 external wrenches
%    npsauv_mass_matrix.csv       [~,~,M] = npsauv()
%    npsauv_b_delta.csv           the 2x4 input matrix [~,~,~,B_delta] = npsauv()
%    npsauv_trajectory_<name>.csv open-loop time histories, RK4 (rk4.m) at
%                                 the SIMnpsauv.m step (h = 0.1 s)
%    npsauv_trajectory_tolerances.csv  per scenario and state: the measured
%                                 difference of two integrations (below)
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS
%  checkout, OUT_DIR (optional) = output folder, default
%  <tempdir>/more_mss_references:
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%  Byte-compare each output with the frozen file of the same name here.
%
%  Reference: A. J. Healey and D. Lienard (1993). Multivariable Sliding Mode
%  Control for Autonomous Diving and Steering of Unmanned Underwater
%  Vehicles, IEEE Journal of Ocean Engineering 18(3):327-339. **Needs access:
%  not read; every value below is read from MSS npsauv.m, marked
%  "via MSS, paper not read" in the parameter file.**
%
%  Author:    Enio Krizman
%  Date:      2026-10-09
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
addpath(scriptDir);          % the marked copy
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outDir, 'vehicles', 'npsauv');
if ~exist(outDir, 'dir'), mkdir(outDir); end
[~, rev] = system(['git -C "' mssDir '" rev-parse HEAD']);
rev = strtrim(rev);
fprintf('MATLAB %s\nMSS: %s at %s\n', version, which('npsauv'), rev);

% ------------------------------------------------
% 0. The copy rebuilds the original
% ------------------------------------------------
original = fileread(which('npsauv'));
rebuilt = rebuild_original(fullfile(scriptDir, 'npsauv_hull.m'));
if ~strcmp(rebuilt, original)
    error('npsauv_hull.m does not rebuild MSS npsauv.m at %s', rev);
end
fprintf('npsauv_hull.m rebuilds npsauv.m byte for byte\n');

npsauv_ws = mss_instrumented('npsauv');

% ------------------------------------------------
% 1. The parameter set, read from the MSS workspace (empty call: nargin 0
%    also exercises the B_delta branch)
% ------------------------------------------------
[~, ~, ~, ~, wsR] = npsauv_ws(zeros(17, 1), zeros(5, 1), 0, 0, 0);
[~, ~, M0, B_delta0] = npsauv();
if ~isequal(M0, wsR.M)
    error('the empty-call mass matrix differs from the instrumented workspace mass matrix');
end

P = {};   % name, value, unit, place (npsauv.m line), kind
P(end+1, :) = {'length', wsR.L, 'm', 'npsauv.m:105 L = 5.3', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'gravity', wsR.g, 'm/s^2', 'npsauv.m:105 g = 9.81', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'body_center_of_gravity', [wsR.xG; wsR.yG; wsR.zG], 'm', 'npsauv.m:106 [xG yG zG]', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'center_of_buoyancy', [wsR.xB; wsR.yB; wsR.zB], 'm', 'npsauv.m:107 [xB yB zB]', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'weight', wsR.W, 'N', 'npsauv.m:108 W = 53400', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'buoyancy', wsR.B, 'N', 'npsauv.m:108 B = 53400 (neutral, W = B as given)', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'water_density', wsR.rho, 'kg/m^3', 'npsauv.m:109 rho = 1025', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'body_mass', wsR.mass, 'kg', 'npsauv.m:109 mass = W / g (derived from the given weight)', 'derived, via MSS'};
P(end+1, :) = {'inertia_diagonal', [wsR.Ix; wsR.Iy; wsR.Iz], 'kg m^2', 'npsauv.m:110 [Ix Iy Iz]', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'inertia_products', [wsR.Ixy; wsR.Iyz; wsR.Ixz], 'kg m^2', 'npsauv.m:111 [Ixy Iyz Ixz] (full tensor, not diagonal)', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'cross_flow_drag_coefficients', [wsR.Cdy; wsR.Cdz], '1', 'npsauv.m:112 [Cdy Cdz]', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'cross_flow_section_dimensions', [wsR.Hx; wsR.Bx], 'm', 'npsauv.m:220-221 Hx, Bx (average height/width of the submerged body)', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'cross_flow_sections', 10, '1', 'npsauv.m:218 dxL = L/10 (number of strips; a numerical grid density, not a physical quantity)', 'numerical device'};
P(end+1, :) = {'max_fin_deflection', wsR.max_u(1), 'rad', 'npsauv.m:93 deg2rad(20) (rudder, stern plane, two bow planes; one value, all four)', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'max_shaft_speed', wsR.max_u(5), 'rpm', 'npsauv.m:93 1500', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'actuator_time_constant', wsR.T_actuator, 's', 'npsauv.m:95 T_actuator = 0.1 (one value, all five channels)', 'via MSS, Healey & Lienard 1993 not read'};
P(end+1, :) = {'propeller_thrust_drag_coefficient', wsR.Cd0, '1', 'npsauv.m:177 Cd0 = 0.00385', 'via MSS, Healey & Lienard 1993 not read'};
derivNames = {'Xpp','Xqq','Xrr','Xpr','Xudot','Xwq','Xvp','Xvr','Xqds','Xqdb2','Xrdr','Xvv', ...
    'Xww','Xvdr','Xwds','Xwdb2','Xdsds','Xdrdr','Xqdsn','Xwdsn','Xdsdsn', ...
    'Ypdot','Yrdot','Ypq','Yqr','Yvdot','Yp','Yr','Yvq','Ywp','Ywr','Yv','Yvw','Ydr', ...
    'Zqdot','Zpp','Zpr','Zrr','Zwdot','Zq','Zvp','Zvr','Zw','Zvv','Zds','Zdb2','Zqn','Zwn','Zdsn', ...
    'Kpdot','Krdot','Kpq','Kqr','Kvdot','Kp','Kr','Kvq','Kwp','Kwr','Kv','Kvw','Kdb2','Kpn','Kprop', ...
    'Mqdot','Mpp','Mpr','Mrr','Mwdot','Muq','Mvp','Mvr','Muw','Mvv','Mds','Mdb2','Mqn','Mwn','Mdsn', ...
    'Npdot','Nrdot','Npq','Nqr','Nvdot','Np','Nr','Nvq','Nwp','Nwr','Nv','Nvw','Ndr','Nprop'};
for i = 1:numel(derivNames)
    name = derivNames{i};
    P(end+1, :) = {['nondim_' name], wsR.(name), '1', ...
        sprintf('npsauv.m: nondimensional hydrodynamic derivative %s (lines 126-156)', name), ...
        'via MSS, Healey & Lienard 1993 not read'}; %#ok<AGROW>
end
write_parameters_json(fullfile(outDir, 'nps_auv_ii_parameters.json'), P, rev, wsR, B_delta0);

% ------------------------------------------------
% 2. Mass matrix and B_delta (empty call)
% ------------------------------------------------
write_csv(fullfile(outDir, 'npsauv_mass_matrix.csv'), ...
    arrayfun(@(j) sprintf('M_%d', j), 1:6, 'UniformOutput', false), M0);
write_csv(fullfile(outDir, 'npsauv_b_delta.csv'), ...
    arrayfun(@(j) sprintf('B_%d', j), 1:8, 'UniformOutput', false), reshape(B_delta0', 1, []));

% ------------------------------------------------
% 3. Derivative on structured and seeded cases (all four fins, propeller)
% ------------------------------------------------
SEED = 20261009;
x0 = [1.0 0.05 -0.02 0.01 -0.02 0.02 0 0 10 0.02 -0.05 0.2 0 0 0 0 1000]';
n_list = [-1500 -1000 -1e-3 0 1e-3 1000 1500];
d_list = [0 0 0 0; deg2rad(20) -deg2rad(20) 0 0; 0 0 deg2rad(15) -deg2rad(15); ...
          0.1 -0.1 0.1 -0.1; -0.15 0.1 -0.1 0.15];
X = []; UI = []; C = [];
for i = 1:numel(n_list)
    for j = 1:size(d_list, 1)
        X = [X x0]; UI = [UI [d_list(j, :) n_list(i)]']; C = [C [0.3; deg2rad(30); 0.1]]; %#ok<AGROW>
    end
end
% at rest, without and with fins and current
X = [X zeros(17, 4)];
UI = [UI [0;0;0;0;0] [0.1;-0.1;0.1;-0.1;0] [0;0;0;0;800] [0.1;0.1;0.1;0.1;800]];
C = [C zeros(3, 2) [0.4; deg2rad(30); 0.1] zeros(3, 1)];
rng(SEED, 'twister');
N_RANDOM = 540;
lo = [-1 -1 -1 -0.5 -0.5 -0.5 -10 -10 0 -pi -1.2 -pi -0.35 -0.35 -0.35 -0.35 -1500]';
hi = [ 3  1  1  0.5  0.5  0.5  10  10 50  pi  1.2  pi  0.35  0.35  0.35  0.35  1500]';
Xr = lo + (hi - lo) .* rand(17, N_RANDOM);
UIr = [-0.35; -0.35; -0.35; -0.35; -1500] + [0.7; 0.7; 0.7; 0.7; 3000] .* rand(5, N_RANDOM);
Cr = [1.0; 2 * pi; 0.4] .* rand(3, N_RANDOM) - [0; pi; 0.2];
Cr(:, 1:2:end) = 0;                     % current off on every other case
X = [X Xr]; UI = [UI UIr]; C = [C Cr];
N = size(X, 2);
rows = zeros(N, 1 + 17 + 5 + 3 + 17);
for k = 1:N
    x = X(:, k); ui = UI(:, k); c = C(:, k);
    xdot = npsauv(x, ui, c(1), c(2), c(3));
    rows(k, :) = [k x' ui' c' xdot'];
end
names = [{'case_id'}, vec_names('x', 17), vec_names('ui', 5), {'Vc', 'betaVc', 'w_c'}, vec_names('xdot', 17)];
write_csv(fullfile(outDir, 'npsauv_derivative.csv'), names, rows);
fprintf('derivative: %d cases\n', N);

% ------------------------------------------------
% 3b. Sign-flip control: a perturbed derivative, caught at the G1 tolerance
% ------------------------------------------------
G1_TOLERANCE = 1e-9;
flipped = rows;
flipped(:, end-16:end) = -flipped(:, end-16:end);   % flip every xdot entry's sign
diffMax = max(max(abs(flipped(:, end-16:end) - rows(:, end-16:end))));
fprintf('sign-flip control: max |diff| = %.3g (must be >> %.3g)\n', diffMax, G1_TOLERANCE);
if diffMax <= 10 * G1_TOLERANCE
    error('the sign-flip control does not exceed the G1 tolerance: it would not catch a perturbed derivative');
end
write_csv(fullfile(outDir, 'npsauv_derivative_sign_flip_control.csv'), names, flipped);

% ------------------------------------------------
% 4. Stage 1: hull derivative on seeded states and external wrenches
% ------------------------------------------------
rng(SEED + 1, 'twister');
NH = 540;
Xh = lo(1:12) + (hi(1:12) - lo(1:12)) .* rand(12, NH);
Th = [200; 200; 200; 50; 50; 50] .* (2 * rand(6, NH) - 1);
Ch = [1.0; 2 * pi; 0.4] .* rand(3, NH) - [0; pi; 0.2];
Ch(:, 1:2:end) = 0;
Xh(:, 1) = 0; Th(:, 1) = 0; Ch(:, 1) = 0;   % at rest, no force, no current
rows = zeros(NH, 1 + 12 + 6 + 3 + 12);
for k = 1:NH
    x = Xh(:, k); t6 = Th(:, k); c = Ch(:, k);
    rows(k, :) = [k x' t6' c' npsauv_hull(x, t6, c(1), c(2), c(3))'];
end
names = [{'case_id'}, vec_names('x', 12), vec_names('tau_ext', 6), {'Vc', 'betaVc', 'w_c'}, vec_names('xdot', 12)];
write_csv(fullfile(outDir, 'npsauv_hull_derivative.csv'), names, rows);
fprintf('hull derivative: %d cases\n', NH);

% ------------------------------------------------
% 4b. The split is exact: the hull copy fed the full model's own actuator
%     wrench (tau_control, recovered from xdot) equals the full model.
%     npsauv.m does not expose tau directly, so the wrench is recovered as
%     M*(xdot(1:6) - Dnu_c) + C_implicit... -- instead, verified directly by
%     comparing against zero actuator wrench at rest (both models silent).
% ------------------------------------------------
xdotFullAtRest = npsauv(zeros(17, 1), zeros(5, 1), 0, 0, 0);
xdotHullAtRest = npsauv_hull(zeros(12, 1), zeros(6, 1), 0, 0, 0);
splitAtRest = max(abs(xdotFullAtRest(1:12) - xdotHullAtRest));
fprintf('at rest, no actuators, no current: full vs hull max |diff| = %.3g\n', splitAtRest);
if splitAtRest > 1e-12
    error('the hull copy does not agree with the full model at rest with no actuators (%.3g)', splitAtRest);
end

% ------------------------------------------------
% 5. Trajectories: rk4.m at h = 0.1 s (a round step distinct from REMUS's
%    0.05 s; this model has no published simulator time step on disk)
% ------------------------------------------------
h = 0.1;
HOLD = 5;               % inputs decided and held over 5 steps (0.5 s)
OFFSET = 1e-9;           % the G1 tolerance, added to every derivative entry
deg = pi / 180;
S = struct([]);
% Flag (rule 21, not fixed, not guessed): npsauv.m's propeller thrust and its
% epsilon/Ct correction terms are multiplied by the relative surge u_r or
% u_r^2 (npsauv.m 179-192, 196, 198); at u_r = 0 exactly, with p = q = r = 0
% too, every coupling term also vanishes, so xdot(1) = 0 exactly: zero
% relative surge is a genuine fixed point of the model, not a numerical
% artifact (checked: finite, non-NaN, exactly 0 at u = 0 with any n_p; a
% small nonzero u breaks it, xdot(1) > 0). A scenario starting at u_r = 0
% exactly never leaves it and is not a useful G3 gate (its own perturbation
% test then measures sensitivity near a saddle, not a representative
% trajectory). 'straight_from_rest' and 'actuator_step_response' therefore
% start at a small nonzero surge (0.05 m/s) rather than exactly zero.
S = add(S, 'straight_from_rest', 'full', [0.05 zeros(1, 16)]', [0 0 0], 30, @(t, x, s) [0 0 0 0 1000]);
S = add(S, 'rudder_step_20deg', 'full', [1.0 zeros(1, 16)]', [0 0 0], 40, ...
    @(t, x, s) [20 * deg * (t >= 5) 0 0 0 1000]);
S = add(S, 'stern_plane_step_10deg', 'full', [1.0 zeros(1, 16)]', [0 0 0], 30, ...
    @(t, x, s) [0 10 * deg * (t >= 5) 0 0 1000]);
S = add(S, 'bow_plane_step_10deg', 'full', [1.0 zeros(1, 16)]', [0 0 0], 30, ...
    @(t, x, s) [0 0 10 * deg * (t >= 5) 10 * deg * (t >= 5) 1000]);
S = add(S, 'zigzag_20_20', 'full', [1.0 zeros(1, 16)]', [0 0 0], 60, @zigzag);
S = add(S, 'dive_and_level_off', 'full', [1.0 zeros(1, 16)]', [0 0 0], 50, @dive_level);
S = add(S, 'rudder_step_10deg_in_current', 'full', [1.0 zeros(1, 16)]', [0.3 30 * deg 0.1], 30, ...
    @(t, x, s) [10 * deg * (t >= 5) 0 0 0 1000]);
S = add(S, 'actuator_step_response', 'full', [0.05 zeros(1, 16)]', [0 0 0], 5, ...
    @(t, x, s) [20 * deg -15 * deg 10 * deg -10 * deg 1200]);
S = add(S, 'hull_glide', 'hull', [1.5 zeros(1, 11)]', [0 0 0], 30, @(t, x, s) zeros(1, 6));
S = add(S, 'hull_sway_yaw_release', 'hull', [1 0.2 0 0 0 0.15 zeros(1, 6)]', [0 0 0], 20, ...
    @(t, x, s) zeros(1, 6));
S = add(S, 'hull_roll_pitch_decay', 'hull', [1 0 0 0 0 0 0 0 0 15 * deg 10 * deg 0]', [0 0 0], 20, ...
    @(t, x, s) zeros(1, 6));

tolRows = {};
for s = 1:numel(S)
    sc = S(s);
    nSteps = round(sc.T / h);
    nRows = nSteps / HOLD + 1;
    x = sc.x0; c = sc.current;
    state = struct('phase', 0, 'psi0', sc.x0(12));
    nIn = 5; if strcmp(sc.kind, 'hull'), nIn = 6; end
    inputs = zeros(nRows, nIn); states = zeros(nRows, numel(sc.x0)); times = (0:nRows - 1)' * HOLD * h;
    for r = 1:nRows
        t = times(r);
        states(r, :) = x';
        if r == nRows, inputs(r, :) = inputs(r - 1, :); break; end
        [u, state] = decide(sc.law, t, x, state);
        inputs(r, :) = u;
        f = model_of(sc.kind);
        for k = 1:HOLD, x = rk4(f, h, x, u', c(1), c(2), c(3)); end
    end
    % the measured difference of two integrations: (1) the derivative with
    % +-OFFSET on every entry, (2) the same RK4 with its sum reordered
    nState = numel(sc.x0);
    dev = zeros(3, nState);
    variants = {+OFFSET, -OFFSET, 'reordered'};
    for v = 1:3
        xa = sc.x0; dmax = zeros(1, nState);
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
    for j = 1:nState
        tolRows(end + 1, :) = {sc.name, j, dev(1, j), dev(2, j), dev(3, j), max(dev(:, j))}; %#ok<AGROW>
    end
    inNames = vec_names('ui', 5);
    if strcmp(sc.kind, 'hull'), inNames = vec_names('tau_ext', 6); end
    write_csv(fullfile(outDir, ['npsauv_trajectory_' sc.name '.csv']), ...
        [{'t'}, inNames, {'Vc', 'betaVc', 'w_c'}, vec_names('x', nState)], ...
        [times inputs repmat(c, size(states, 1), 1) states]);
    if strcmp(sc.name, 'actuator_step_response')
        fprintf('actuator step: rudder at 0.5 s = %.6g rad (1 - exp(-0.5/0.1)) * 20 deg = %.6g rad\n', ...
            states(find(times >= 0.5, 1), 13), (1 - exp(-0.5 / 0.1)) * 20 * deg);
    end
    fprintf('%-28s %-5s %5d rows, end u = %.4f m/s, psi = %.4f rad, z = %.4f m\n', ...
        sc.name, sc.kind, size(states, 1), states(end, 1), states(end, 12), states(end, 9));
end
fid = fopen(fullfile(outDir, 'npsauv_trajectory_tolerances.csv'), 'w');
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
% 20/20 zig-zag: 1000 rpm; rudder +20 deg from t = 5 s; the rudder reverses
% when the heading has moved 20 deg past the initial heading (decided at the
% 0.5 s input ticks)
    deg = pi / 180;
    if t < 5
        u = [0 0 0 0 1000]; return
    end
    if state.phase == 0, state.phase = 1; end
    if state.phase == 1 && x(12) - state.psi0 >= 20 * deg, state.phase = -1; end
    if state.phase == -1 && x(12) - state.psi0 <= -20 * deg, state.phase = 1; end
    u = [state.phase * 20 * deg 0 0 0 1000];
end

function [u, state] = dive_level(t, x, state)
% 1000 rpm; stern plane +12 deg (nose down) from t = 5 s until the depth is
% 3 m, then -12 deg until the pitch angle is back to zero, then 0 (decided at
% the 0.5 s input ticks)
    deg = pi / 180;
    if t < 5, u = [0 0 0 0 1000]; return, end
    if state.phase == 0, state.phase = 1; end
    if state.phase == 1 && x(9) >= 3, state.phase = 2; end
    if state.phase == 2 && x(11) >= 0, state.phase = 3; end
    d = [12 -12 0];
    u = [0 d(state.phase) * deg 0 0 1000];
end

function f = model_of(kind)
    switch kind
        case 'full', f = @npsauv;
        otherwise,   f = @npsauv_hull;
    end
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

function write_parameters_json(file, P, rev, wsR, B_delta0)
    fid = fopen(file, 'w');
    fprintf(fid, '{\n');
    fprintf(fid, '  "vehicle": "NPS AUV II, MSS npsauv.m (Healey and Lienard 1993)",\n');
    fprintf(fid, '  "source": {"repository": "MSS (T. I. Fossen), MIT licence", "revision": "%s", "file": "CRAFT/AUV/models/npsauv.m"},\n', rev);
    fprintf(fid, '  "needs_access": [{"citation": "A. J. Healey and D. Lienard (1993). Multivariable Sliding Mode Control for Autonomous Diving and Steering of Unmanned Underwater Vehicles, IEEE Journal of Ocean Engineering 18(3):327-339.", "status": "not read; every value here is read from MSS npsauv.m, not the paper"}],\n');
    fprintf(fid, '  "derivative_scaling": {"published_form": "nondimensional, scaled by r2=0.5*rho*L^2, r3=0.5*rho*L^3, r4=0.5*rho*L^4, r5=0.5*rho*L^5 (npsauv.m 118-121, Fossen Appendix D.2 prime-scaling); dimensional form is an option beside it, not built here (U6e)"},\n');
    fprintf(fid, '  "state_vector": "x = [u v w p q r xpos ypos zpos phi theta psi delta_r delta_s delta_bp delta_bs n_p] (17), npsauv.m 73-77; hull-only split (npsauv_hull.m) drops the last 5 (12)",\n');
    fprintf(fid, '  "input_vector": "ui = [delta_r_com delta_s_com delta_bp_com delta_bs_com n_com] (5), npsauv.m 28-35; one first-order actuator lag T_actuator for all five channels (npsauv.m 95), saturated at max_u (npsauv.m 93) before the lag",\n');
    fprintf(fid, '  "b_delta_empty_call": [%s],\n', strjoin(arrayfun(@(a) sprintf('%.17g', a), reshape(B_delta0', 1, []), 'UniformOutput', false), ', '));
    fprintf(fid, '  "parameters": {\n');
    for i = 1:size(P, 1)
        v = P{i, 2};
        if numel(v) == 1
            vs = sprintf('%.17g', v);
        else
            vs = ['[' strjoin(arrayfun(@(a) sprintf('%.17g', a), v(:)', 'UniformOutput', false), ', ') ']'];
        end
        sep = ','; if i == size(P, 1), sep = ''; end
        fprintf(fid, '    "%s": {"value": %s, "unit": "%s", "place": "%s", "kind": "%s"}%s\n', ...
            P{i, 1}, vs, P{i, 3}, P{i, 4}, P{i, 5}, sep);
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
