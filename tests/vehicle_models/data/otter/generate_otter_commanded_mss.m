%% ================================================================
%  Otter (MSS otter_legacy.m) commanded runs: the two propeller shaft speeds as inputs
%
%  Open-loop time histories of the unmodified CRAFT/USV/models/otter_legacy.m
%  (T. I. Fossen, MSS, MIT; revision 6e78e98) with commanded shaft speeds
%  n = [n_left n_right]' (rad/s, otter_legacy.m 32-35), a payload of 25 kg at
%  r_bP = [0.05 0 -0.35]' m (the payload of the Otter reference set) and no
%  current: rk4.m at the step of SIMotter.m (h = 0.05 s), the command decided
%  every 5 steps (0.25 s) and held over them, a row per decision. No file of
%  MSS is changed or copied: otter_legacy() is called as written, with the payload
%  trim it carries (259-269).
%  Runs: a thrust step from rest, a differential turn, an astern step from
%  speed, a spin in place, commands beyond the shaft-speed limits (otter_legacy.m
%  230: the saturation acts on both signs) and a zigzag of differential steps.
%  The two straight runs last 16 s: the Otter's open-loop yaw is unstable (the
%  Munk moment of m2c.m, otter_legacy.m 156-162), so a straight run of 40 s amplifies a
%  derivative difference of 1e-9 to 0.5 rad of heading and 2.8 m of sway, which
%  leaves no gate; at 16 s the measured sensitivity is below 1e-5 in every state.
%
%  Tolerances (as the REMUS 100 runs): per run and state the measured
%  difference of two integrations of the same run, (1) the derivative with
%  +-1e-9 added to every entry (the G1 tolerance) and (2) the same RK4 with its
%  weighted sum written in another order; the tolerance is the largest.
%
%  Outputs, under <OUT_DIR>/vehicles/otter/ (never this folder), revision tag
%  in the name:
%    otter_trajectory_<run>_mss_<revision>.csv   t, ui_01..02 (rad/s), Vc, beta_c, x_01..12
%    otter_trajectory_tolerances_mss_<revision>.csv   scenario, state, offset_plus,
%                                                offset_minus, reordered, tolerance
%    otter_commanded_meta_mss_<revision>.json    revision, step, hold, payload, shaft-speed limits
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS checkout,
%  OUT_DIR (optional) = output folder, default <tempdir>/more_mss_references;
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%  Byte-compare each output with the frozen file of the same name here.
%
%  Author:    Enio Krizman
%  Date:      2026-10-26
% ================================================================
clear functions;
format long g
format compact

mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout root.');
end
addpath(genpath(mssDir));
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outDir, 'vehicles', 'otter');
if ~exist(outDir, 'dir'), mkdir(outDir); end
[~, rev] = system(['git -C "' mssDir '" rev-parse --short HEAD']);
rev = strtrim(rev);
% MSS otter_legacy.m is the first-principles Otter (the code of otter.m up to ed2ce7b, same numbers); since 6e78e98 otter.m is the Capytaine-based model.
fprintf('MATLAB %s\nMSS: %s at %s\n', version, which('otter_legacy'), rev);

h = 0.05;               % SIMotter.m: sampling time
HOLD = 5;               % the command is decided every 5 steps (0.25 s) and held over them
OFFSET = 1e-9;          % the G1 tolerance, added to every derivative entry
mP = 25;                % payload mass (kg)
rP = [0.05 0 -0.35]';   % payload position from the CO (m)
Vc = 0; betaC = 0;      % no current
[~, ~, ~, ~, n_min, n_max] = otter_legacy();   % shaft-speed limits (rad/s), otter_legacy.m 151-152

S = struct([]);
S = add(S, 'straight_from_rest_thrust_step', zeros(12, 1), 16, ...
    @(t, x) [0 0] + [80 80] * (t >= 2));
S = add(S, 'differential_turn', [1.0 zeros(1, 11)]', 60, ...
    @(t, x) [70 70] + [30 -10] * (t >= 5));
S = add(S, 'astern_from_speed', [2.0 zeros(1, 11)]', 16, ...
    @(t, x) [70 70] * (t < 3) + [-80 -80] * (t >= 3));
S = add(S, 'spin_in_place', zeros(12, 1), 40, ...
    @(t, x) [60 -60] * (t >= 2));
S = add(S, 'commands_beyond_the_limits', zeros(12, 1), 30, ...
    @(t, x) [150 -150] * (t >= 2));
S = add(S, 'zigzag_differential_steps', [1.5 zeros(1, 11)]', 60, @zigzag);

tolRows = {};
for s = 1:numel(S)
    sc = S(s);
    nSteps = round(sc.T / h);
    nRows = nSteps / HOLD + 1;
    x = sc.x0;
    inputs = zeros(nRows, 2); states = zeros(nRows, 12); times = (0:nRows - 1)' * HOLD * h;
    for r = 1:nRows
        t = times(r);
        states(r, :) = x';
        if r == nRows, inputs(r, :) = inputs(r - 1, :); break; end
        n = sc.law(t, x);
        inputs(r, :) = n;
        for k = 1:HOLD
            x = rk4(@otter_legacy, h, x, n', mP, rP, Vc, betaC);
        end
    end
    % the measured difference of two integrations of the same run
    dev = zeros(3, 12);
    variants = {+OFFSET, -OFFSET, 'reordered'};
    for v = 1:3
        xa = sc.x0; dmax = zeros(1, 12);
        for r = 1:nRows
            dmax = max(dmax, abs(xa' - states(r, :)));
            if r == nRows, break; end
            n = inputs(r, :)';
            for k = 1:HOLD
                if ischar(variants{v})
                    xa = rk4_reordered(h, xa, n, mP, rP, Vc, betaC);
                else
                    g = @(xx, nn, a1, a2, a3, a4) otter_legacy(xx, nn, a1, a2, a3, a4) + variants{v};
                    xa = rk4(g, h, xa, n, mP, rP, Vc, betaC);
                end
            end
        end
        dev(v, :) = dmax;
    end
    for j = 1:12
        tolRows(end + 1, :) = {sc.name, j, dev(1, j), dev(2, j), dev(3, j), max(dev(:, j))}; %#ok<AGROW>
    end
    write_csv(fullfile(outDir, ['otter_trajectory_' sc.name '_mss_' rev '.csv']), ...
        [{'t'}, vec_names('ui', 2), {'Vc', 'beta_c'}, vec_names('x', 12)], ...
        [times inputs repmat([Vc betaC], nRows, 1) states]);
    fprintf('%-34s %5d rows, end u = %.4f m/s, psi = %.4f rad, z = %.4f m, max |n| = %.1f rad/s\n', ...
        sc.name, nRows, states(end, 1), states(end, 12), states(end, 9), max(abs(inputs(:))));
end
fid = fopen(fullfile(outDir, ['otter_trajectory_tolerances_mss_' rev '.csv']), 'w');
fprintf(fid, 'scenario,state,offset_plus,offset_minus,reordered,tolerance\n');
for i = 1:size(tolRows, 1)
    fprintf(fid, '%s,%d,%.17g,%.17g,%.17g,%.17g\n', tolRows{i, :});
end
fclose(fid);
fid = fopen(fullfile(outDir, ['otter_commanded_meta_mss_' rev '.json']), 'w');
fprintf(fid, '{"mss_revision": "%s", "step": %.17g, "hold_steps": %d, "payload_mass": %.17g, "payload_position": [%.17g, %.17g, %.17g], "n_min": %.17g, "n_max": %.17g}\n', ...
    rev, h, HOLD, mP, rP(1), rP(2), rP(3), n_min, n_max);
fclose(fid);
fprintf('Saved to %s\n', outDir);

%% ================================================================
% Helpers
% ================================================================
function S = add(S, name, x0, T, law)
    sc = struct('name', name, 'x0', x0, 'T', T, 'law', law);
    if isempty(S), S = sc; else, S(end + 1) = sc; end
end

function n = zigzag(t, ~)
% differential steps of alternating sign every 10 s about a mean of 70 rad/s
    sgn = (-1)^floor(t / 10);
    n = [70 70] + 25 * sgn * [1 -1];
end

function x = rk4_reordered(h, x, n, mP, rP, Vc, betaC)
% The same Runge-Kutta 4 as rk4.m, the weighted sum written in another order
    k1 = otter_legacy(x, n, mP, rP, Vc, betaC);
    k2 = otter_legacy(x + 0.5 * h * k1, n, mP, rP, Vc, betaC);
    k3 = otter_legacy(x + 0.5 * h * k2, n, mP, rP, Vc, betaC);
    k4 = otter_legacy(x + h * k3, n, mP, rP, Vc, betaC);
    x = x + h * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6);
end

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
