%% ================================================================
%  Payload trim of a floating craft: references from MSS otter_legacy.m
%
%  The equilibrium shift of the restoring force that otter_legacy.m gives a point
%  payload (CRAFT/USV/models/otter_legacy.m 265-275, MSS revision 6e78e98):
%
%      f_payload = Rzyx(phi,theta,psi)' * [0 0 mP*g]'        (259)
%      m_payload = Smtrx(r_bP) * f_payload                    (260)
%      g_0       = [f_payload; m_payload]                     (261-262)
%      eta_0     = [0; 0; inv(G(3:5,3:5)) * g_0(3:5); 0]      (265)
%      eta       = eta - eta_0;   J = eulerang(eta(4:6))      (266, 269)
%      xdot      = [... - G*eta ...; J*nu]                    (272-)
%
%  Two forms are written out for every case:
%    mss        : otter_legacy.m as it is. Its xdot is the unmodified function's;
%                 eta_0, g_0, the restoring vector G*eta (after the shift)
%                 and the shifted-attitude kinematics are read from a copy of
%                 otter_legacy.m instrumented at run time (in tempdir, never saved).
%    corrected  : the payload weight enters the heave equation once. The
%                 displacement nabla = (m+mP)/rho already carries the weight
%                 of hull and payload (otter_legacy.m 127-128), so g_0 keeps the
%                 lever-arm moment only, g_0 = [0 0 0 m_payload']', and the
%                 kinematics are evaluated at the actual attitude. The
%                 result is built from the instrumented workspace of the same
%                 call (G, M, C, tau, damping, cross-flow), line by line the
%                 equation of the unmodified file with those two terms changed.
%  Cases: payload masses 0, 10, 25 and 50 kg at two payload positions
%  (one on the centreline, one off it), seeded random states with small
%  attitudes, propeller speeds inside the limits of otter_legacy.m and a current.
%  The instrumented copy is checked against the unmodified otter_legacy.m (xdot) on
%  every case.
%
%  Output, under <OUT_DIR>/restoring/ (never this folder):
%    payload_trim_mss_<revision>.csv   one row per case:
%      case_id, mp, rp_1..3, x_1..12, n_1..2, Vc, beta_c,
%      eta0_1..6, g0_1..6, g_mss_1..6, G_01..36 (row-major),
%      xdot_mss_1..12, eta0c_1..6, g_corr_1..6, xdot_corr_1..12
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS checkout,
%  OUT_DIR (optional) = output folder, default <tempdir>/more_mss_references;
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%  Byte-compare the output with the frozen file of the same name here.
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
outDir = fullfile(outDir, 'restoring');
if ~exist(outDir, 'dir'), mkdir(outDir); end
[~, rev] = system(['git -C "' mssDir '" rev-parse --short HEAD']);
rev = strtrim(rev);
% MSS otter_legacy.m is the first-principles Otter (the code of otter.m up to ed2ce7b, same numbers); since 6e78e98 otter.m is the Capytaine-based model.
fprintf('MATLAB %s\nMSS: %s at %s\n', version, which('otter_legacy'), rev);

otter_ws = mss_instrumented('otter_legacy');
[~, ~, ~, ~, n_min, n_max] = otter_legacy();   % propeller speed limits (rad/s), otter_legacy.m 151-152

rng(20261026, 'twister');
masses = [0 10 25 50];
positions = [0.05 0 -0.35
            -0.30 0.12 -0.20];
perCase = 6;                            % seeded states per (mass, position)

rows = [];
maxSelfCheck = 0;
id = 0;
for mp = masses
    for ip = 1:size(positions, 1)
        rp = positions(ip, :)';
        for k = 1:perCase
            id = id + 1;
            nu  = [2.0*rand-0.5, 0.4*rand-0.2, 0.2*rand-0.1, 0.1*rand-0.05, 0.1*rand-0.05, 0.3*rand-0.15]';
            eta = [20*rand-10, 20*rand-10, 0.05*rand-0.025, ...
                   0.2*rand-0.1, 0.2*rand-0.1, 2*pi*rand-pi]';
            x = [nu; eta];
            n = n_min + (n_max - n_min) * rand(2, 1);
            Vc = 0.5 * rand;
            beta_c = 2*pi*rand - pi;

            xdot0 = otter_legacy(x, n, mp, rp, Vc, beta_c);
            [xdot1, ~, ~, ~, ~, ~, ws] = otter_ws(x, n, mp, rp, Vc, beta_c);
            maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0)]);

            g_mss = ws.G * ws.eta;                            % restoring vector of otter_legacy.m (ws.eta is shifted)

            % corrected form: lever-arm moment only in g_0, kinematics at the actual attitude
            g0c = [0; 0; 0; ws.m_payload];
            eta0c = [0; 0; ws.G(3:5, 3:5) \ g0c(3:5); 0];
            g_corr = ws.G * (x(7:12) - eta0c);
            J_actual = eulerang(x(10), x(11), x(12));
            nu_dot = ws.nu_c_dot + ws.M \ (ws.tau + ws.tau_damp + ws.tau_crossflow ...
                - ws.C * ws.nu_r - g_corr);
            xdot_corr = [nu_dot; J_actual * ws.nu];

            G_rowmajor = ws.G.'; G_rowmajor = G_rowmajor(:);
            rows(end + 1, :) = [id, mp, rp', x', n', Vc, beta_c, ws.eta_0', ws.g_0', g_mss', G_rowmajor', ...
                xdot0', eta0c', g_corr', xdot_corr']; %#ok<AGROW>
        end
    end
end
fprintf('Instrumented vs unmodified otter (xdot), max |diff| = %.3g\n', maxSelfCheck);
if maxSelfCheck ~= 0
    error('the instrumented copy differs from otter_legacy.m');
end

names = [{'case_id', 'mp'}, vec_names('rp', 3), vec_names('x', 12), vec_names('n', 2), {'Vc', 'beta_c'}, ...
    vec_names('eta0', 6), vec_names('g0', 6), vec_names('g_mss', 6), vec_names('G', 36), ...
    vec_names('xdot_mss', 12), vec_names('eta0c', 6), vec_names('g_corr', 6), vec_names('xdot_corr', 12)];
write_csv(fullfile(outDir, ['payload_trim_mss_' rev '.csv']), names, rows);
fprintf('Saved %d rows to %s\n', size(rows, 1), outDir);

%% ================================================================
% Helpers
% ================================================================
function names = vec_names(prefix, count)
    names = arrayfun(@(i) sprintf('%s_%02d', prefix, i), 1:count, 'UniformOutput', false);
end

function write_csv(path, names, rows)
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', strjoin(names, ','));
    fmt = [strjoin(repmat({'%.17g'}, 1, size(rows, 2)), ',') '\n'];
    fprintf(fid, fmt, rows.');
    fclose(fid);
end

function fh = mss_instrumented(name)
% Copy of the MSS file <name>.m, written to tempdir at run time, that also
% returns its workspace as a last output. Only the function line changes and
% one capture line is added before the final 'end'; the body is MSS's text.
    src = fileread(which(name));
    hdr = ['function\s*\[([^\]]*)\]\s*=\s*' name '\s*\('];
    if isempty(regexp(src, hdr, 'once'))
        error('%s.m: function line not found', name);
    end
    src = regexprep(src, hdr, ['function [$1,ws__] = ' name '_ws('], 'once');
    if numel(regexp(src, '(^|\n)\s*function\s')) ~= 1
        error('%s.m: expected exactly one function', name);
    end
    capture = [newline 'ws__ = struct(); vl__ = who; ' ...
        'for i__ = 1:numel(vl__), ws__.(vl__{i__}) = eval(vl__{i__}); end' newline];
    k = regexp(src, '\n[ \t]*end[ \t\r\n]*$');
    if isempty(k)
        src = [src capture];          % function without a closing 'end'
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
