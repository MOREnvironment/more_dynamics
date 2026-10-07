%% ================================================================
%  Surface hydrostatics reference data from CURRENT MSS (2026-10-06)
%
%  MSS line numbers below are those of MSS ac77394 (2026-10-06), where the
%  CSVs were first made; re-run at 72656d1 (release 2.0.2, 2026-10-07): both
%  outputs byte-identical (2.0.2 moved otter.m's lines up by one).
%
%  Two CSVs beside this file, every number computed by MSS itself:
%
%  1. surface_gmtrx_mss_current.csv -- LIBRARY/modeling/Gmtrx.m on 50
%     random small-craft argument sets (nabla, A_wp, GMT, GML, x_F, r_bP),
%     r_bP off the CO, so the CF -> CO -> P transforms are covered. Gmtrx.m fixes
%     rho = 1025 and g = 9.81 (lines 25-26).
%
%  2. surface_chain_mss_current.csv -- the hydrostatic coefficient chain
%     (A_wp, I_T, I_L, KB, BM, GM, G) as MSS computes it for:
%       kind 1: CRAFT/USV/models/otter.m lines 121-194 (twin hull), 20
%               payload cases (mp, rp); case 1 is the payload of the catamaran reference
%               (mp = 25, rp = [0.05 0 -0.35]).
%       kind 2: mssExamples/exShipHydrostatics.m (monohull, script).
%       kind 3: CRAFT/SHIP/models/osv.m lines 79-102 (monohull).
%     Intermediate terms are read from the function workspace through a
%     copy instrumented at run time in tempdir (never saved), as the catamaran
%     generator test_dynamics_consistency_mss_current.m does;
%     the copy is checked against the unmodified function (xdot, M).
%
%  Numbers are written with %.17g (exact doubles).
%  Run: MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%  (MSS_DIR is required; MATLAB_BIN is the matlab executable.)
% ================================================================
clear functions;
format long g

scriptDir = fileparts(mfilename('fullpath'));
mssDir = getenv('MSS_DIR');   % no path relative to a machine
if isempty(mssDir)
    error('MSS_DIR is not set: point it at an MSS checkout (e.g. MSS_DIR=/path/to/MSS)');
end
addpath(genpath(mssDir));
fprintf('MATLAB %s\nMSS: %s\n', version, which('Gmtrx'));

%% 1. Gmtrx.m on random arguments
rng(20261006, 'twister');
N = 50;
rows = zeros(N, 8 + 36);
for k = 1:N
    % small-craft range (|G| <= ~1e6), where the frozen G1 tolerance of
    % 1e-9 absolute is reachable in double precision
    nabla = 0.01 + 1.99 * rand;
    A_wp  = 0.1 + 9.9 * rand;
    GMT   = 0.05 + 1.95 * rand;
    GML   = 0.5 + 19.5 * rand;
    x_F   = -1 + 2 * rand;
    r_bP  = -1 + 2 * rand(3, 1);
    G = Gmtrx(nabla, A_wp, GMT, GML, x_F, r_bP);
    Gt = G.';                                   % row-major, as the catamaran reference
    rows(k, :) = [nabla, A_wp, GMT, GML, x_F, r_bP', Gt(:)'];
end
names = [{'nabla','A_wp','GMT','GML','x_F','r_bP_1','r_bP_2','r_bP_3'}, ...
         arrayfun(@(i) sprintf('G_%02d', i), 1:36, 'UniformOutput', false)];
write_csv(fullfile(scriptDir, 'surface_gmtrx_mss_current.csv'), names, rows);
fprintf('Saved %d Gmtrx cases\n', N);

%% 2. Coefficient chain
chainNames = [{'kind','hull_count','L','B_hull','T','nabla','Cw','r_bg_3', ...
    'y_hull','I_L_factor','x_F','rho','g','mp','rp_1','rp_2','rp_3', ...
    'A_wp','I_T','I_L','KB','BM_T','BM_L','GM_T','GM_L'}, ...
    arrayfun(@(i) sprintf('G_%02d', i), 1:36, 'UniformOutput', false)];
chain = [];

% kind 1: otter.m, 20 payloads
otter_ws = mss_instrumented('otter');
maxSelfCheck = 0;
x = zeros(12, 1); n = [0; 0];
for k = 1:20
    if k == 1
        mp = 25; rp = [0.05; 0; -0.35];
    else
        mp = 45 * rand;                         % otter.m line 32: 0..45 kg
        rp = [-0.4 + 0.8 * rand; -0.2 + 0.4 * rand; -0.5 * rand];
    end
    [xdot0, ~, M0] = otter(x, n, mp, rp, 0, 0);
    [xdot1, ~, M1, ~, ~, ~, ws] = otter_ws(x, n, mp, rp, 0, 0);
    maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0); abs(M1(:) - M0(:))]);
    Gt = ws.G.';
    chain = [chain; 1, 2, ws.L, ws.B_pont, ws.T, ws.nabla, ws.Cw_pont, ws.rg(3), ...
        ws.y_pont, 0.8, ws.LCF, ws.rho, ws.g, mp, rp', ...
        2 * ws.Aw_pont, ws.I_T, ws.I_L, ws.KB, ws.BM_T, ws.BM_L, ws.GM_T, ws.GM_L, ...
        Gt(:)']; %#ok<AGROW>
end
% I_L_factor 0.8: otter.m line 177 (I_L = 0.8 * 2 * (1/12) * B_pont * L^3)

% kind 2: exShipHydrostatics.m (script)
s = run_script('exShipHydrostatics');
Gt = s.G.';
chain = [chain; 2, 1, s.L, s.B, s.T, s.nabla, s.Cw, s.r_bG(3), 0, 0.7, s.LCF, ...
    1025, 9.81, 0, 0, 0, 0, s.Awp, s.I_T, s.I_L, s.KB, s.BM_T, s.BM_L, s.GM_T, s.GM_L, Gt(:)'];
% rho, g: Gmtrx.m lines 25-26; I_L_factor 0.7: exShipHydrostatics.m line 49

% kind 3: osv.m (persistent vessel struct)
clear osv
osv_ws = mss_instrumented('osv');
xo = zeros(12, 1); uo = zeros(6, 1);
[xdot0, ~, M0] = osv(xo, uo);
clear osv_ws
[xdot1, ~, M1, ws] = osv_ws(xo, uo);
maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0); abs(M1(:) - M0(:))]);
v = ws.vessel;
Gt = v.G.';
chain = [chain; 3, 1, v.L, v.B, v.T, v.nabla, v.Cw, v.r_bg(3), 0, 0.7, v.LCF, ...
    1025, 9.81, 0, 0, 0, 0, v.Awp, v.I_T, v.I_L, v.KB, v.BM_T, v.BM_L, v.GM_T, v.GM_L, Gt(:)'];
% I_L_factor 0.7: osv.m line 93

write_csv(fullfile(scriptDir, 'surface_chain_mss_current.csv'), chainNames, chain);
fprintf('Instrumented vs unmodified (xdot, M), max |diff| = %.3g\n', maxSelfCheck);
fprintf('Saved %d chain cases\n', size(chain, 1));

%% ================================================================
% Helpers
% ================================================================
function write_csv(path, names, data)
% Full double precision (%.17g round-trips every double), unlike writetable.
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', strjoin(names, ','));
    fmt = [repmat('%.17g,', 1, size(data, 2) - 1) '%.17g\n'];
    fprintf(fid, fmt, data.');
    fclose(fid);
end

function s = run_script(name)
% Runs an MSS example script in this function's workspace (output captured,
% not printed) and returns its variables.
    evalc(name);
    s = struct(); vl__ = who;
    for i__ = 1:numel(vl__)
        if ~strcmp(vl__{i__}, 'name') && ~strcmp(vl__{i__}, 's')
            s.(vl__{i__}) = eval(vl__{i__});
        end
    end
end

function fh = mss_instrumented(name)
% Copy of the MSS file <name>.m, written to tempdir at run time, that also
% returns its workspace as a last output. Only the function line changes and
% one capture line is added before the final 'end'; the body is MSS's text.
% (Same helper as the catamaran and spheroid generators.)
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
