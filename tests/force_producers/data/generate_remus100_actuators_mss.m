%% ================================================================
%  remus100.m actuator forces (propeller, rudder, stern planes), computed by MSS
%
%  MSS line numbers below are those of MSS ac77394 (2026-10-06), where the
%  CSV was first made; byte-identical at 72656d1 and cc07579 (release 2.0.2
%  shifted the remus100.m lines, see the data folder's SOURCE.md).
%
%  Reference data for the force-producer blocks `propeller` ("linearized")
%  and `fins` of more_dynamics:
%    - remus100.m (T. I. Fossen, MSS, MIT) is run as written; a copy
%      instrumented at run time (in tempdir, never saved) only returns its
%      workspace, so the actuator terms X_prop, K_prop, X_r, X_s, Y_r, Z_s
%      and tau (lines 164-252) can be read. No MSS line is typed here.
%    - the instrumented copy is checked against the unmodified remus100.m
%      (xdot and M) on every case.
%  Inputs (test constructions, seeded): states, rudder / stern-plane angles
%  and propeller speeds, inside and beyond the saturation limits of
%  remus100.m lines 110-116; constant current Vc = 0.3 m/s, betaVc = 30 deg,
%  w_c = 0.1 m/s so that nu_r differs from nu in all three axes.
%  Output: <OUT_DIR>/force_producers/remus100_actuators_mss_current.csv
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS checkout,
%  MATLAB_BIN = the MATLAB executable, OUT_DIR (optional) = output folder,
%  default <tempdir>/more_mss_references; never this folder, so a run cannot
%  overwrite a frozen file. Byte-compare each output with the frozen file of
%  the same name (cmp <OUT_DIR>/<block>/<file> <this data folder>/../<block>/<file>):
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%
%  Author:    Enio Krizman
%  Date:      2026-10-06
% ================================================================
clear functions;
format long g
format compact

scriptDir = fileparts(mfilename('fullpath'));
% MSS is found only through MSS_DIR (no path relative to this repository).
% Run with, e.g.,
%   MSS_DIR=/path/to/MSS "$MATLAB_BIN" -batch "run('<this file>')"
mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout root.');
end
addpath(genpath(mssDir));
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
for sub = {'force_producers'}
    if ~exist(fullfile(outDir, sub{1}), 'dir'), mkdir(fullfile(outDir, sub{1})); end
end
fprintf('MATLAB %s\nMSS: %s\n', version, which('remus100'));

remus100_ws = mss_instrumented('remus100');

% ------------------------------------------------
% Inputs (test constructions)
% ------------------------------------------------
SEED   = 20261006;
Vc     = 0.3;
beta_c = deg2rad(30);
w_c    = 0.1;

% (a) structured cases at one state: propeller speeds at and around zero
%     and the limit, fin angles at, inside and beyond the limit
x0 = [1.5 0.1 -0.05 0.02 -0.03 0.05 0 0 5 0.05 -0.1 0.3]';
n_list = [-2000 -1525 -800 -1e-3 0 1e-3 800 1525 2000];
d_list = [0 0; deg2rad(20) -deg2rad(20); 0.6 -0.6; -0.2 0.1];
X = []; UI = [];
for i = 1:numel(n_list)
    for j = 1:size(d_list, 1)
        X  = [X  x0];                                  %#ok<AGROW>
        UI = [UI [d_list(j, :) n_list(i)]'];           %#ok<AGROW>
    end
end

% (b) seeded random cases
rng(SEED, 'twister');
N_RANDOM = 1000;
lo = [-1 -1 -1 -0.5 -0.5 -0.5 -10 -10 -10 -0.5 -0.5 -pi]';
hi = [ 3  1  1  0.5  0.5  0.5  10  10  10  0.5  0.5  pi]';
X  = [X  lo + (hi - lo) .* rand(12, N_RANDOM)];
UI = [UI [-0.6; -0.6; -2000] + [1.2; 1.2; 4000] .* rand(3, N_RANDOM)];

N = size(X, 2);
allRecords = table();
maxSelfCheck = 0;

for k = 1:N
    x  = X(:, k);
    ui = UI(:, k);

    [xdot0, ~, M0] = remus100(x, ui, Vc, beta_c, w_c);
    [xdot1, ~, M1, ws] = remus100_ws(x, ui, Vc, beta_c, w_c);
    maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0); abs(M1(:) - M0(:))]);

    record = struct();
    record.case_id = k;
    for i = 1:12, record.(sprintf('x%d', i)) = x(i); end
    for i = 1:3,  record.(sprintf('ui%d', i)) = ui(i); end
    record = flatten_vector(ws.nu_r, "nu_r", record);
    record.U_r     = ws.U_r;
    record.n_p     = ws.n_p;
    record.delta_r = ws.delta_r;
    record.delta_s = ws.delta_s;
    record.X_prop  = ws.X_prop;
    record.K_prop  = ws.K_prop;
    record.X_r     = ws.X_r;
    record.X_s     = ws.X_s;
    record.Y_r     = ws.Y_r;
    record.Z_s     = ws.Z_s;
    record = flatten_vector(ws.tau, "tau", record);

    allRecords = [allRecords; struct2table(record)]; %#ok<AGROW>
end

fprintf('Instrumented vs unmodified remus100 (xdot, M), max |diff| = %.3g\n', maxSelfCheck);
fprintf('rho = %g, delta_max = %.15g rad, n_max = %g rpm, t_prop = %g, x_r = %.15g, x_s = %.15g\n', ...
    ws.rho, ws.delta_max, ws.n_max, ws.t_prop, ws.x_r, ws.x_s);

outFile = fullfile(outDir, 'force_producers', 'remus100_actuators_mss_current.csv');
writetable(allRecords, outFile);
fprintf('Saved %d rows x %d columns to %s\n', height(allRecords), width(allRecords), outFile);

%% ================================================================
% Helpers
% ================================================================
function S = flatten_vector(vec, prefix, S)
    vec = vec(:);
    for i = 1:numel(vec)
        S.(sprintf('%s_%02d', prefix, i)) = vec(i);
    end
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
