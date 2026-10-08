%% ================================================================
%  remus100.m rudder and stern planes, computed by MSS, propeller at rest
%
%  Reference data for the fin parts and the fin compositions of
%  more_dynamics (servo "ideal", inflow "translational", flow angle "none",
%  section "quadratic_drag", and the deflection-only composition):
%    - remus100.m (T. I. Fossen, MSS, MIT) is run as written; a copy
%      instrumented at run time (in tempdir, never saved) only returns its
%      workspace, so the intermediate quantities of the fin lines can be
%      read: the saturated angles delta_r, delta_s (lines 113-114), the
%      section-plane speeds U_rh, U_rv (234-235), the forces X_r, X_s, Y_r,
%      Z_s (238-245), the constants rho, A_r, A_s, CL_delta_r, CL_delta_s,
%      x_r, x_s, delta_max (98, 109, 179-189) and tau (248-254).
%      No MSS line is typed here.
%    - the propeller command ui(3) is 0 in every case, so remus100.m's own
%      X_prop and K_prop are 0 (lines 164-176) and its tau is the fins'
%      share alone; both are written so a reader can check it.
%    - the instrumented copy is checked against the unmodified remus100.m
%      (xdot and M) on every case.
%  Inputs (test constructions, seeded): states and rudder / stern-plane
%  angles inside and beyond the limit of line 109; structured cases with no
%  current and zero, forward, astern, sideways and vertical speed; random
%  cases with a constant current Vc = 0.3 m/s, betaVc = 30 deg, w_c = 0.1 m/s
%  so that nu_r differs from nu in all three axes.
%  MSS line numbers are those of MSS cc07579 (2026-10-07).
%  Output: <OUT_DIR>/force_producers/fin_parts/remus100_fins_mss_cc07579.csv
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS checkout,
%  MATLAB_BIN = the MATLAB executable, OUT_DIR (optional) = output folder,
%  default <tempdir>/more_mss_references; never this folder, so a run cannot
%  overwrite a frozen file. Byte-compare the output with the frozen file:
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%    cmp <OUT_DIR>/force_producers/fin_parts/remus100_fins_mss_cc07579.csv <this folder>/remus100_fins_mss_cc07579.csv
%
%  Author:    Enio Krizman
%  Date:      2026-10-08
% ================================================================
clear functions;
format long g
format compact

% MSS is found only through MSS_DIR (no path relative to this repository).
mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout root.');
end
addpath(genpath(mssDir));
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
target = fullfile(outDir, 'force_producers', 'fin_parts');
if ~exist(target, 'dir'), mkdir(target); end
fprintf('MATLAB %s\nMSS: %s\n', version, which('remus100'));

remus100_ws = mss_instrumented('remus100');

% ------------------------------------------------
% Inputs (test constructions)
% ------------------------------------------------
SEED = 20261008;

% (a) structured cases, no current (Vc = 0, w_c = 0): nu = nu_r
nu_list = [ 0    0    0    0    0    0;      % at rest: U_rh = U_rv = 0
            1.5  0    0    0    0    0;      % straight ahead
           -1.0  0    0    0    0    0;      % straight astern
            0    0.4  0    0    0    0;      % sideways only (u = 0)
            0    0    0.3  0    0    0;      % vertical only (u = 0)
            1.2  0.2 -0.1  0.05 -0.1 0.2];   % general
d_list  = [ 0         0;
            deg2rad(20) -deg2rad(20);        % at the limit
            0.6       -0.6;                  % beyond the limit
           -0.2        0.1];
X = []; UI = []; CUR = [];
for i = 1:size(nu_list, 1)
    for j = 1:size(d_list, 1)
        x = zeros(12, 1);
        x(1:6) = nu_list(i, :)';
        X   = [X x];                                   %#ok<AGROW>
        UI  = [UI [d_list(j, :) 0]'];                  %#ok<AGROW>
        CUR = [CUR [0; 0; 0]];                         %#ok<AGROW>
    end
end

% (b) seeded random cases with a current
rng(SEED, 'twister');
N_RANDOM = 1000;
lo = [-1 -1 -1 -0.5 -0.5 -0.5 -10 -10 -10 -0.5 -0.5 -pi]';
hi = [ 3  1  1  0.5  0.5  0.5  10  10  10  0.5  0.5  pi]';
X   = [X  lo + (hi - lo) .* rand(12, N_RANDOM)];
UI  = [UI [[-0.6; -0.6] + [1.2; 1.2] .* rand(2, N_RANDOM); zeros(1, N_RANDOM)]];
CUR = [CUR repmat([0.3; deg2rad(30); 0.1], 1, N_RANDOM)];

N = size(X, 2);
allRecords = table();
maxSelfCheck = 0;

for k = 1:N
    x  = X(:, k);
    ui = UI(:, k);
    c  = CUR(:, k);

    [xdot0, ~, M0] = remus100(x, ui, c(1), c(2), c(3));
    [xdot1, ~, M1, ws] = remus100_ws(x, ui, c(1), c(2), c(3));
    maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0); abs(M1(:) - M0(:))]);

    record = struct();
    record.case_id = k;
    for i = 1:12, record.(sprintf('x%d', i)) = x(i); end
    for i = 1:3,  record.(sprintf('ui%d', i)) = ui(i); end
    record.Vc     = c(1);
    record.betaVc = c(2);
    record.w_c    = c(3);
    record = flatten_vector(ws.nu_r, "nu_r", record);
    record.delta_r = ws.delta_r;
    record.delta_s = ws.delta_s;
    record.U_rh    = ws.U_rh;
    record.U_rv    = ws.U_rv;
    record.X_r     = ws.X_r;
    record.X_s     = ws.X_s;
    record.Y_r     = ws.Y_r;
    record.Z_s     = ws.Z_s;
    record.X_prop  = ws.X_prop;
    record.K_prop  = ws.K_prop;
    record = flatten_vector(ws.tau, "tau", record);
    record.rho        = ws.rho;
    record.delta_max  = ws.delta_max;
    record.A_r        = ws.A_r;
    record.A_s        = ws.A_s;
    record.CL_delta_r = ws.CL_delta_r;
    record.CL_delta_s = ws.CL_delta_s;
    record.x_r        = ws.x_r;
    record.x_s        = ws.x_s;

    allRecords = [allRecords; struct2table(record)]; %#ok<AGROW>
end

fprintf('Instrumented vs unmodified remus100 (xdot, M), max |diff| = %.3g\n', maxSelfCheck);
fprintf('max |X_prop| = %.3g, max |K_prop| = %.3g (propeller at rest)\n', ...
    max(abs(allRecords.X_prop)), max(abs(allRecords.K_prop)));

outFile = fullfile(target, 'remus100_fins_mss_cc07579.csv');
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
