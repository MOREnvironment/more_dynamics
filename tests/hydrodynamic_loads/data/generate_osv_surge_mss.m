%% ================================================================
%  Surge damping of MSS osv.m at release 2.0.2 (72656d1) and after the
%  fix of the persistent surge coefficient (cc07579), 2026-10-07
%
%  MSS 2.0.2 removed LIBRARY/modeling/forceSurgeDamping.m; CRAFT/SHIP/models/
%  osv.m now damps surge with (osv.m 178-185 @ 72656d1)
%
%    D_nl_11 = exp(-k_u*|u_r|) * D11 - Xuu * |u_r|,   k_u = 3,
%    X       = -D_nl_11 * u_r     (surge row of -(CRB + CA + D_nonlinear)*nu_r,
%                                  osv.m 203-207; D from Dmtrx.m is diagonal)
%
%  with D11 = M(1,1)/T1 (LIBRARY/modeling/Dmtrx.m 62) and Xuu from
%  LIBRARY/modeling/XuuITTC.m (ITTC-1957 line, Rn bounded at 1e5, Mumford
%  wetted area). The same formula is evaluated here on the speed grid of
%  surge_damping_mss_ac77394.csv (sets 1-5, branch 0, all 44 speeds), with
%  B, T, C_B of each set from xuu_ittc_mss_ac77394.csv:
%    set 1   = osv.m itself: D11 = vessel.D(1,1) and the call of the
%              unmodified osv.m are read from the function (instrumented
%              copy in tempdir, checked equal to osv.m);
%    sets 2-5 = D11 = Dmtrx(...)(1,1) with M(1,1) = m + A11 (A11 of
%              addedMassSurge.m, the column of the grid file).
%
%  At 72656d1 osv.m keeps its parameters in a persistent struct and sets
%  vessel.D(1,1) = 0 (line 188) after line 185 has read it, so the linear
%  term is present only on the first call after 'clear osv'; MSS commit
%  de9a316 (2026-10-07) removed that line, so from cc07579 on both calls
%  give the same value. For set 1 the file holds the first and the second
%  call of osv.m on the same state (columns *_first, *_second); NaN for
%  sets 2-5. The script detects which osv.m it runs and names the output.
%
%  Numbers are written with %.17g (exact doubles).
%  Run (MSS_DIR required, MATLAB_BIN = the MATLAB executable; OUT_DIR
%  optional, default <tempdir>/more_mss_references; files land in its hydrodynamics/ subfolder, never this folder, so a
%  run cannot overwrite a frozen file):
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%  Output: osv_surge_damping_mss_2_0_2.csv (osv.m zeroes the coefficient,
%  72656d1) or osv_surge_damping_mss_cc07579.csv (it does not, cc07579 and
%  later); byte-compare it with the frozen file of the same name.
%
%  Author:    Enio Krizman
%  Date:      2026-10-07
% ================================================================
clear functions;
format long g

scriptDir = fileparts(mfilename('fullpath'));
mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: point it at an MSS checkout (export MSS_DIR=<path>)');
end
addpath(genpath(mssDir));
if exist('forceSurgeDamping', 'file')
    error('forceSurgeDamping.m is on the path: this file needs MSS 2.0.2 or later');
end
fprintf('MATLAB %s\nMSS: %s\n', version, which('osv'));

k_u = 3;                                       % osv.m 183

sg = readmatrix(fullfile(scriptDir, 'surge_damping_mss_ac77394.csv'));
sg = sg(sg(:, 2) == 0, :);               % branch 0: one row per speed
ittc = readmatrix(fullfile(scriptDir, 'xuu_ittc_mss_ac77394.csv'));
% grid columns: set,branch,probe,m,S,L,T1,rho,u_max,thrust_max,u_r,X,Xuu,Xu,A11
% ittc columns: set,rho,L,B,T,C_B,u_r,Xuu

osvW = mss_instrumented('osv');
clear osv_ws
[~, ~, ~, ws0] = osvW(zeros(12, 1), zeros(6, 1));
% vessel.D(1,1) is already zeroed in the captured workspace (line 188), so
% D11 is recomputed with osv.m's own call (lines 133-134)
v = ws0.vessel;
Dosv = Dmtrx([v.T1, v.T2, v.T6], [v.zeta4, v.zeta5], v.MRB, v.MA, v.G);
D11_osv = Dosv(1, 1);
if v.D(1, 1) == 0
    tag = '2_0_2';        % osv.m zeroes the persistent coefficient (72656d1)
elseif v.D(1, 1) == D11_osv
    tag = 'cc07579';      % the zeroing line is gone (de9a316, cc07579)
else
    error('vessel.D(1,1) after a call is neither 0 nor Dmtrx(1,1): revisit osv.m');
end
fprintf('osv.m variant: %s\n', tag);

hdr = {'set', 'probe', 'm', 'L', 'B', 'T', 'C_B', 'rho', 'T1', 'M11', 'D11', ...
       'k_u', 'u_r', 'Xuu', 'D_nl_11', 'X', ...
       'D_nl_11_osv_first', 'D_nl_11_osv_second', ...
       'nu_dot_1_osv_first', 'nu_dot_1_osv_second'};
rows = zeros(size(sg, 1), numel(hdr));
maxSelf = 0; maxFirst = 0;
for k = 1:size(sg, 1)
    s = sg(k, 1); probe = sg(k, 3); m = sg(k, 4); L = sg(k, 6);
    T1 = sg(k, 7); rho = sg(k, 8); u = sg(k, 11); A11 = sg(k, 15);
    p = ittc(find(ittc(:, 1) == s, 1), :);
    B = p(4); T = p(5); C_B = p(6);
    if p(2) ~= rho || p(3) ~= L
        error('set %d: rho or L differ between the two grid files', s);
    end
    if s == 1
        M11 = ws0.vessel.M(1, 1);
        D11 = D11_osv;
    else
        M11 = m + A11;
        Dfull = Dmtrx([T1 T1 T1], [0.15 0.3], m * eye(6), A11 * eye(6), eye(6));
        D11 = Dfull(1, 1);
    end
    Xuu = XuuITTC(u, rho, L, B, T, C_B);
    D_nl_11 = exp(-k_u * abs(u)) * D11 - Xuu * abs(u);   % osv.m 185
    X = -D_nl_11 * u;
    first = NaN; second = NaN; ndf = NaN; nds = NaN;
    if s == 1
        x = [u; zeros(11, 1)];
        ui = zeros(6, 1);
        clear osv osv_ws
        xd1 = osv(x, ui);  xd2 = osv(x, ui);          % unmodified osv.m
        [xw1, ~, ~, w1] = osvW(x, ui); [xw2, ~, ~, w2] = osvW(x, ui);
        maxSelf = max([maxSelf; abs(xw1 - xd1); abs(xw2 - xd2)]);
        first = w1.D_nonlinear(1, 1); second = w2.D_nonlinear(1, 1);
        ndf = xd1(1); nds = xd2(1);
        maxFirst = max(maxFirst, abs(first - D_nl_11));
    end
    rows(k, :) = [s, probe, m, L, B, T, C_B, rho, T1, M11, D11, k_u, u, Xuu, ...
                  D_nl_11, X, first, second, ndf, nds];
end
fprintf('Instrumented vs unmodified osv (xdot), max |diff| = %.3g\n', maxSelf);
fprintf('Set 1: formula vs osv.m first call D_nl_11, max |diff| = %.3g\n', maxFirst);
fprintf('Set 1: D11 = %.17g, M11 = %.17g\n', D11_osv, ws0.vessel.M(1, 1));

outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outDir, 'hydrodynamics');
if ~exist(outDir, 'dir'), mkdir(outDir); end
out = fullfile(outDir, ['osv_surge_damping_mss_' tag '.csv']);
fid = fopen(out, 'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
fmt = [strjoin(repmat({'%.17g'}, 1, numel(hdr)), ','), '\n'];
for k = 1:size(rows, 1)
    fprintf(fid, fmt, rows(k, :));
end
fclose(fid);
fprintf('Saved %d rows x %d columns to %s\n', size(rows, 1), numel(hdr), out);

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
