%% ================================================================
%  Otter differential thruster sweep, computed by MSS
%
%  Every constant and the thrust law come from MSS (T. I. Fossen, MIT):
%    - B_prop, n_min, n_max from [~,~,~,B_prop,n_min,n_max] = otter()
%      (CRAFT/USV/models/otter.m, nargin == 0 branch), so g = 9.81.
%    - allocation as CRAFT/USV/SIMotter.m (lines 95, 183-184 at 99bf0b3):
%      Binv = invQR(B_prop); u = Binv*[tau_X; tau_N]; n_c = sign(u).*sqrt(abs(u))
%    - saturation and thrust from otter.m itself (lines 220-232 at 99bf0b3),
%      read from a copy instrumented at run time (in tempdir, never saved)
%      that only returns its workspace; checked against the unmodified otter.m.
%  The sweep (requested tau_X, tau_N over 120 s in 0.1 s steps: ramps, steps
%  and reversals, ahead and astern) is a test construction.
%  First made at MSS 99bf0b3 (2026-10-05); byte-identical at 72656d1 and
%  cc07579 (otter.m and SIMotter.m changed only comments and the KB line,
%  which this sweep does not read).
%  Output: <OUT_DIR>/force_producers/differential_thruster_mss_current.csv
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
mssDir = getenv('MSS_DIR');   % MSS is found only through MSS_DIR, never a sibling folder
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout root.');
end
addpath(genpath(mssDir));
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
for sub = {'force_producers'}
    if ~exist(fullfile(outDir, sub{1}), 'dir'), mkdir(fullfile(outDir, sub{1})); end
end
fprintf('MATLAB %s\nMSS: %s\n', version, which('otter'));

otter_ws = mss_instrumented('otter');

[~, ~, ~, B_prop, n_min, n_max] = otter();
Binv = invQR(B_prop);
fprintf('n_max = %.12f rad/s, n_min = %.12f rad/s\n', n_max, n_min);
disp('B_prop ='); disp(B_prop);

% ------------------------------------------------
% Sweep (test construction)
% ------------------------------------------------
dt = 0.1;
t_phase = 10;
t_total = 3 * 4 * t_phase;
t = 0:dt:t_total;

ramp_profile = @(t, T, A) arrayfun(@(x) ...
    (x<T)*A*(x/T) + ...
    (x>=T && x<2*T)*(A - 2*A*((x-T)/T)) + ...
    (x>=2*T && x<3*T)*(-A + 2*A*((x-2*T)/T)), t);

tau_X_ramp = ramp_profile(t, t_phase, 150);
tau_N_ramp = ramp_profile(t, t_phase, 150);

tau_X = zeros(size(t));
tau_N = zeros(size(t));

phase1 = round(4*t_phase/dt);
phase2 = 2*phase1;

tau_X(1:phase1) = tau_X_ramp(1:phase1);
tau_N(phase1+1:phase2) = tau_N_ramp(1:phase1);
tau_X(phase2+1:end) = tau_X_ramp(1:length(t)-phase2);
tau_N(phase2+1:end) = tau_N_ramp(1:length(t)-phase2);

% ------------------------------------------------
% Allocation and thrust
% ------------------------------------------------
x0 = zeros(12,1); mp = 25; rp = zeros(3,1);   % do not enter the thrust law
nT = numel(t);
n_left = zeros(1,nT);  n_right = zeros(1,nT);
thrust_left = zeros(1,nT); thrust_right = zeros(1,nT);
tau_X_out = zeros(1,nT); tau_N_out = zeros(1,nT);
maxSelfCheck = 0;

for i = 1:nT
    u = Binv * [tau_X(i); tau_N(i)];
    n_c = sign(u) .* sqrt(abs(u));

    xdot0 = otter(x0, n_c, mp, rp, 0, 0);
    [xdot1, ~, ~, ~, ~, ~, ws] = otter_ws(x0, n_c, mp, rp, 0, 0);
    maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0)]);

    n_left(i) = ws.n(1);   n_right(i) = ws.n(2);      % after satlim (line 220)
    thrust_left(i) = ws.Thrust(1); thrust_right(i) = ws.Thrust(2);
    tau_X_out(i) = ws.tau(1);  tau_N_out(i) = ws.tau(6);
end

fprintf('Instrumented vs unmodified otter (xdot), max |diff| = %.3g\n', maxSelfCheck);
fprintf('gravity used by otter: %.10f m/s2\n', ws.g);

T = table(t', tau_X', tau_N', n_left', n_right', thrust_left', thrust_right', ...
          tau_X_out', tau_N_out', ...
          'VariableNames', {'time','tau_X_cmd','tau_N_cmd','n_cmd_left','n_cmd_right', ...
                            'thrust_left','thrust_right','tau_X','tau_N'});
outFile = fullfile(outDir, 'force_producers', 'differential_thruster_mss_current.csv');
writetable(T, outFile);
fprintf('Saved %d rows x %d columns to %s\n', height(T), width(T), outFile);

%% ================================================================
% Helper
% ================================================================
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
