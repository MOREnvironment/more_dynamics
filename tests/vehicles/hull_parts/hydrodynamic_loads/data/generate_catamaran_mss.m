%% ================================================================
%  Otter-based catamaran reference data, computed by MSS otter.m
%
%  MSS line numbers below are those of MSS 99bf0b3 (2026-10-05), where the
%  CSV was first made; byte-identical at 72656d1 (release 2.0.2). At cc07579
%  otter.m's KB uses the pontoon waterplane area (commit cf349d4), which
%  changes G44, G55 and through them Kp, Mq (see the data folder's SOURCE.md).
%
%  Every term comes from CRAFT/USV/models/otter.m (T. I. Fossen, MSS, MIT):
%    - otter.m is run as written; a copy instrumented at run time (in
%      tempdir, never saved) only returns its workspace so the intermediate
%      terms (MRB, MA, G, Xu..Nr, tau_damp, ...) can be read.
%    - the instrumented copy is checked against the unmodified otter.m
%      (xdot and M) on every case.
%    - otter.m takes propeller speeds, not tau; with n = [0 0] its propeller
%      tau is zero, so the external tau of the inputs is added in the otter.m
%      state equation (lines 261-263).
%    - D = diag([Xu Yv Zw Kp Mq Nr]) of the terms otter.m forms (lines
%      202-207); otter.m builds no D matrix.
%  Inputs: inputs.csv beside this file (50 frozen states and wrenches),
%  payload mp = 25 kg at rp = [0.05 0 -0.35] m, current 0.3 m/s at 30 deg.
%  Outputs, under <OUT_DIR>:
%    full/catamaran_dynamics_full_debug_mss.csv            (every column)
%    hydrodynamics/catamaran_matlab_reference_mss_<tag>.csv (its columns;
%        tag = current for the L*B_pont KB, cc07579 for the Aw_pont KB)
%    rigid_body/matlab_reference_mss_current.csv           (its columns)
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
for sub = {'full', 'hydrodynamics', 'rigid_body'}
    if ~exist(fullfile(outDir, sub{1}), 'dir'), mkdir(fullfile(outDir, sub{1})); end
end
fprintf('MATLAB %s\nMSS: %s\n', version, which('otter'));

otter_ws = mss_instrumented('otter');

% ------------------------------------------------
% Inputs and constants (payload and current of the frozen reference)
% ------------------------------------------------
inputs = readtable(fullfile(scriptDir, 'inputs.csv'));
N = height(inputs);
mp = 25;
rp = [0.05 0 -0.35]';
Vc = 0.3;
beta_c = deg2rad(30);
n = [0; 0];   % propellers off

allRecords = table();
maxSelfCheck = 0;

for k = 1:N
    x   = table2array(inputs(k, 1:12))';
    tau = table2array(inputs(k, 13:18))';

    [xdot0, ~, M0] = otter(x, n, mp, rp, Vc, beta_c);
    [xdot1, ~, M1, ~, ~, ~, ws] = otter_ws(x, n, mp, rp, Vc, beta_c);
    if any(ws.tau ~= 0)
        error('otter propeller tau is not zero with n = 0');
    end
    maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0); abs(M1(:) - M0(:))]);

    tau_total = tau + ws.tau_damp + ws.tau_crossflow;
    D = diag([ws.Xu, ws.Yv, ws.Zw, ws.Kp, ws.Mq, ws.Nr]);
    eta_trim = ws.eta;   % otter.m line 255: eta = eta - eta_0

    % otter.m lines 261-263 with the external tau of the inputs
    nu_dot  = ws.nu_c_dot + ws.M \ (tau_total - ws.C * ws.nu_r - ws.G * eta_trim);
    eta_dot = ws.J * ws.nu;

    record = struct();
    for i = 1:12, record.(sprintf('x%d', i)) = x(i); end
    for i = 1:6,  record.(sprintf('tau%d', i)) = tau(i); end

    record = flatten_vector(ws.nu_c,     "nu_c",     record);
    record = flatten_vector(ws.nu_c_dot, "nu_c_dot", record);
    record = flatten_vector(ws.nu_r,     "nu_r",     record);
    record.U = ws.U;
    record.u_c_scalar = ws.u_c;
    record.v_c_scalar = ws.v_c;

    record = flatten_vector(ws.tau_damp,      "tau_damp",      record);
    record = flatten_vector(ws.tau_crossflow, "tau_crossflow", record);
    record = flatten_vector(tau_total,        "tau_total",     record);

    mats = {ws.MRB, ws.MA, ws.M, ws.CRB, ws.CA, ws.C, ws.G, D};
    prefixes = {'M_RB','M_A','M_total','C_RB','C_A','C_total','G','D'};
    for j = 1:numel(mats)
        record = flatten_matrix(mats{j}, prefixes{j}, record);
    end

    record = flatten_vector(eta_trim, "eta_trim", record);
    record = flatten_vector(nu_dot,   "nu_dot",   record);
    record = flatten_vector(eta_dot,  "eta_dot",  record);

    allRecords = [allRecords; struct2table(record)]; %#ok<AGROW>
end

fprintf('Instrumented vs unmodified otter (xdot, M), max |diff| = %.3g\n', maxSelfCheck);
fprintf('gravity used by otter: %.10f m/s2, rho in otter: %g\n', ws.g, ws.rho);
fprintf('Ig diag (kg m2): %.6f %.6f %.6f\n', diag(ws.Ig));

outFile = fullfile(outDir, 'full', 'catamaran_dynamics_full_debug_mss.csv');
writetable(allRecords, outFile);
fprintf('Saved %d rows x %d columns to %s\n', height(allRecords), width(allRecords), outFile);
otterText = fileread(which('otter'));
if contains(otterText, '0.5*nabla/(L*B_pont)')
    catTag = 'current';      % KB with L*B_pont (MSS up to 72656d1)
elseif contains(otterText, '0.5*nabla/Aw_pont')
    catTag = 'cc07579';      % KB with the pontoon waterplane area (cf349d4)
else
    error('otter.m: KB line not recognised; revisit the output names');
end
dataDir = fileparts(scriptDir);
extract_like(outFile, fullfile(dataDir, 'hydrodynamics', 'catamaran_matlab_reference_mss_current.csv'), ...
    fullfile(outDir, 'hydrodynamics', ['catamaran_matlab_reference_mss_' catTag '.csv']));
extract_like(outFile, fullfile(dataDir, 'rigid_body', 'matlab_reference_mss_current.csv'), ...
    fullfile(outDir, 'rigid_body', 'matlab_reference_mss_current.csv'));

%% ================================================================
% Helpers
% ================================================================
function S = flatten_vector(vec, prefix, S)
    vec = vec(:);
    for i = 1:numel(vec)
        S.(sprintf('%s_%02d', prefix, i)) = vec(i);
    end
end

function S = flatten_matrix(M, prefix, S)
    M = M.';   % row-major
    M = M(:);
    for i = 1:numel(M)
        S.(sprintf('%s_%02d', prefix, i)) = M(i);
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

function extract_like(fullCsv, frozenCsv, outCsv)
% Columns named in the header of frozenCsv, cut as text from fullCsv in
% fullCsv's column order (no number is re-written).
fid = fopen(frozenCsv, 'r'); want = strsplit(strtrim(fgetl(fid)), ','); fclose(fid);
L = splitlines(strtrim(fileread(fullCsv)));
idx = find(ismember(strsplit(L{1}, ','), want));
if numel(idx) ~= numel(want)
    error('%s: %d of %d columns found', fullCsv, numel(idx), numel(want));
end
fid = fopen(outCsv, 'w');
for k = 1:numel(L)
    f = strsplit(L{k}, ',');
    fprintf(fid, '%s\n', strjoin(f(idx), ','));
end
fclose(fid);
fprintf('Extracted %d columns to %s\n', numel(idx), outCsv);
end
