%% ================================================================
%  Spheroid AUV reference data, computed by MSS remus100.m
%
%  MSS line numbers below are those of MSS 99bf0b3 (2026-10-05), where the
%  CSV was first made. From ac77394 on only tau_crossflow and what follows
%  from it change, by the cross-flow Reynolds number on the diameter
%  (HYDRO/cylinderDrag.m); 72656d1 and cc07579 change nothing on this path.
%
%  Every term comes from CRAFT/AUV/models/remus100.m (T. I. Fossen, MSS,
%  MIT) and the functions it calls:
%    - remus100.m is run as written; a copy instrumented at run time (in
%      tempdir, never saved) only returns its workspace so the intermediate
%      terms (M_RB, M_A, C, D, tau_liftdrag, ...) can be read.
%    - the instrumented copy is checked against the unmodified remus100.m
%      (xdot and M) on every case.
%    - remus100.m takes actuator commands, not tau; with ui = [0 0 0] its
%      actuator tau is zero, so the external tau of the inputs is added in
%      the remus100.m state equation (lines 255-256).
%  MSS's own constants on this path: rho 1026 in imlay61/forceLiftDrag,
%  1025 in spheroid/crossFlowDrag; surge-only damping fade; crossFlowDrag
%  with 20 strip midpoints.
%  Inputs: inputs.csv beside this file (50 frozen states and wrenches).
%  Outputs, under <OUT_DIR>:
%    full/spheroid_auv_dynamics_full_debug_mss.csv                (every column)
%    hydrodynamics/spheroid_matlab_reference_mss_ac77394.csv      (its columns)
%    hydrostatics/spheroid_matlab_reference_mss_current.csv       (its columns)
%    rigid_body/spheroid_matlab_reference_mss_current.csv         (its columns)
%  (hydrodynamics/spheroid_matlab_reference_mss_current.csv is the pre-ac77394
%  cross-flow and cannot be made from MSS ac77394 or later.)
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
for sub = {'full', 'hydrodynamics', 'hydrostatics', 'rigid_body'}
    if ~exist(fullfile(outDir, sub{1}), 'dir'), mkdir(fullfile(outDir, sub{1})); end
end
fprintf('MATLAB %s\nMSS: %s\n', version, which('remus100'));

remus100_ws = mss_instrumented('remus100');

% ------------------------------------------------
% Inputs and constants (same as test_dynamics_consitency.m lines 15, 21-23)
% ------------------------------------------------
inputs = readtable(fullfile(scriptDir, 'inputs.csv'));
N = height(inputs);
Vc     = 0.3;
beta_c = deg2rad(30);
w_c    = 0.0;
ui     = zeros(3,1);   % no rudder, no stern plane, no propeller

allRecords = table();
maxSelfCheck = 0;

for k = 1:N
    x   = table2array(inputs(k, 1:12))';
    tau = table2array(inputs(k, 13:18))';

    [xdot0, ~, M0] = remus100(x, ui, Vc, beta_c, w_c);
    [xdot1, ~, M1, ws] = remus100_ws(x, ui, Vc, beta_c, w_c);
    if any(ws.tau ~= 0)
        error('remus100 actuator tau is not zero with ui = 0');
    end
    maxSelfCheck = max([maxSelfCheck; abs(xdot1 - xdot0); abs(M1(:) - M0(:))]);

    % remus100.m lines 255-256 with the external tau of the inputs
    nu_dot  = ws.Dnu_c + ws.M \ (tau + ws.tau_liftdrag + ws.tau_crossflow ...
              - ws.C * ws.nu_r - ws.D * ws.nu_r - ws.g);
    eta_dot = ws.J * ws.nu;

    record = struct();
    for i = 1:12, record.(sprintf('x%d', i)) = x(i); end
    for i = 1:6,  record.(sprintf('tau%d', i)) = tau(i); end

    record = flatten_vector(ws.nu_c,  "nu_c",  record);
    record = flatten_vector(ws.Dnu_c, "Dnu_c", record);
    record = flatten_vector(ws.nu_r,  "nu_r",  record);
    record.alpha = ws.alpha;
    record.U_r   = ws.U_r;
    record.U     = ws.U;

    record = flatten_vector(ws.tau_liftdrag,  "tau_lift_drag", record);
    record = flatten_vector(ws.tau_crossflow, "tau_crossflow", record);
    record = flatten_vector(tau + ws.tau_liftdrag + ws.tau_crossflow, "tau_total", record);

    mats = {ws.MRB, ws.MA, ws.M, ws.CRB, ws.CA, ws.C, ws.D};
    prefixes = {'M_RB','M_A','M','C_RB','C_A','C','D'};
    for j = 1:numel(mats)
        record = flatten_matrix(mats{j}, prefixes{j}, record);
    end

    record = flatten_vector(ws.g,    "g",       record);
    record = flatten_vector(nu_dot,  "nu_dot",  record);
    record = flatten_vector(eta_dot, "eta_dot", record);

    allRecords = [allRecords; struct2table(record)]; %#ok<AGROW>
end

fprintf('Instrumented vs unmodified remus100 (xdot, M), max |diff| = %.3g\n', maxSelfCheck);
fprintf('gravity used by remus100: %.10f m/s2, rho in remus100: %g\n', ws.g_mu, ws.rho);

outFile = fullfile(outDir, 'full', 'spheroid_auv_dynamics_full_debug_mss.csv');
writetable(allRecords, outFile);
fprintf('Saved %d rows x %d columns to %s\n', height(allRecords), width(allRecords), outFile);
dataDir = fileparts(scriptDir);
for target = {'hydrodynamics/spheroid_matlab_reference_mss_ac77394.csv', ...
              'hydrostatics/spheroid_matlab_reference_mss_current.csv', ...
              'rigid_body/spheroid_matlab_reference_mss_current.csv'}
    extract_like(outFile, fullfile(dataDir, target{1}), fullfile(outDir, target{1}));
end

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
