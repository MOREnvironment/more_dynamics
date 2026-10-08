%% ================================================================
%  Rigid-body + added-mass forms, computed by MSS
%
%  MSS line numbers below are those of MSS ac77394 (2026-10-06), where the
%  CSVs were first made; byte-identical at 72656d1 (release 2.0.2) and at
%  cc07579 (m2c.m's 3-DOF branch was corrected in a3406cf, but the 3-DOF
%  mass matrices here have no surge coupling, M(1,2) = M(1,3) = 0, where the
%  old and the corrected branch agree; see the data folder's SOURCE.md).
%
%  Reference data for every rigid-body / added-mass form of more_dynamics
%  that has an MSS counterpart. Every matrix comes from an MSS function called as
%  written (LIBRARY/modeling and LIBRARY/kinematics):
%    rbody.m          MRB, CRB(nu2) about CO           (generic, hull)
%    m2c.m            C(nu) from M, 6-DOF and 3-DOF    (Lagrangian forms, CA)
%    spheroid.m       MRB, CRB(nu2) of a prolate spheroid (rho = 1025, line 35)
%    imlay61.m        MA, CA, Lamb k-factors           (rho = 1026, line 31)
%    addedMassSurge.m A11 (Soeding 1982)               (hull surge added mass)
%  Transcribed MSS lines (no MSS function holds them):
%    remus100.m 205-208  the CA entries zeroed in remus100 ("stabilized" CA)
%    otter.m 153-158,160 the scaled-derivative added mass, MA = -diag(...),
%                        rotational terms scaled by Ig about the CG
%    hull mass m = rho*Cb*L*B*T (definition of the block coefficient)
%  Parameter sets are TEST INPUTS (not vehicle data); each output row echoes
%  its inputs so the Python tests read them from the CSV, never retype them.
%  States: cases 1-50 = nu of the frozen 50 inputs (x1..x6 of
%  inputs.csv beside this file); cases 51-200 = seeded uniform(-3,3)
%  (rng(20261006,'twister')). The same nu is used for C_RB and C_A.
%  Numbers are written with %.17g (round-trip double precision).
%  Outputs, under <OUT_DIR>/rigid_body/:
%    rigid_body_rbody_mss_current.csv      (rbody + m2c(MRB))
%    rigid_body_spheroid_mss_current.csv   (spheroid + imlay61)
%    rigid_body_hull_mss_current.csv       (rbody + addedMassSurge + m2c)
%    rigid_body_m2c_3dof_mss_current.csv   (m2c 3-DOF)
%    rigid_body_m2c_3dof_coupled_mss_cc07579.csv
%                                          (m2c 3-DOF on fully coupled M3;
%                                           section 5, made at cc07579 only:
%                                           before a3406cf the branch dropped
%                                           the surge coupling)
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
% MSS is found only through MSS_DIR (no path relative to this checkout);
% unset -> stop with a message naming the variable.
mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: export MSS_DIR=<path to the MSS checkout>');
end
addpath(genpath(mssDir));
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
for sub = {'rigid_body'}
    if ~exist(fullfile(outDir, sub{1}), 'dir'), mkdir(fullfile(outDir, sub{1})); end
end
fprintf('MSS_DIR: %s\n', mssDir);
fprintf('MATLAB %s\n', version);
fprintf('rbody: %s\nm2c: %s\nspheroid: %s\nimlay61: %s\naddedMassSurge: %s\n', ...
    which('rbody'), which('m2c'), which('spheroid'), which('imlay61'), which('addedMassSurge'));

% ------------------------------------------------
% States
% ------------------------------------------------
inputs = readmatrix(fullfile(scriptDir, 'inputs.csv'));
nu_frozen = inputs(:, 1:6);
assert(size(nu_frozen, 1) == 50);
rng(20261006, 'twister');
nu_seeded = -3 + 6 * rand(150, 6);
NU = [nu_frozen; nu_seeded];
NCASE = size(NU, 1);

% ------------------------------------------------
% 1) Generic rigid body: rbody.m, m2c(MRB, nu)
% ------------------------------------------------
rb_sets = struct( ...
    'm',   {120.0, 31.9, 55.0}, ...
    'R',   {[0.45 0.70 0.75], [0.06 0.41 0.41], [0.432 0.50 0.50]}, ...
    'r_bG',{[0.12 -0.04 0.21], [0.0 0.0 0.02], [0.0 0.0 0.0]});
rows = {};
for s = 1:numel(rb_sets)
    p = rb_sets(s);
    for k = 1:NCASE
        nu = NU(k, :)';
        [MRB, CRB] = rbody(p.m, p.R(1), p.R(2), p.R(3), nu(4:6), p.r_bG');
        CRB_lag = m2c(MRB, nu);
        rows{end+1} = [s, k, p.m, p.R, p.r_bG, nu', rowmajor(MRB), rowmajor(CRB), rowmajor(CRB_lag)]; %#ok<SAGROW>
    end
end
hdr = [{'set', 'case', 'm', 'R44', 'R55', 'R66', 'r_bG_1', 'r_bG_2', 'r_bG_3'}, ...
       vec_names('nu', 6), mat_names('M_RB', 36), mat_names('C_RB', 36), mat_names('C_RB_m2c', 36)];
write_csv(fullfile(outDir, 'rigid_body', 'rigid_body_rbody_mss_current.csv'), hdr, cell2mat(rows'));

% ------------------------------------------------
% 2) Prolate spheroid: spheroid.m, imlay61.m (+ remus100.m 205-208)
% ------------------------------------------------
L_auv = 1.6; D_auv = 0.19;               % remus100.m lines 132-133
sp_sets = struct( ...
    'a',   {1.0096 * L_auv/2, 1.2}, ...  % remus100.m line 135
    'b',   {1.0096 * D_auv/2, 0.15}, ... % remus100.m line 136
    'r44', {0.3, 0.0}, ...               % remus100.m line 137; set 2: imlay61 default (nargin 3)
    'r_bG',{[0 0 0.02], [0.03 -0.01 0.04]});   % remus100.m line 138; set 2 test input
rows = {};
for s = 1:numel(sp_sets)
    p = sp_sets(s);
    for k = 1:NCASE
        nu = NU(k, :)';
        [MRB, CRB] = spheroid(p.a, p.b, nu(4:6), p.r_bG');
        if s == 2
            [MA, CA] = imlay61(p.a, p.b, nu);          % nargin == 3: r44 = 0
        else
            [MA, CA] = imlay61(p.a, p.b, nu, p.r44);
        end
        CA_stab = remus100_zeroing(CA);
        rows{end+1} = [s, k, p.a, p.b, p.r44, p.r_bG, nu', ...
                       rowmajor(MRB), rowmajor(CRB), rowmajor(MA), rowmajor(CA), rowmajor(CA_stab)]; %#ok<SAGROW>
    end
end
hdr = [{'set', 'case', 'a', 'b', 'r44', 'r_bG_1', 'r_bG_2', 'r_bG_3'}, vec_names('nu', 6), ...
       mat_names('M_RB', 36), mat_names('C_RB', 36), mat_names('M_A', 36), mat_names('C_A', 36), ...
       mat_names('C_A_stab', 36)];
write_csv(fullfile(outDir, 'rigid_body', 'rigid_body_spheroid_mss_current.csv'), hdr, cell2mat(rows'));

% ------------------------------------------------
% 3) Displacement hull: m = rho*Cb*L*B*T, rbody.m, addedMassSurge.m,
%    otter.m 153-158/160 pattern, m2c.m (+ remus100.m 205-208)
% ------------------------------------------------
hull_sets = struct( ...
    'rho', {1025.0, 1000.0}, 'L', {5.2, 3.0}, 'B', {2.15, 1.2}, 'T', {0.3, 0.4}, ...
    'Cb',  {0.233, 0.45}, ...
    'Rs',  {[0.35 0.25 0.25], [0.38 0.27 0.24]}, ...
    'r_bG',{[0 0 0.025], [0.15 0 -0.10]}, ...
    'c',   {[-1.0 -1.5 -1.0 -0.2 -0.8 -1.2], [-1.2 -1.1 -0.9 -0.3 -0.7 -1.0]});
rows = {};
for s = 1:numel(hull_sets)
    p = hull_sets(s);
    m = p.rho * p.Cb * p.L * p.B * p.T;
    R44 = p.Rs(1) * p.B; R55 = p.Rs(2) * p.L; R66 = p.Rs(3) * p.L;
    Ig = m * diag([R44^2, R55^2, R66^2]);        % about the CG (rbody.m line 35)
    A11 = addedMassSurge(m, p.L, p.rho);
    d = p.c .* [A11, m, m, Ig(1,1), Ig(2,2), Ig(3,3)];
    MA = -diag(d);                               % otter.m line 160 pattern
    for k = 1:NCASE
        nu = NU(k, :)';
        [MRB, CRB] = rbody(m, R44, R55, R66, nu(4:6), p.r_bG');
        CA = m2c(MA, nu);
        CA_stab = remus100_zeroing(CA);
        rows{end+1} = [s, k, p.rho, p.L, p.B, p.T, p.Cb, p.Rs, p.r_bG, p.c, nu', m, A11, ...
                       rowmajor(MRB), rowmajor(CRB), rowmajor(MA), rowmajor(CA), rowmajor(CA_stab)]; %#ok<SAGROW>
    end
end
hdr = [{'set', 'case', 'rho', 'L', 'B', 'T', 'Cb', 'Rs_1', 'Rs_2', 'Rs_3', 'r_bG_1', 'r_bG_2', 'r_bG_3'}, ...
       vec_names('c', 6), vec_names('nu', 6), {'m', 'A11'}, ...
       mat_names('M_RB', 36), mat_names('C_RB', 36), mat_names('M_A', 36), mat_names('C_A', 36), ...
       mat_names('C_A_stab', 36)];
write_csv(fullfile(outDir, 'rigid_body', 'rigid_body_hull_mss_current.csv'), hdr, cell2mat(rows'));

% ------------------------------------------------
% 4) 3-DOF: m2c.m with a 3x3 M = -[[X_du 0 0]; [0 Y_dv Y_dr]; [0 Y_dr N_dr]]
%    (nu3 = [u v r] = nu([1 2 6]))
% ------------------------------------------------
dof3_sets = struct('X_du', {-5.0, -12.0}, 'Y_dv', {-60.0, -150.0}, ...
                   'Y_dr', {-4.0, 7.5}, 'N_dr', {-30.0, -220.0});
rows = {};
for s = 1:numel(dof3_sets)
    p = dof3_sets(s);
    M3 = -[p.X_du 0 0; 0 p.Y_dv p.Y_dr; 0 p.Y_dr p.N_dr];
    for k = 1:NCASE
        nu3 = NU(k, [1 2 6])';
        C3 = m2c(M3, nu3);
        rows{end+1} = [s, k, p.X_du, p.Y_dv, p.Y_dr, p.N_dr, nu3', rowmajor(M3), rowmajor(C3)]; %#ok<SAGROW>
    end
end
hdr = [{'set', 'case', 'X_du', 'Y_dv', 'Y_dr', 'N_dr'}, vec_names('nu3', 3), ...
       mat_names('M3', 9), mat_names('C3', 9)];
write_csv(fullfile(outDir, 'rigid_body', 'rigid_body_m2c_3dof_mss_current.csv'), hdr, cell2mat(rows'));

% ------------------------------------------------
% 5) 3-DOF: m2c.m on fully coupled symmetric M3 (every off-diagonal entry
%    non-zero), MSS cc07579 lines 52-57: p = M*nu, C = [0 0 -p2; 0 0 p1;
%    p2 -p1 0]. Set 1 is M3 = [10 2 1; 2 20 3; 1 3 30]; set 2 a second
%    positive-definite test matrix. Case 0 = nu3 = [1 2 3] (C*nu3 =
%    [-153; 51; 17] for set 1); cases 1-200 = NU(k, [1 2 6]).
% ------------------------------------------------
coupled_sets = {[10 2 1; 2 20 3; 1 3 30], [25 -1.5 4; -1.5 60 7.5; 4 7.5 220]};
rows = {};
for s = 1:numel(coupled_sets)
    M3 = coupled_sets{s};
    nu3 = [1; 2; 3];
    rows{end+1} = [s, 0, nu3', rowmajor(M3), rowmajor(m2c(M3, nu3))]; %#ok<SAGROW>
    for k = 1:NCASE
        nu3 = NU(k, [1 2 6])';
        rows{end+1} = [s, k, nu3', rowmajor(M3), rowmajor(m2c(M3, nu3))]; %#ok<SAGROW>
    end
end
hdr = [{'set', 'case'}, vec_names('nu3', 3), mat_names('M3', 9), mat_names('C3', 9)];
write_csv(fullfile(outDir, 'rigid_body', 'rigid_body_m2c_3dof_coupled_mss_cc07579.csv'), hdr, cell2mat(rows'));

fprintf('done: %d cases per set\n', NCASE);

% ================================================================
function CA = remus100_zeroing(CA)
% remus100.m lines 205-208, transcribed (no MSS function holds them).
CA(5,3) = 0; CA(3,5) = 0;
CA(5,1) = 0; CA(1,5) = 0;
CA(6,1) = 0; CA(1,6) = 0;
CA(6,2) = 0; CA(2,6) = 0;
end

function v = rowmajor(M)
M = M.';
v = M(:)';
end

function names = vec_names(prefix, n)
names = arrayfun(@(i) sprintf('%s_%d', prefix, i), 1:n, 'UniformOutput', false);
end

function names = mat_names(prefix, n)
names = arrayfun(@(i) sprintf('%s_%02d', prefix, i), 1:n, 'UniformOutput', false);
end

function write_csv(path, header, data)
assert(numel(header) == size(data, 2), 'header/data width mismatch: %s', path);
fid = fopen(path, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
fmt = [strjoin(repmat({'%.17g'}, 1, size(data, 2)), ','), '\n'];
fprintf(fid, fmt, data.');
fclose(fid);
fprintf('%s: %d x %d\n', path, size(data, 1), size(data, 2));
end
