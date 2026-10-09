%% ================================================================
%  Otter catamaran template reference (MATLAB, full debug)
%
%  Template generator: an inline copy of the Otter model (otter_debug below,
%  written from MSS otter.m as it stood in 2026-02) that also returns every
%  intermediate term; MSS functions it calls (Smtrx, Hmtrx, Rzyx, gravity,
%  crossFlowDrag, m2c, ...) come from MSS_DIR. Latitude gravity
%  gravity(63.446827 deg), the inertia about the CO shifted a second time by H,
%  20 strip end points in the cross-flow sum (the template settings).
%  Produces full/catamaran_dynamics_full_debug_template.csv and the column
%  extracts hydrodynamics/catamaran_matlab_reference.csv and
%  rigid_body/matlab_reference.csv (the frozen files' columns, cut as text).
%  The frozen files were made with MSS of 2026-02; at MSS cc07579 the
%  hydrodynamics extract differs in tau_crossflow (crossFlowDrag.m and
%  cylinderDrag.m changed since) and the rigid-body extract in C_RB at 1e-13
%  (see the data folder's SOURCE.md); they stay frozen as the template record.
%  Inputs: inputs.csv beside this file. Output under OUT_DIR (default
%  <tempdir>/more_mss_references; never this folder, so a run cannot overwrite a
%  frozen file). Run: MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%
%  Author:    Enio Krizman
%  Date:      2026-02-18
% ================================================================
clear functions;
rehash toolboxcache;
mssDir = getenv('MSS_DIR'); if isempty(mssDir), error('MSS_DIR is not set: point it at the MSS checkout root.'); end; addpath(genpath(mssDir));


disp("⚙️ Running OtterUSV cross-validation (full debug mode)...");

% --- Load the frozen test inputs ---
inputs = readtable(fullfile(fileparts(mfilename('fullpath')), 'inputs.csv'))
N = height(inputs);

mp = 25;
rp = [0.05 0 -0.35]';
Vc = 0.3;
beta_c = deg2rad(30);

% --- Prepare storage ---
allRecords = table();

for k = 1:N
    x = table2array(inputs(k, 1:12))';
    tau = table2array(inputs(k, 13:18))';

    % --- Run Otter debug model ---
    [xdot, data] = otter_debug(x, [0;0], mp, rp, Vc, beta_c, tau);

    record = struct();

    % --- Inputs ---
    for i = 1:12, record.(sprintf('x%d', i)) = x(i); end
    for i = 1:6,  record.(sprintf('tau%d', i)) = tau(i); end

    % --- Dynamics parameters ---
    record = flatten_vector(data.nu_c, "nu_c", record);
    record = flatten_vector(data.nu_c_dot, "nu_c_dot", record);
    record = flatten_vector(data.nu_r, "nu_r", record);
    record.U = data.U;
    record.u_c_scalar = data.u_c;
    record.v_c_scalar = data.v_c_;

    % --- Forces ---
    record = flatten_vector(data.tau_damp, "tau_damp", record);
    record = flatten_vector(data.tau_crossflow, "tau_crossflow", record);
    record = flatten_vector(data.tau_total, "tau_total", record);
    
    % --- Matrices ---
    matrices = {'MRB','MA','M','CRB','CA','C','G','D'};
    prefixes = {'M_RB','M_A','M_total','C_RB','C_A','C_total','G','D'};
    for j = 1:numel(matrices)
        Mval = data.(matrices{j});
        Mflat = flatten_matrix(Mval, prefixes{j});
        fn = fieldnames(Mflat);
        for f = 1:numel(fn)
            record.(fn{f}) = Mflat.(fn{f});
        end
    end

    % --- Trim and outputs ---
    record = flatten_vector(data.eta_trim, "eta_trim", record);
    record = flatten_vector(xdot(1:6), "nu_dot", record);
    record = flatten_vector(xdot(7:12), "eta_dot", record);

    % --- Append ---
    allRecords = [allRecords; struct2table(record)];

    % === PRINT FIRST CASE DEBUG ===
if k == 1
    format long g; % ensure full precision output (up to 15–16 digits)

    disp("==================================================")
    disp("🔍 FIRST CASE DEBUG OUTPUT (before CSV flattening, full precision)")
    disp("==================================================")

    fprintf("Case index: %d\n", k);

    fprintf("nu = [");
    fprintf(' %.15f', x(1:6));
    fprintf(" ]\n");

    fprintf("eta = [");
    fprintf(' %.15f', x(7:12));
    fprintf(" ]\n");

    fprintf("tau = [");
    fprintf(' %.15f', tau);
    fprintf(" ]\n\n");

    disp("---- x_total ----");
    disp(num2str(xdot, '%.15f'));

    disp("==================================================")
end

end

% --- Output folder: parent/compare_results/compare_full_test (relative to this .m file) ---
thisFile  = mfilename('fullpath');     % full path to this script
scriptDir = fileparts(thisFile);       % folder where this script lives
outRoot = getenv('OUT_DIR'); if isempty(outRoot), outRoot = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outRoot, 'full');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

outFileMat = fullfile(outDir, 'catamaran_dynamics_full_debug_template.csv');

% --- Write MATLAB CSV ---
writetable(allRecords, outFileMat);
disp("Saved MATLAB results to: " + string(outFileMat));
dataDir = fileparts(scriptDir);
for target = {'hydrodynamics/catamaran_matlab_reference.csv', 'rigid_body/matlab_reference.csv'}
    [d, ~] = fileparts(target{1}); if ~exist(fullfile(outRoot, d), 'dir'), mkdir(fullfile(outRoot, d)); end
    extract_like(outFileMat, fullfile(dataDir, target{1}), fullfile(outRoot, target{1}));
end


%% === Helper: Flatten a vector into struct fields ===
function S = flatten_vector(vec, prefix, S)
    if nargin < 3, S = struct(); end
    vec = vec(:);
    for i = 1:numel(vec)
        S.(sprintf('%s_%02d', prefix, i)) = vec(i);
    end
end

%% === Helper: Flatten a matrix ===
function S = flatten_matrix(M, prefix)
    M = M';      
    M = M(:);
    S = struct();
    for i = 1:numel(M)
        S.(sprintf('%s_%02d', prefix, i)) = M(i);
    end
end

%% === otter_debug (extended) ===
function [xdot, data] = otter_debug(x, n, mp, rp, V_c, beta_c, tau)
% Extended version of otter_debug that returns all dynamics parameters.

% --- Constants ---
mu = deg2rad(63.446827); % Lattitude for Trondheim, Norway (deg)
g = gravity(mu);      % Gravity vector (m/s2)
rho = 1025; L = 2.0; B = 1.08; m = 55.0;
rg = [0.2 0 -0.2]';
R44 = 0.4 * B; R55 = 0.25 * L; R66 = 0.25 * L;
T_sway = 1; T_yaw = 1; Umax = 6 * 0.5144;
B_pont = 0.25; y_pont = 0.395; Cw_pont = 0.75; Cb_pont = 0.4;

% --- States ---
nu = x(1:6); eta = x(7:12); nu2 = nu(4:6);
U = sqrt(nu(1)^2 + nu(2)^2 + nu(3)^2);
u_c = V_c * cos(beta_c - eta(6)); v_c = V_c * sin(beta_c - eta(6));
nu_c = [u_c v_c 0 0 0 0]';
nu_r = nu - nu_c;
nu_c_dot = [-Smtrx(nu2)*nu_c(1:3); zeros(3,1)];

% --- Geometry and inertia ---
nabla = (m+mp)/rho; T = nabla / (2*Cb_pont*B_pont*L);
Ig_CG = m * diag([R44^2, R55^2, R66^2]);
rg = (m*rg + mp*rp)/(m+mp);
Ig = Ig_CG - m*Smtrx(rg)^2 - mp*Smtrx(rp)^2;
I3 = eye(3); O3 = zeros(3);
MRB_CG = [(m+mp)*I3 O3; O3 Ig];
CRB_CG = [(m+mp)*Smtrx(nu2) O3; O3 -Smtrx(Ig*nu2)];
H = Hmtrx(rg);
MRB = H'*MRB_CG*H;
CRB = H'*CRB_CG*H;

% --- Added mass and Coriolis ---
Xudot = -addedMassSurge(m,L,rho);
Yvdot = -1.5*m; Zwdot = -1.0*m;
Kpdot = -0.2*Ig(1,1); Mqdot = -0.8*Ig(2,2); Nrdot = -1.7*Ig(3,3);
MA = -diag([Xudot, Yvdot, Zwdot, Kpdot, Mqdot, Nrdot]);
CA = m2c(MA, nu_r);
M = MRB + MA; C = CRB + CA;

% --- Hydrostatics and restoring ---
Aw_pont = Cw_pont * L * B_pont;
I_T = 2*(1/12)*L*B_pont^3*(6*Cw_pont^3/((1+Cw_pont)*(1+2*Cw_pont))) + 2*Aw_pont*y_pont^2;
I_L = 0.8*2*(1/12)*B_pont*L^3;
KB = (1/3)*(5*T/2 - 0.5*nabla/(L*B_pont));
BM_T = I_T/nabla; BM_L = I_L/nabla;
KM_T = KB + BM_T; KM_L = KB + BM_L; KG = T - rg(3);
GM_T = KM_T - KG; GM_L = KM_L - KG;
G33 = rho*g*(2*Aw_pont); G44 = rho*g*nabla*GM_T; G55 = rho*g*nabla*GM_L;
G_CF = diag([0 0 G33 G44 G55 0]);
LCF = -0.2; H = Hmtrx([LCF 0 0]);
G = H'*G_CF*H;

% --- Damping ---
w3 = sqrt(G33/M(3,3)); w4 = sqrt(G44/M(4,4)); w5 = sqrt(G55/M(5,5));
Xu = -24.4*g/Umax; Yv = -M(2,2)/T_sway; Zw = -2*0.3*w3*M(3,3);
Kp = -2*0.2*w4*M(4,4); Mq = -2*0.4*w5*M(5,5); Nr = -M(6,6)/T_yaw;
D = diag([Xu, Yv, Zw, Kp, Mq, Nr]);

% --- Forces ---
tau_damp = [Xu*nu_r(1), Yv*nu_r(2), Zw*nu_r(3), Kp*nu_r(4), Mq*nu_r(5), Nr*(1+10*abs(nu_r(6)))*nu_r(6)]';
tau_crossflow = crossFlowDrag(L,B_pont,T,nu_r);
tau_total = tau + tau_damp + tau_crossflow;

% --- Trim and payload ---
f_payload = Rzyx(eta(4),eta(5),eta(6))'*[0 0 mp*g]';
m_payload = Smtrx(rp)*f_payload;
g_0 = [f_payload; m_payload];
eta_0 = [0;0;inv(G(3:5,3:5))*g_0(3:5);0];
eta_trim = eta - eta_0;

J = eulerang(eta_trim(4),eta_trim(5),eta_trim(6));
xdot = [nu_c_dot + M\(tau_total - C*nu_r - G*eta_trim);
        J*nu];

% --- Return all data ---
data.nu_c = nu_c;
data.nu_c_dot = nu_c_dot;
data.nu_r = nu_r;
data.U = U;
data.u_c = u_c;
data.v_c_ = v_c;
data.tau_damp = tau_damp;
data.tau_crossflow = tau_crossflow;
data.tau_total = tau_total;
data.MRB = MRB;
data.MA = MA;
data.M = M;
data.CRB = CRB;
data.CA = CA;
data.C = C;
data.G = G;
data.D = D;
data.eta_trim = eta_trim;
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
