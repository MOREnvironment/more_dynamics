%% ================================================================
%  Spheroid AUV template reference (MATLAB, full debug)
%
%  Template generator: an inline copy of the REMUS 100 model (remus100_core,
%  imlay611, forceLiftDragg below, written from MSS remus100.m, imlay61.m and
%  forceLiftDrag.m as they stood in 2026-02) that also returns every
%  intermediate term; other MSS functions come from MSS_DIR. Template
%  settings: rho = 1025 everywhere, surge and sway damping fade, cross-flow
%  through MSS crossFlowDrag.m.
%  Produces full/spheroid_auv_dynamics_full_debug_template.csv and the column
%  extracts hydrodynamics/spheroid_matlab_reference.csv and
%  hydrostatics/spheroid_matlab_reference.csv (the frozen files' columns, cut
%  as text). The frozen files were made with MSS of 2026-02; at MSS cc07579
%  the hydrodynamics extract differs in tau_crossflow only (crossFlowDrag.m
%  and cylinderDrag.m changed since); the hydrostatics extract is identical.
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
format long g
format compact

disp("⚙️ Running SpheroidAUV cross-validation (FULL debug mode)...");

% ------------------------------------------------
% Load the frozen test inputs
% ------------------------------------------------
inputs = readtable(fullfile(fileparts(mfilename('fullpath')), 'inputs.csv'));
N = height(inputs);

% ------------------------------------------------
% Constants (the template settings)
% ------------------------------------------------
Vc      = 0.3;
beta_c = deg2rad(30);
w_c     = 0.0;

% ------------------------------------------------
% Storage
% ------------------------------------------------
allRecords = table();

% ================================================================
% Main loop (NO integration, pure algebraic comparison)
% ================================================================
for k = 1:N
    x   = table2array(inputs(k, 1:12))';
    tau = table2array(inputs(k, 13:18))';

    % --- Run REMUS-style core ---
    [xdot, ~, ~, dbg] = remus100_core(x, tau, Vc, beta_c, w_c);

    nu_dot  = xdot(1:6);
    eta_dot = xdot(7:12);

    record = struct();

    % ------------------------------------------------------------
    % Inputs
    % ------------------------------------------------------------
    for i = 1:12
        record.(sprintf('x%d', i)) = x(i);
    end
    for i = 1:6
        record.(sprintf('tau%d', i)) = tau(i);
    end

    % ------------------------------------------------------------
    % Flow parameters
    % ------------------------------------------------------------
    record = flatten_vector(dbg.nu_c,   "nu_c",   record);
    record = flatten_vector(dbg.Dnu_c,  "Dnu_c",  record);
    record = flatten_vector(dbg.nu_r,   "nu_r",   record);

    record.alpha = dbg.alpha;
    record.U_r   = dbg.U_r;
    record.U     = dbg.U;

    % ------------------------------------------------------------
    % Forces
    % ------------------------------------------------------------
    record = flatten_vector(dbg.tau_lift_drag,  "tau_lift_drag", record);
    record = flatten_vector(dbg.tau_cross_flow, "tau_crossflow", record);
    record = flatten_vector(dbg.tau_total,      "tau_total",    record);

    % ------------------------------------------------------------
    % Matrices (FULL SET — THIS WAS THE BUG)
    % ------------------------------------------------------------
    matrices = {'M_RB','M_A','M','C_RB','C_A','C','D'};
    prefixes = {'M_RB','M_A','M','C_RB','C_A','C','D'};

    for j = 1:numel(matrices)
        Mval = dbg.(matrices{j});
        Mflat = flatten_matrix(Mval, prefixes{j});
        fn = fieldnames(Mflat);
        for f = 1:numel(fn)
            record.(fn{f}) = Mflat.(fn{f});
        end
    end

    % Restoring
    record = flatten_vector(dbg.g, "g", record);

    % ------------------------------------------------------------
    % Outputs
    % ------------------------------------------------------------
    record = flatten_vector(nu_dot,  "nu_dot",  record);
    record = flatten_vector(eta_dot, "eta_dot", record);

    % ------------------------------------------------------------
    % Append
    % ------------------------------------------------------------
    allRecords = [allRecords; struct2table(record)];

    % ------------------------------------------------------------
    % First-case debug print
    % ------------------------------------------------------------
    if k == 1
        disp("==================================================")
        disp("🔍 FIRST CASE DEBUG (RAW, PRE-FLATTEN)")
        disp("==================================================")

        fprintf("nu   = ["); fprintf(' %.15f', x(1:6)); fprintf(" ]\n");
        fprintf("eta  = ["); fprintf(' %.15f', x(7:12)); fprintf(" ]\n");
        fprintf("tau  = ["); fprintf(' %.15f', tau); fprintf(" ]\n\n");

        disp("nu_dot ="); disp(num2str(nu_dot,  '%.15f'));
        disp("eta_dot ="); disp(num2str(eta_dot,'%.15f'));

        disp("==================================================")
    end
end

% ================================================================
% Save CSV
% ================================================================
thisFile  = mfilename('fullpath');
scriptDir = fileparts(thisFile);
outRoot = getenv('OUT_DIR'); if isempty(outRoot), outRoot = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outRoot, 'full');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

outFile = fullfile(outDir, 'spheroid_auv_dynamics_full_debug_template.csv');
writetable(allRecords, outFile);
dataDir = fileparts(scriptDir);
for target = {'hydrodynamics/spheroid_matlab_reference.csv', 'hydrostatics/spheroid_matlab_reference.csv'}
    [d, ~] = fileparts(target{1}); if ~exist(fullfile(outRoot, d), 'dir'), mkdir(fullfile(outRoot, d)); end
    extract_like(outFile, fullfile(dataDir, target{1}), fullfile(outRoot, target{1}));
end

disp("📁 Saved MATLAB results to: " + string(outFile));

%% ================================================================
% Helper: Flatten vector
% ================================================================
function S = flatten_vector(vec, prefix, S)
    if nargin < 3, S = struct(); end
    vec = vec(:);
    for i = 1:numel(vec)
        S.(sprintf('%s_%02d', prefix, i)) = vec(i);
    end
end

%% ================================================================
% Helper: Flatten matrix (row-major)
% ================================================================
function S = flatten_matrix(M, prefix)
    M = M.';   % transpose -> row-major flatten
    M = M(:);
    S = struct();
    for i = 1:numel(M)
        S.(sprintf('%s_%02d', prefix, i)) = M(i);
    end
end
%% remus100
function [xdot,U,M, dbg] = remus100_core(x,tau,Vc,betaVc,w_c)
% REMUS100_CORE  Core hydrodynamic REMUS 100 model (no internal actuation)
% Compatible with MATLAB/Octave.
%
% Simulates passive 6-DOF AUV motion with hydrodynamic forces, 
% restoring forces, damping, and ocean currents.
%
% Inputs:
%   x     : [12×1] or [13×1] state vector
%            x = [u v w p q r x y z phi theta psi]'
%            or quaternion form [u v w p q r x y z eta eps1 eps2 eps3]'
%   tau   : [6×1] generalized control forces/moments [X Y Z K M N]'
%   Vc    : Current velocity magnitude [m/s] (optional)
%   betaVc: Current direction [rad] (optional)
%   w_c   : Vertical current [m/s] (optional)
%
% Outputs:
%   xdot  : State derivative
%   U     : Vehicle total speed [m/s]
%   M     : 6×6 total inertia matrix
%
% Notes:
%   - No rudder, stern-plane, or propeller model inside.
%   - tau is an external input (from actuator or controller).
%
% Author: Adapted from Thor I. Fossen (2025)
% Date:   2025-12-06

if nargin == 0
    x = zeros(12,1); tau = zeros(6,1); Vc = 0; betaVc = 0; w_c = 0;
end
if nargin == 2, Vc = 0; betaVc = 0; w_c = 0; end
if nargin == 4, w_c = 0; end

if length(x) ~= 12 && length(x) ~= 13
    error('x-vector must have dimension 12 or 13');
end
if length(tau) ~= 6
    error('tau-vector must have dimension 6');
end

% -------------------------------------------------------------------------
% Constants
% -------------------------------------------------------------------------
mu  = deg2rad(63.446827);  % Latitude for Trondheim
g_mu = gravity(mu);        % Gravity vector (m/s²)
rho = 1025;                % Water density (kg/m³)

% -------------------------------------------------------------------------
% Kinematic and dynamic variables
% -------------------------------------------------------------------------
nu = x(1:6); 
if length(x) == 12
    psi = x(12);
else
    psi = atan2(2*(x(10)*x(13) + x(11)*x(12)), 1 - 2*(x(12)^2 + x(13)^2));
end

% Ocean currents (BODY frame)
u_c = Vc * cos(betaVc - psi);
v_c = Vc * sin(betaVc - psi);
nu_c = [u_c v_c w_c 0 0 0]';
Dnu_c = [nu(6)*v_c -nu(6)*u_c 0 0 0 0]';

% Relative velocities
nu_r = nu - nu_c;
alpha = atan2(nu_r(3), nu_r(1));
U_r = norm(nu_r(1:3));
U   = norm(nu(1:3));

% -------------------------------------------------------------------------
% Geometry and hydrodynamics
% -------------------------------------------------------------------------
L_auv = 1.6; D_auv = 0.19;
S = 0.7 * L_auv * D_auv;
a = 1.0096 * L_auv/2;
b = 1.0096 * D_auv/2;
r44 = 0.3;                      
r_bG = [0 0 0.02]';
r_bB = [0 0 0]';

Cd = 0.42; 
CD_0 = Cd * pi * b^2 / S;

% -------------------------------------------------------------------------
% Mass, Coriolis, Damping
% -------------------------------------------------------------------------
[MRB,CRB] = spheroid(a,b,nu(4:6),r_bG);
[MA,CA] = imlay611(a,b,nu_r,r44);

% CA-terms in roll, pitch and yaw can destabilize the model if quadratic
% rotational damping is missing. These terms are assumed to be zero
CA(5,3) = 0; CA(3,5) = 0;  % Quadratic velocity terms due to pitching
CA(5,1) = 0; CA(1,5) = 0;  
CA(6,1) = 0; CA(1,6) = 0;  % Munk moment in yaw 
CA(6,2) = 0; CA(2,6) = 0;

M = MRB + MA;
C = CRB + CA;
m = MRB(1,1);
W = m * g_mu;
B = W;

% Damping matrix
T1 = 20; T2 = 20; T6 = 1;
zeta4 = 0.3; zeta5 = 0.8;
D = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bG' r_bB']);
D(1,1) = D(1,1) * exp(-3*U_r);
D(2,2) = D(2,2) * exp(-3*U_r);

% -------------------------------------------------------------------------
% Hydrodynamic forces
% -------------------------------------------------------------------------
tau_liftdragg  = forceLiftDragg(D_auv,S,CD_0,alpha,U_r);
tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r,'cylinder');

% Orientation and restoring forces
if length(x) == 13
    [J,R] = quatern(x(10:13));
else
    [J,R] = eulerang(x(10),x(11),x(12));
end
g = gRvect(W,B,R,r_bG,r_bB);

% -------------------------------------------------------------------------
% State-space dynamics
% -------------------------------------------------------------------------
xdot = [ Dnu_c + M \ (tau + tau_liftdragg + tau_crossflow - C*nu_r - D*nu_r - g)
         J * nu ];
% --- Debug snapshot (extended for initialization and validation) ---
dbg = struct();

% --- State & Kinematics ---
dbg.nu       = nu;          % body velocities
dbg.eta      = x(7:12);     % position/orientation
dbg.nu_c     = nu_c;        % ocean current
dbg.Dnu_c    = Dnu_c;       % derivative of current
dbg.nu_r     = nu_r;        % relative velocity
dbg.alpha    = alpha;       % angle of attack [rad]
dbg.U_r      = U_r;         % relative flow speed [m/s]
dbg.U        = U;           % total body speed [m/s]

% --- Matrices ---
dbg.M_RB = MRB;
dbg.M_A  = MA;
dbg.M    = M;

dbg.C_RB = CRB;  
dbg.C_A  = CA;    
dbg.C    = C;

dbg.D    = D;
dbg.g    = g;

% --- Forces & Moments ---
dbg.tau_total      = tau + tau_liftdragg + tau_crossflow;
dbg.tau_actuator   = tau;
dbg.tau_lift_drag  = tau_liftdragg;
dbg.tau_cross_flow = tau_crossflow;

% --- Added-Mass Derivatives (for validation) ---
if exist('MA', 'var')
    diag_derivs = [MA(1,1) MA(2,2) MA(3,3) MA(4,4) MA(5,5) MA(6,6)];
else
    diag_derivs = nan(1,6);
end

% --- Principal inertias from MRB (3×3 rotational block) ---
% Moment of inertia
Ix = (2/5) * m * b^2;
Iy = (1/5) * m * (a^2 + b^2);
Iz = Iy;

% --- Damping coefficients (for display only) ---
T1 = 20; T2 = 20; T6 = 1;
zeta4 = 0.3; zeta5 = 0.8;
damping_coeff = [ ...
    M(1,1)/T1;
    M(2,2)/T2;
    M(3,3)/T2;  % T3 = T2
    M(4,4)*2*zeta4*sqrt(W*(r_bG(3)-r_bB(3))/M(4,4));
    M(5,5)*2*zeta5*sqrt(W*(r_bG(3)-r_bB(3))/M(5,5));
    M(6,6)/T6 ];

% --- Parameters (expanded for printing) ---
dbg.params = struct( ...
    'L', L_auv, ...
    'D', D_auv, ...
    'S', S, ...
    'a', a, ...
    'b', b, ...
    'r44', r44, ...
    'r_bG', r_bG, ...
    'r_bB', r_bB, ...
    'rho', rho, ...
    'm', MRB(1,1), ...
    'Volume', MRB(1,1)/rho, ...
    'W', W, ...
    'B', B, ...
    'g_mu', g_mu, ...
    'Cd', Cd, ...
    'CD_0', CD_0, ...
    'Ix', Ix, ...
    'Iy', Iy, ...
    'Iz', Iz, ...
    'diag_derivs', diag_derivs, ...
    'damping_coeff', damping_coeff, ...
    'zeta4', zeta4, ...
    'zeta5', zeta5, ...
    'T1', T1, ...
    'T2', T2, ...
    'T6', T6 ...
);

end

%%
function [MA,CA] = imlay611(a,b,nu,r44)
% [MA,CA] = imlay61(a,b,nu,r44) computes the 6x6 hydrodynamic added mass  
% system matrix MA and the 6x6 added mass Coriolis and centripetal matrix  
% CA for a prolate spheroid with semiaxes a > b using the Lamb's 
% k-factors k1, k2 and k_prime (Fossen 2021, Section 8.4.2). The matrix MA 
% is assumed to  be diagonal when the CO is chosen on the centerline 
% midtships. The length of the AUV is L = 2*a while the diamater is D = 2*b. 
%
% Inputs: a, b: spheroid semiaxes a > b
%         nu = [u, v, w, p, q, r]': generalized velocity vector
%         r44: hydrodynamic added moment MA(4,4) = r44 * Ix in roll.
%              If r44 is not specified, MA(4,4) = 0.
%              Typicaly values for r44 are 0.2-0.4.
%
% Output: MA: 6x6 diagonal hydrodynamic added mass system matrix
%         CA: 6x6 hydrodynamic added Coriolis and centripetal matrix
%
% Example: [MA,CA] = imlay61(a, b, [u,v,w,p,q,r]')
%          [MA,CA] = imlay61(a, b, [u,v,w,p,q,r]', r44)
%
% Refs: Lamb, H. (1932). Hydrodynamics. Cambridge University Press. London.
%       Imlay, F. H. (1961). The Complete Expressions for Added Mass of a 
%          Rigid Body Moving in an Ideal Fluid. Technical Report DTMB 1528. 
%          David Taylor Model Basin. Washington D.C.
%
% Author:     Thor I. Fossen 
% Date:       24 Apr 2021
% Revisions:  
     
% prolate spheroid formulas
rho = 1025;
m = 4/3 * pi * rho * a * b^2;
Ix = (2/5) * m * b^2;
Iy = (1/5) * m * (a^2 + b^2);      

% Imlay (1961) gives a zero added moment in roll for a spheroid. This is
% compensated by manually specifying a nonzero r44 for other hull effects.  
% MA(4,4) will in practise be nonzero due to control surfaces, propellers etc. 

if (a < 0), error('a must be larger than 0'); end
if (b < 0), error('b must be larger than 0'); end
if (a <= b), error('a must be larger than b'); end
    
if (nargin == 3)
   r44 = 0;
end
MA_44 = r44 * Ix; 

% Lamb's k-factors
e = sqrt(1-(b/a)^2);
alpha_0 = ( 2 * (1-e^2)/e^3 ) * ( 0.5 * log((1+e)/(1-e)) - e );  
beta_0  = 1/e^2 - (1-e^2)/(2*e^3) * log((1+e)/(1-e)); 

k1 = alpha_0 / (2 - alpha_0);
k2 = beta_0  / (2 - beta_0);
k_prime = e^4*(beta_0-alpha_0) / ((2-e^2)*(2*e^2-(2-e^2)*(beta_0-alpha_0)));    

% Added mass system matrix expressed in the CO
MA = diag([m*k1 m*k2 m*k2 MA_44 k_prime*Iy k_prime*Iy]);

% Added mass Coriolis and centripetal matrix expressed in the CO
CA = m2c(MA,nu);
end

%%
function tau_liftdragg = forceLiftDragg(b,S,CD_0,alpha,U_r)
% tau_liftdrag = forceLiftDrag(b,S,CD_0,alpha,Ur) computes the hydrodynamic
% lift and drag forces of a submerged "wing profile" for varying angle of
% attack (Beard and McLain 2012). Application:
%
%  M d/dt nu_r + C(nu_r)*nu_r + D*nu_r + g(eta) = tau + tau_liftdrag
%
% Output:
%  tau_liftdrag:  6x1 generalized force vector
%
% Inputs:
%  b:       wing span (m)
%  S:       wing area (m^2)
%  CD_0:    parasitic drag (alpha = 0), typically 0.1-0.2 for a streamlined body
%  alpha:   angle of attack, scalar or vector (rad)
%  U_r:     relative speed (m/s)
%
% Example:
%
% Cylinder-shaped AUV with length L = 1.8, diameter D = 0.2 and CD_0 = 0.1:
%    tau_liftdrag = forceLiftDrag(0.2, 1.8*0.2, 0.1, alpha, U_r)
% 
% Author:    Thor I. Fossen
% Date:      25 April 2021 

rho = 1025;

[CL,CD] = coeffLiftDrag(b,S,CD_0,alpha,0);

F_drag = 1/2 * rho * U_r^2 * S * CD;    % drag force
F_lift = 1/2 * rho * U_r^2 * S * CL;    % lift force

% transform from FLOW axes to BODY axes using angle of attack
% transform from FLOW axes to BODY axes using angle of attack
tau_liftdragg = [...
    cos(alpha) * (-F_drag) - sin(alpha) * (-F_lift)
    0
    sin(alpha) * (-F_drag) + cos(alpha) * (-F_lift)
    0
    0
    0 ];
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
