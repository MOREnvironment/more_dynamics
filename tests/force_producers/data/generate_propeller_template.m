%% ================================================================
%  REMUS 100 propeller template reference (MATLAB): tau, X_prop, K_prop, Ja
%
%  Template generator: the linearised propeller of MSS remus100.m written
%  out (remus_propulsion_reference below, as remus100.m stood in 2026-02),
%  with rho = 1025 (remus100.m uses 1026) on a speed / rpm grid.
%  Produces force_producers/propeller_matlab.csv (byte-identical at MSS
%  cc07579).
%  Inputs: inputs.csv beside this file. Output under OUT_DIR (default
%  <tempdir>/more_mss_references; never this folder, so a run cannot overwrite a
%  frozen file). Run: MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%
%  Author:    Enio Krizman
%  Date:      2026-02-20
% ================================================================

clear; clc;

%% --- Output folder ---
mssDir = getenv('MSS_DIR'); if isempty(mssDir), error('MSS_DIR is not set: point it at the MSS checkout root.'); end; addpath(genpath(mssDir));
outRoot = getenv('OUT_DIR'); if isempty(outRoot), outRoot = fullfile(tempdir, 'more_mss_references'); end
outdir = fullfile(outRoot, 'force_producers');
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% --- Sweep ranges ---
rpm_vals = -1500:100:1500;   % include reverse
U_vals = 0:0.25:2.5;

%% --- Preallocate ---
RPM = [];
U = [];
Ja_all = [];

Xprop = [];
Kprop = [];

tau1 = [];
tau2 = [];
tau3 = [];
tau4 = [];
tau5 = [];
tau6 = [];

%% ============================================================
% ====================== SWEEP ================================
%% ============================================================

for i = 1:length(rpm_vals)
    for j = 1:length(U_vals)

        rpm = rpm_vals(i);
        u = U_vals(j);

        [tau,X,K,Ja] = remus_propulsion_reference(u,rpm);

        RPM(end+1,1) = rpm;
        U(end+1,1) = u;
        Ja_all(end+1,1) = Ja;

        Xprop(end+1,1) = X;
        Kprop(end+1,1) = K;

        tau1(end+1,1) = tau(1);
        tau2(end+1,1) = tau(2);
        tau3(end+1,1) = tau(3);
        tau4(end+1,1) = tau(4);
        tau5(end+1,1) = tau(5);
        tau6(end+1,1) = tau(6);

    end
end

%% ============================================================
% ====================== SAVE CSV =============================
%% ============================================================

T = table(RPM,U,Ja_all,Xprop,Kprop, ...
          tau1,tau2,tau3,tau4,tau5,tau6);

csvfile = fullfile(outdir, "propeller_matlab.csv");
writetable(T,csvfile);

fprintf("CSV saved to: %s\n", csvfile);



function [tau, X_prop, K_prop, Ja] = remus_propulsion_reference(U, rpm)

% ============================================================
% Minimal REMUS100 propulsion model (Fossen)
% Extracted directly from remus100.m
% ============================================================

%% --- Constants from REMUS model ---
rho = 1025;               % Match remus100
n_max = 1525;             % RPM saturation

% Saturate RPM
rpm = max(min(rpm,n_max),-n_max);

% Convert to revolutions per second
n = rpm / 60;

D_prop = 0.14;
t_prop = 0.1;
Va = 0.944 * U;

Ja_max = 0.6632;

KT_0 = 0.4566;
KQ_0 = 0.0700;

KT_max = 0.1798;
KQ_max = 0.0312;

% Compute advance number safely
if abs(n) > 1e-6
    Ja = Va / (abs(n) * D_prop);
else
    Ja = 0;
end

%% --- Thrust model ---
if n > 0
    X_prop = rho * D_prop^4 * ( ...
        KT_0 * abs(n) * n + ...
        (KT_max-KT_0)/Ja_max * (Va/D_prop) * abs(n) );

    K_prop = rho * D_prop^5 * ( ...
        KQ_0 * abs(n) * n + ...
        (KQ_max-KQ_0)/Ja_max * (Va/D_prop) * abs(n) );

else
    X_prop = rho * D_prop^4 * KT_0 * abs(n) * n;
    K_prop = rho * D_prop^5 * KQ_0 * abs(n) * n;
end

% Apply thrust deduction like remus100
X_prop = (1 - t_prop) * X_prop;

% Roll moment scaling
K_prop = K_prop / 10;

tau = zeros(6,1);                                
tau(1) =X_prop;
tau(4) = K_prop;

end