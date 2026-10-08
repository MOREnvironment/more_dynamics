%% ================================================================
%  Otter differential thruster template sweep (MATLAB)
%
%  Template generator: the Otter thruster constants and thrust law written out
%  (as MSS otter.m and SIMotter.m stood in 2026-02), with the latitude
%  gravity gravity(63.446827 deg) in the thrust limits; first part: the
%  constants and single-point checks, second part (after the ===== line): the
%  sweep of requested tau_X, tau_N whose rows are the frozen file.
%  Produces force_producers/differential_thruster_matlab.csv (byte-identical
%  at MSS cc07579).
%  Inputs: inputs.csv beside this file. Output under OUT_DIR (default
%  <tempdir>/more_mss_references; never this folder, so a run cannot overwrite a
%  frozen file). Run: MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%
%  Author:    Enio Krizman
%  Date:      2026-02-19
% ================================================================

clear; clc; close all;
mssDir = getenv('MSS_DIR'); if isempty(mssDir), error('MSS_DIR is not set: point it at the MSS checkout root.'); end; addpath(genpath(mssDir));
outRoot = getenv('OUT_DIR'); if isempty(outRoot), outRoot = fullfile(tempdir, 'more_mss_references'); end
outDir = fullfile(outRoot, 'force_producers'); if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Parameters (same as otter.m)
g      = gravity(deg2rad(63.446827));
y_pont = 0.395;

l1 = -y_pont;                      % left lever arm (m)
l2 =  y_pont;                      % right lever arm (m)

k_pos = 0.02216/2;                 % positive bollard coeff (one thruster)
k_neg = 0.01289/2;                 % negative bollard coeff (one thruster)

n_max =  sqrt((0.5*24.4*g)/k_pos); % max shaft speed (rad/s)
n_min = -sqrt((0.5*13.6*g)/k_neg); % min shaft speed (rad/s)

fprintf('================ Thruster Logic Test (otter.m) ================\n');
fprintf('k_pos = %.8f   k_neg = %.8f\n', k_pos, k_neg);
fprintf('n_max = %.6f rad/s   n_min = %.6f rad/s\n\n', n_max, n_min);

%% Propulsion matrix (B_prop of otter.m)
B_prop = k_pos * [...
        1  1; ...
        y_pont  -y_pont ];

disp('B_prop = ');
disp(B_prop);

%% ========== BASIC THRUSTER LOGIC TEST (n -> Thrust -> tau) ==========
tests = { ...
    struct('name','Symmetric +',     'n',[  10;   10]), ...
    struct('name','Differential +',  'n',[  10;    5]), ...
    struct('name','Symmetric -',     'n',[ -10;  -10]), ...
    struct('name','Mixed (+/-)',     'n',[  10;  -10]), ...
    struct('name','Saturate +',      'n',[1e6; 1e6]), ...
    struct('name','Saturate -',      'n',[-1e6;-1e6]) ...
};

for k = 1:numel(tests)
    
    name = tests{k}.name;
    n_in = tests{k}.n;
    n = n_in;
    Thrust = zeros(2,1);

    % Saturation + thrust law
    for i = 1:2
        if n(i) > n_max
            n(i) = n_max;
        elseif n(i) < n_min
            n(i) = n_min;
        end
        if n(i) > 0
            Thrust(i) = k_pos * n(i) * abs(n(i));
        else
            Thrust(i) = k_neg * n(i) * abs(n(i));
        end
    end

    % Generalized force (6DOF form)
    tau = [ Thrust(1)+Thrust(2); 0; 0; 0; 0; -l1*Thrust(1) - l2*Thrust(2) ];

    fprintf('--- %s ---\n', name);
    fprintf('n_in      = [% .6f, % .6f] rad/s\n', n_in(1), n_in(2));
    fprintf('n_sat     = [% .6f, % .6f] rad/s\n', n(1), n(2));
    fprintf('Thrust    = [% .6f, % .6f] N\n', Thrust(1), Thrust(2));
    fprintf('tau1 (X)  = % .6f N\n', tau(1));
    fprintf('tau6 (N)  = % .6f N*m\n\n', tau(6));
end

%% ========== THRUSTER ALLOCATION TEST (tau -> n_cmd -> Thrust -> tau) ==========
fprintf('================ Thruster Allocation Test (otter.m) ================\n');

% Invert B_prop for allocation (invQR, as SIMotter.m)
Binv = inv(B_prop);

disp('B_prop inverse (Binv):');
disp(Binv);

% Test cases for surge + yaw
alloc_tests = { ...
    struct('name','Forward thrust', 'tau',[103; 0]), ...
    struct('name','Forward thrust satutrate', 'tau',[500; 0]), ...
    struct('name','Reverse thrust', 'tau',[-101; 0]), ...
    struct('name','Reverse thrust saturate', 'tau',[-500; 0]), ...
    struct('name','Yaw only 1',       'tau',[0; -10]), ...
    struct('name','Yaw only1',       'tau',[0; 10]), ...
    struct('name','Yaw only saturate',       'tau',[0; -500]), ...
    struct('name','Combined',       'tau',[80; 3]), ...
    struct('name','Mixed torque',   'tau',[30; -5]) ...
};

for k = 1:numel(alloc_tests)
    name = alloc_tests{k}.name;
    tau_2d = alloc_tests{k}.tau;

    % Allocation (u = Binv * tau)
    u = Binv * tau_2d;

    % Compute commanded n from u = n|n|
    n = sign(u) .* sqrt(abs(u));



    % Compute thrust from n_cmd
    for i = 1:1:2
        if n(i) > n_max           % saturation, physical limits
           n(i) = n_max; 
        elseif n(i) < n_min
           n(i) = n_min; 
       end
        
       if n(i) > 0                          
         Thrust(i) = k_pos * n(i) * abs(n(i));    % positive thrust (N) 
       else
         Thrust(i) = k_neg * n(i) * abs(n(i));    % negative thrust (N) 
       end
    end

    % Map thrust -> tau (check closure)
    tau_check = [Thrust(1)+Thrust(2); -l1*Thrust(1) - l2*Thrust(2)];

    fprintf('--- %s ---\n', name);
    fprintf('tau_2d (input) = [tau_X=%.3f, tau_N=%.3f]\n', tau_2d(1), tau_2d(2));
    fprintf('u (n|n|)       = [% .6f, % .6f]\n', u(1), u(2));
    fprintf('n_cmd          = [% .6f, % .6f] rad/s\n', n(1), n(2));
    fprintf('Thrust         = [% .6f, % .6f] N\n', Thrust(1), Thrust(2));
    fprintf('tau1 (X)       = % .6f N\n', tau_check(1));
    fprintf('tau6 (N·m)     = % .6f N·m\n\n', tau_check(2));
end

fprintf('✅ MATLAB allocation and thrust test completed successfully.\n');

%% ========== EXTENDED CONTINUOUS TEST (sweep) ==========
fprintf('================ Extended Thruster Sweep Test ================\n');

dt = 0.1;          % time step [s]
t_phase = 10;      % seconds per ramp section
t_total = 3 * 4 * t_phase;
t = 0:dt:t_total;

% Ramp profile generator
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

% Phase 1: Surge only
tau_X(1:phase1) = tau_X_ramp(1:phase1);
% Phase 2: Yaw only
tau_N(phase1+1:phase2) = tau_N_ramp(1:phase1);
% Phase 3: Both combined
tau_X(phase2+1:end) = tau_X_ramp(1:length(t)-phase2);
tau_N(phase2+1:end) = tau_N_ramp(1:length(t)-phase2);

% Preallocate
n_left = zeros(size(t));
n_right = zeros(size(t));
thrust_left = zeros(size(t));
thrust_right = zeros(size(t));
tau_X_out = zeros(size(t));
tau_N_out = zeros(size(t));

for i = 1:length(t)
    tau_2d = [tau_X(i); tau_N(i)];
    u = Binv * tau_2d;
    n_cmd = sign(u) .* sqrt(abs(u));

    for k = 1:2
        if n_cmd(k) > n_max
            n_cmd(k) = n_max;
        elseif n_cmd(k) < n_min
            n_cmd(k) = n_min;
        end

        if n_cmd(k) > 0
            Thrust(k) = k_pos * n_cmd(k) * abs(n_cmd(k));
        else
            Thrust(k) = k_neg * n_cmd(k) * abs(n_cmd(k));
        end
    end

    tau_check = [Thrust(1)+Thrust(2); -l1*Thrust(1) - l2*Thrust(2)];

    n_left(i) = n_cmd(1);
    n_right(i) = n_cmd(2);
    thrust_left(i) = Thrust(1);
    thrust_right(i) = Thrust(2);
    tau_X_out(i) = tau_check(1);
    tau_N_out(i) = tau_check(2);
end

% Save results
T = table(t', tau_X', tau_N', n_left', n_right', thrust_left', thrust_right', ...
          tau_X_out', tau_N_out', ...
          'VariableNames', {'time','tau_X_cmd','tau_N_cmd','n_cmd_left','n_cmd_right','thrust_left','thrust_right','tau_X','tau_N'});
writetable(T, fullfile(outDir, 'differential_thruster_matlab.csv'));
fprintf('✅ Results saved to results_thruster_test_matlab.csv\n');

% Plot
figure;
subplot(3,1,1);
plot(t, tau_X, 'b', t, tau_N, 'r');
ylabel('Commanded τ [N, N·m]');
legend('τ_X','τ_N'); grid on;

subplot(3,1,2);
plot(t, n_left, 'b', t, n_right, 'r');
ylabel('n [rad/s]');
legend('Left','Right'); grid on;

subplot(3,1,3);
plot(t, thrust_left, 'b', t, thrust_right, 'r');
ylabel('Thrust [N]');
xlabel('Time [s]');
legend('Left','Right'); grid on;

fprintf('✅ MATLAB continuous sweep test completed successfully.\n');
