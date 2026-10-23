% generate_grethe_mss.m  Frozen MATLAB references of the Grethe model file
% grethe.m (MSS style, CRAFT/USV/models/grethe.m), for the G1 and G3 gates of
% the library's Monohull. The model file and its runs are not in this
% repository: they are read from the MSS checkout that holds them.
%
%  Outputs (written under OUT_DIR, default <tempdir>/more_grethe_references):
%    grethe_m_derivative.csv         G1: 300 seeded states and commands,
%                                    xdot of grethe.m (thrust vector only),
%                                    the outboard wrench tau_out, and the
%                                    same with the leg force and with the
%                                    restoring matrix given at x = -2.24 m
%    grethe_m_mass_matrix.csv        [~,~,M] = grethe()
%    grethe_m_payload_derivative.csv G1: 12 seeded states with a payload of
%                                    25 kg at a seeded position
%    grethe_m_outboard.json          the outboard limits and position of
%                                    grethe.m (read, not typed)
%    grethe_m_trajectory_<name>.csv  G3 open-loop runs of exGrethe.m (the
%                                    'gate' set), and with suffixes _leg and
%                                    _x_ref the two other sets
%    grethe_m_trajectory_tolerances*.csv  per run and state: the largest
%                                    difference to the run with every
%                                    derivative shifted by +-1e-9 and to the
%                                    RK4 with its sum reordered
%
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS checkout
%  with grethe.m, OUT_DIR (optional) = output folder:
%    MSS_DIR=<MSS checkout> OUT_DIR=<folder> "$MATLAB_BIN" -batch "run('<this file>')"
%  with TMPDIR unset. The seed of the cases is 20261023.
%
% Author:    Enio Krizman
% Date:      2026-10-23

mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout that holds CRAFT/USV/models/grethe.m.');
end
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_grethe_references'); end
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(genpath(mssDir));
[status, rev] = system(['git -C "' mssDir '" rev-parse --short HEAD']);
if status ~= 0, rev = 'unknown'; end
fprintf('MATLAB %s\nMSS: %s at %s', version, which('grethe'), rev);

% --- mass matrix -------------------------------------------------------------
[~,~,M] = grethe();
fid = fopen(fullfile(outDir,'grethe_m_mass_matrix.csv'),'w');
fprintf(fid,'%s\n',strjoin(arrayfun(@(j) sprintf('M_%d',j),1:6,'UniformOutput',false),','));
fprintf(fid,'%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n',M');
fclose(fid);

% --- outboard values of the model file --------------------------------------
[~,~,~,lim] = grethe();
fid = fopen(fullfile(outDir,'grethe_m_outboard.json'),'w');
fprintf(fid,'%s\n',jsonencode(struct('max_thrust',lim.T_max,'max_propulsive_power',lim.P_prop, ...
    'position',lim.r_P(:)','steering_convention',lim.steering_convention),'PrettyPrint',true));
fclose(fid);

% --- derivative cases (G1) ---------------------------------------------------
rng(20261023,'twister');
N = 300;
opt_gate = struct('leg',false,'x_ref',0);
opt_leg  = struct('leg',true,'x_ref',0);
opt_xref = struct('leg',false,'x_ref',-2.24);
rows = zeros(N, 1 + 12 + 2 + 6 + 12 + 6 + 12 + 12);
for k = 1:N
    if k == 1
        x = zeros(12,1); ui = [0; 0];                            % at rest
    elseif k == 2
        x = zeros(12,1); ui = [300; 0.3];                        % thrust at rest
    else
        x = [ -1 + 4*rand            % u in [-1,3] m/s
              -0.5 + rand            % v
              -0.3 + 0.6*rand        % w
              -0.3 + 0.6*rand        % p
              -0.3 + 0.6*rand        % q
              -0.3 + 0.6*rand        % r
              -10 + 20*rand          % x
              -10 + 20*rand          % y
              -0.2 + 0.4*rand        % z (m)
              -0.15 + 0.3*rand       % phi (rad)
              -0.15 + 0.3*rand       % theta (rad)
              -pi + 2*pi*rand ];     % psi
        ui = [ -511.5 + 1023*rand; -0.7 + 1.4*rand ];           % thrust (N), angle (rad)
    end
    [xd, ~, ~, ~, tau_out]   = grethe(x,ui,0,zeros(3,1),0,0,opt_gate);
    [xd_leg, ~, ~, ~, tau_leg] = grethe(x,ui,0,zeros(3,1),0,0,opt_leg);
    xd_xref = grethe(x,ui,0,zeros(3,1),0,0,opt_xref);
    rows(k,:) = [k x' ui' tau_out' xd' tau_leg' xd_leg' xd_xref'];
end
names = ['case_id', arrayfun(@(i) sprintf('x_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('ui_%02d',i),1:2,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('tau_out_%02d',i),1:6,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('tau_out_leg_%02d',i),1:6,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_leg_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_x_ref_%02d',i),1:12,'UniformOutput',false)];
fid = fopen(fullfile(outDir,'grethe_m_derivative.csv'),'w');
fprintf(fid,'%s\n',strjoin(names,','));
fprintf(fid,[repmat('%.17g,',1,size(rows,2)-1) '%.17g\n'],rows');
fclose(fid);
fprintf('%d derivative cases written\n',N);

% --- payload cases (G1) --------------------------------------------------------
rng(20261024,'twister');
Np = 12;
prow = zeros(Np, 1 + 4 + 12 + 2 + 12);
for k = 1:Np
    mP = 25; r_bP = [-0.5 + rand; -0.3 + 0.6*rand; -0.3 + 0.3*rand];
    if k == 1, x = zeros(12,1); ui = [0;0]; else
        x = [ -1 + 4*rand; -0.5 + rand; -0.3 + 0.6*rand; -0.3 + 0.6*rand; -0.3 + 0.6*rand; -0.3 + 0.6*rand; ...
              -10 + 20*rand; -10 + 20*rand; -0.2 + 0.4*rand; -0.15 + 0.3*rand; -0.15 + 0.3*rand; -pi + 2*pi*rand ];
        ui = [ -511.5 + 1023*rand; -0.7 + 1.4*rand ];
    end
    xd = grethe(x,ui,mP,r_bP,0,0,opt_gate);
    prow(k,:) = [k mP r_bP' x' ui' xd'];
end
names = ['case_id','payload_mass','payload_position_01','payload_position_02','payload_position_03', ...
         arrayfun(@(i) sprintf('x_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('ui_%02d',i),1:2,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_%02d',i),1:12,'UniformOutput',false)];
fid = fopen(fullfile(outDir,'grethe_m_payload_derivative.csv'),'w');
fprintf(fid,'%s\n',strjoin(names,','));
fprintf(fid,[repmat('%.17g,',1,size(prow,2)-1) '%.17g\n'],prow');
fclose(fid);

% --- open-loop runs (G3) -----------------------------------------------------
exGrethe(outDir,'gate');
exGrethe(outDir,'leg');
exGrethe(outDir,'x_ref');
fprintf('Saved to %s\n',outDir);
