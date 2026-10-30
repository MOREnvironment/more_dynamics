% generate_grethe_mss.m  Frozen MATLAB references of the Grethe model file
% grethe.m (MSS style, CRAFT/USV/models/grethe.m), for the G1 and G3 gates of
% the library's Monohull. The model file and its runs are not in this
% repository: they are read from the MSS checkout that holds them.
%
%  Outputs (written under OUT_DIR, default <tempdir>/more_grethe_references):
%    grethe_m_derivative.csv         G1: 300 seeded states and commands,
%                                    xdot of grethe.m (thrust vector only),
%                                    the outboard wrench tau_out, and the
%                                    same with the leg force, with the
%                                    restoring matrix given at x = -2.24 m,
%                                    and with the smoothstep residual (ramp)
%    grethe_m_mass_matrix.csv        [~,~,M] = grethe()
%    grethe_m_payload_derivative.csv G1: 12 seeded states with a payload of
%                                    25 kg at a seeded position
%    grethe_m_outboard.json          the outboard limits and position of
%                                    grethe.m (read, not typed)
%    grethe_m_trajectory_<name>.csv  G3 open-loop runs of exGrethe.m (the
%                                    'gate' set), and with suffixes _leg and
%                                    _x_ref the two other sets
%    ..._x_ref.csv, xdot_x_ref_*     the restoring matrix given at x = -2.24 m
%                                    is computed here by a marked copy of
%                                    grethe.m (see 'The shift to the point P'
%                                    below); the output of the unmodified file
%                                    is kept beside it: xdot_x_ref_as_written_*
%                                    columns and *_x_ref_as_written.csv files
%    grethe_m_trajectory_tolerances*.csv  per run and state: the largest
%                                    difference to the run with every
%                                    derivative shifted by +-1e-9 and to the
%                                    RK4 with its sum reordered
%    grethe_m_plant.json             the actuator values of grethe.m's lim
%                                    (thrust and steering time constants,
%                                    steering rate and angle limits)
%    grethe_m_plant_derivative.csv   G1: 300 seeded states of the plant of
%                                    the model file with its defaults (leg
%                                    force, ramp, restoring at the origin) and
%                                    the actuator equations of its caller,
%                                    SIMgrethe.m 257-259 in continuous form,
%                                    xadot = [Th_dot; delta_dot; xdot]
%    grethe_m_plant_trajectory_<name>.csv  G3: the six runs of exGrethe.m on
%                                    that plant, thrust and steering commands
%                                    held 10 steps, actuators at rest at the
%                                    start; columns t, ui_01..02 (commands),
%                                    x_01..12, act_01..02 (thrust, steering)
%    grethe_m_plant_trajectory_tolerances.csv  per run and state, as above
%
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS checkout
%  with grethe.m, OUT_DIR (optional) = output folder:
%    MSS_DIR=<MSS checkout> OUT_DIR=<folder> "$MATLAB_BIN" -batch "run('<this file>')"
%  with TMPDIR unset. The seed of the cases is 20261023.
%
% The shift to the point P. grethe.m moves the restoring matrix from the CO to the
% point P with H_ref' G_CO H_ref and H_ref = Hmtrx([opt.x_ref 0 0]) (line 349),
% the form MSS Gmtrx.m carried up to ed2ce7b. The elastic energy of a small
% rigid motion is the same at P and at the CO, so the matrix at P is
% H(-r)' G_CO H(-r), the inverse shift (MSS Gmtrx.m 40-41 since 6e78e98).
% The x_ref cases are therefore computed by a copy of grethe.m written at run
% time into a folder that is put first on the path for those cases only, with
% that one line changed and marked; the unmodified file's output is written
% beside it as history. The other cases (x_ref = 0) are the unmodified file's.
%
% Author:    Enio Krizman
% Date:      2026-10-23
% Revisions: none by others

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

% --- grethe.m with the inverse CO -> P shift (marked copy, written at run time) ---
srcText = fileread(which('grethe'));
oldLine = 'H_ref = Hmtrx([opt.x_ref 0 0]);';
assert(contains(srcText, oldLine), 'grethe.m: the H_ref line was not found; revisit the marked copy');
newLine = ['H_ref = Hmtrx([-opt.x_ref 0 0]);   % MARKED: inverse shift, ' ...
           'H(r)^-1 = H(-r); the file as written has Hmtrx([opt.x_ref 0 0])'];
corrDir = fullfile(tempdir, 'grethe_p_corrected');
if ~exist(corrDir,'dir'), mkdir(corrDir); end
% the copy keeps the function name grethe so that exGrethe.m calls it; it is the first
% on the path only while the x_ref cases are computed (grethe_p, the helper at the end)
fidc = fopen(fullfile(corrDir,'grethe.m'),'w'); fwrite(fidc, strrep(srcText, oldLine, newLine)); fclose(fidc);

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
opt_gate = struct('leg',false,'x_ref',0,'ramp',false);   % ramp = false: residual zero below Fn_V 0.6
opt_leg  = struct('leg',true,'x_ref',0,'ramp',false);
opt_xref = struct('leg',false,'x_ref',-2.24,'ramp',false);
opt_ramp = struct('leg',false,'x_ref',0,'ramp',true);
rows = zeros(N, 1 + 12 + 2 + 6 + 12 + 6 + 12 + 12 + 12 + 12);
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
        ui = [ -511.5 + 1023*rand; -0.6 + 1.2*rand ];           % thrust (N), angle (rad)
    end
    [xd, ~, ~, ~, tau_out]   = grethe(x,ui,0,zeros(3,1),0,0,opt_gate);
    [xd_leg, ~, ~, ~, tau_leg] = grethe(x,ui,0,zeros(3,1),0,0,opt_leg);
    xd_xref_as_written = grethe(x,ui,0,zeros(3,1),0,0,opt_xref);
    xd_xref = grethe_p(x,ui,0,zeros(3,1),0,0,opt_xref);
    xd_ramp = grethe(x,ui,0,zeros(3,1),0,0,opt_ramp);
    rows(k,:) = [k x' ui' tau_out' xd' tau_leg' xd_leg' xd_xref' xd_ramp' xd_xref_as_written'];
end
names = ['case_id', arrayfun(@(i) sprintf('x_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('ui_%02d',i),1:2,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('tau_out_%02d',i),1:6,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('tau_out_leg_%02d',i),1:6,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_leg_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_x_ref_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_ramp_%02d',i),1:12,'UniformOutput',false), ...
         arrayfun(@(i) sprintf('xdot_x_ref_as_written_%02d',i),1:12,'UniformOutput',false)];
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
        ui = [ -511.5 + 1023*rand; -0.6 + 1.2*rand ];
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
exGrethe(outDir,'ramp');
% x_ref: the unmodified file first (history, renamed), then the marked copy
asWrittenDir = fullfile(tempdir, 'grethe_x_ref_as_written');
if exist(asWrittenDir,'dir'), rmdir(asWrittenDir,'s'); end
mkdir(asWrittenDir);
exGrethe(asWrittenDir,'x_ref');
listing = dir(fullfile(asWrittenDir,'*_x_ref.csv'));
for q = 1:numel(listing)
    [~,base] = fileparts(listing(q).name);
    copyfile(fullfile(asWrittenDir,listing(q).name), fullfile(outDir,[base '_as_written.csv']));
end
rmdir(asWrittenDir,'s');
addpath(corrDir,'-begin'); clear functions;
assert(strcmp(fileparts(which('grethe')), corrDir), 'the marked copy of grethe.m is not first on the path');
exGrethe(outDir,'x_ref');
rmpath(corrDir); clear functions;

% --- the plant: grethe.m with its defaults and its caller's actuators ------------
% grethe.m takes the actual thrust and steering angle as inputs and leaves their
% lags and limits to the caller (grethe.m 38-48; SIMgrethe.m 151-154, 257-259).
% The plant adds them as two states, in continuous form (SIMgrethe.m takes an
% Euler step of them): thrust T' = (T_c - T) / T_thrust, steering
% delta' = satlim((satlim(delta_c) - delta) / T_delta, +-ddelta_max), the command
% limited to +-delta_max (the caller limits the command, SIMgrethe.m 245).
[~,~,~,lim] = grethe();
fid = fopen(fullfile(outDir,'grethe_m_plant.json'),'w');
fprintf(fid,'%s\n',jsonencode(struct('thrust_time_constant',lim.T_thrust,'steering_time_constant',lim.T_delta, ...
    'steering_rate_limit',lim.ddelta_max,'steering_angle_limit',lim.delta_max),'PrettyPrint',true));
fclose(fid);
rng(20261026,'twister');
Nq = 300;
qrows = zeros(Nq, 1 + 14 + 2 + 14);
for k = 1:Nq
    if k == 1
        xa = zeros(14,1); uc = [0; 0];
    elseif k == 2
        xa = zeros(14,1); uc = [300; 0.3];
    else
        x = [ -1 + 4*rand; -0.5 + rand; -0.3 + 0.6*rand; -0.3 + 0.6*rand; -0.3 + 0.6*rand; -0.3 + 0.6*rand; ...
              -10 + 20*rand; -10 + 20*rand; -0.2 + 0.4*rand; -0.15 + 0.3*rand; -0.15 + 0.3*rand; -pi + 2*pi*rand ];
        xa = [ -300 + 811*rand; -lim.delta_max + 2*lim.delta_max*rand; x ];       % actual thrust (N) and angle (rad), then the state
        uc = [ -600 + 1200*rand; -0.9 + 1.8*rand ];                              % commands, some beyond the limits
    end
    qrows(k,:) = [k xa' uc' grethe_plant(xa, uc, lim)'];
end
qnames = ['case_id', arrayfun(@(i) sprintf('xa_%02d',i),1:14,'UniformOutput',false), ...
          arrayfun(@(i) sprintf('uc_%02d',i),1:2,'UniformOutput',false), ...
          arrayfun(@(i) sprintf('xadot_%02d',i),1:14,'UniformOutput',false)];
fid = fopen(fullfile(outDir,'grethe_m_plant_derivative.csv'),'w');
fprintf(fid,'%s\n',strjoin(qnames,','));
fprintf(fid,[repmat('%.17g,',1,size(qrows,2)-1) '%.17g\n'],qrows');
fclose(fid);
fprintf('%d plant derivative cases written\n',Nq);

% the six runs of exGrethe.m (its lines 62-71 and the hold of 10 steps at h = 0.05 s) on the plant
h = 0.05; hold_steps = 10; OFFSET = 1e-9; d2r = pi/180;
x0 = zeros(12,1); x0_fast = zeros(12,1); x0_fast(1) = 2.40;
x0_roll = zeros(12,1); x0_roll(10) = 5*d2r;
x0_pitch = zeros(12,1); x0_pitch(11) = 3*d2r;
scen = {
 'straight_from_rest_thrust_step', x0,       @(t) [400; 0],                                       150
 'coast_down_from_2p4',            x0_fast,  @(t) [0; 0],                                         100
 'turning_circle_20deg',           x0_fast,  @(t) [450; 20*d2r],                                  120
 'steering_square_wave_5deg',      x0_fast,  @(t) [450; 5*d2r*sign(sin(2*pi*t/20)+1e-9)],         120
 'roll_free_decay',                x0_roll,  @(t) [0; 0],                                          30
 'pitch_free_decay',               x0_pitch, @(t) [0; 0],                                          30 };
tolRows = {};
for k = 1:size(scen,1)
    name = scen{k,1}; rule = scen{k,3}; tEnd = scen{k,4};
    nRows = round(tEnd/(h*hold_steps)) + 1;
    xa = [0; 0; scen{k,2}];
    times = (0:nRows-1)' * h * hold_steps;
    U = zeros(nRows,2); XA = zeros(nRows,14);
    for r = 1:nRows
        XA(r,:) = xa';
        uc = rule(times(r)); U(r,:) = uc';
        if r == nRows, break; end
        for j = 1:hold_steps
            xa = rk4(@(xx,uu) grethe_plant(xx,uu,lim), h, xa, uc);
        end
    end
    dev = zeros(3,14);
    variants = {+OFFSET, -OFFSET, 'reordered'};
    for v = 1:3
        xb = [0; 0; scen{k,2}]; dmax = zeros(1,14);
        for r = 1:nRows
            dmax = max(dmax, abs(xb' - XA(r,:)));
            if r == nRows, break; end
            for j = 1:hold_steps
                if ischar(variants{v})
                    k1 = grethe_plant(xb, U(r,:)', lim); k2 = grethe_plant(xb + 0.5*h*k1, U(r,:)', lim);
                    k3 = grethe_plant(xb + 0.5*h*k2, U(r,:)', lim); k4 = grethe_plant(xb + h*k3, U(r,:)', lim);
                    xb = xb + h*(k1/6 + k2/3 + k3/3 + k4/6);
                else
                    xb = rk4(@(xx,uu) grethe_plant(xx,uu,lim) + variants{v}, h, xb, U(r,:)');
                end
            end
        end
        dev(v,:) = dmax;
    end
    % states in the order of the run files: x_01..12 (MSS order), then the two actuator states
    order = [3:14 1 2];
    for jj = 1:14
        tolRows(end+1,:) = {name, jj, dev(1,order(jj)), dev(2,order(jj)), dev(3,order(jj)), max(dev(:,order(jj)))}; %#ok<AGROW>
    end
    rnames = ['t', arrayfun(@(i) sprintf('ui_%02d',i),1:2,'UniformOutput',false), ...
              arrayfun(@(i) sprintf('x_%02d',i),1:12,'UniformOutput',false), ...
              arrayfun(@(i) sprintf('act_%02d',i),1:2,'UniformOutput',false)];
    fid = fopen(fullfile(outDir,['grethe_m_plant_trajectory_' name '.csv']),'w');
    fprintf(fid,'%s\n',strjoin(rnames,','));
    fprintf(fid,[repmat('%.17g,',1,size(rnames,2)-1) '%.17g\n'],[times U XA(:,3:14) XA(:,1:2)]');
    fclose(fid);
    fprintf('plant %-32s %5d rows, end u = %.4f m/s, psi = %.4f rad\n',name,nRows,XA(end,3),XA(end,14));
end
fid = fopen(fullfile(outDir,'grethe_m_plant_trajectory_tolerances.csv'),'w');
fprintf(fid,'scenario,state,offset_plus,offset_minus,reordered,tolerance\n');
for i = 1:size(tolRows,1), fprintf(fid,'%s,%d,%.17g,%.17g,%.17g,%.17g\n',tolRows{i,:}); end
fclose(fid);
fprintf('Saved to %s\n',outDir);

function xadot = grethe_plant(xa, uc, lim)
% The plant: the model file with its defaults and its caller's actuator equations.
    Th = xa(1); delta = xa(2); x = xa(3:14);
    xd = grethe(x, [Th; delta], 0, zeros(3,1), 0, 0);
    Th_dot = (uc(1) - Th) / lim.T_thrust;
    d_dot = satlim((satlim(uc(2), -lim.delta_max, lim.delta_max) - delta) / lim.T_delta, ...
        -lim.ddelta_max, lim.ddelta_max);
    xadot = [Th_dot; d_dot; xd];
end

function varargout = grethe_p(varargin)
% grethe.m with the inverse CO -> P shift (the marked copy of the folder corrDir), called
% for the x_ref derivative cases; the unmodified file is restored on the path afterwards.
    corrDir = fullfile(tempdir, 'grethe_p_corrected');
    addpath(corrDir,'-begin'); clear grethe;
    cleanup = onCleanup(@() cleanup_path(corrDir));
    [varargout{1:max(1,nargout)}] = grethe(varargin{:});
end

function cleanup_path(corrDir)
    rmpath(corrDir); clear grethe;
end
