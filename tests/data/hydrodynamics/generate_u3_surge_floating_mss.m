%% ================================================================
%  Surge damping, ITTC and floating-damping reference data from MSS itself
%  (2026-10-06, MSS ac77394)
%
%  MSS 2.0.2 (72656d1, 2026-10-07) deleted LIBRARY/modeling/forceSurgeDamping.m
%  (commit 108ceda, 2026-10-06), so this file runs only on MSS ac77394 or
%  older. Re-run on 2026-10-07 with MSS 72656d1 on the path and
%  forceSurgeDamping.m taken from ac77394: xuu_ittc and floating_damping
%  byte-identical (XuuITTC.m, Dmtrx.m unchanged in 2.0.2).
%
%  Calls the MSS functions directly, unmodified, on seeded inputs:
%    LIBRARY/modeling/forceSurgeDamping.m  (both branches: ITTC, thrust_max)
%    LIBRARY/modeling/addedMassSurge.m     (A11, through forceSurgeDamping)
%    LIBRARY/modeling/XuuITTC.m            (bounded ITTC-1957 coefficient)
%    LIBRARY/modeling/Dmtrx.m              (surface-craft branch, G matrix)
%  Every number is written with %.17g (exact double round trip).
%
%  Parameter sets: set 1 = osv.m's own vessel numbers (CRAFT/SHIP/models/osv.m
%  lines 62-80, 129), the only MSS vehicle that calls forceSurgeDamping;
%  sets 2-5 = seeded random small craft (test construction, ranges below).
%  Rows with probe = 1 sit at the ITTC singularity Rn = 100 (u_r = 1e-4/L)
%  and are for the physical tests only, never for G1.
%
%  Run from this folder (MSS_DIR required, MATLAB_BIN = the MATLAB executable):
%    MSS_DIR=<an MSS checkout> "$MATLAB_BIN" -batch "run('generate_u3_surge_floating_mss.m')"
%  Outputs beside this file:
%    surge_damping_mss_ac77394.csv, xuu_ittc_mss_ac77394.csv,
%    floating_damping_mss_ac77394.csv
% ================================================================
clear functions;
format long g

scriptDir = fileparts(mfilename('fullpath'));
mssDir = getenv('MSS_DIR');   % no path relative to a workspace
if isempty(mssDir)
    error('MSS_DIR is not set: point it at an MSS checkout');
end
addpath(genpath(mssDir));
fprintf('MATLAB %s\nMSS: %s\n', version, which('forceSurgeDamping'));

rng(20261006, 'twister');

%% ---------------- parameter sets ---------------------------------
% set 1: osv.m lines 62-67, 70, 76-80, 129 (re-typed here from the pinned
% lines; the Python test re-reads them from osv.m and checks this row)
L = 83; B = 18; T = 5; rho = 1025; Cb = 0.65;
K_max = [300e3 300e3 420e3 655e3]';
sets(1) = struct('m', rho * Cb * L * B * T, 'S', L * B + 2 * T * B, 'L', L, ...
    'B', B, 'T', T, 'Cb', Cb, 'T1', 100, 'rho', rho, 'u_max', 7.7, ...
    'thrust_max', K_max(3) + K_max(4));
% sets 2-5: seeded small craft (ranges are a test construction)
for s = 2:5
    Ls = 1 + 19 * rand;                 % 1-20 m
    Bs = Ls * (0.15 + 0.25 * rand);     % beam
    Ts = Bs * (0.1 + 0.4 * rand);       % draft
    Cbs = 0.4 + 0.5 * rand;
    rhos = 1000 + 30 * rand;
    ms = rhos * Cbs * Ls * Bs * Ts;
    sets(s) = struct('m', ms, 'S', Ls * Bs + 2 * Ts * Bs, 'L', Ls, 'B', Bs, ...
        'T', Ts, 'Cb', Cbs, 'T1', 5 + 195 * rand, 'rho', rhos, ...
        'u_max', 1 + 9 * rand, 'thrust_max', 0);
    sets(s).thrust_max = 0.5 * rhos * sets(s).S * 0.01 * sets(s).u_max^2 * (0.5 + rand);
end

u_fixed = [-8 -6 -4 -3 -2 -1 -0.5 -0.1 -1e-3 0 1e-3 0.05 0.1 0.5 1 2 3 4 6 8];

%% ---------------- forceSurgeDamping (both branches) ---------------
fid = fopen(fullfile(scriptDir, 'surge_damping_mss_ac77394.csv'), 'w');
fprintf(fid, 'set,branch,probe,m,S,L,T1,rho,u_max,thrust_max,u_r,X,Xuu,Xu,A11\n');
for s = 1:numel(sets)
    p = sets(s);
    u_list = [u_fixed, -8 + 16 * rand(1, 20)];
    u_sing = 100 * 1e-6 / p.L;                      % Rn = 100 exactly
    probes = [u_sing, u_sing * (1 + 1e-6), u_sing * (1 + 1e-3), -u_sing];
    for branch = 0:1
        for k = 1:(numel(u_list) + numel(probes))
            if k <= numel(u_list)
                u = u_list(k); probe = 0;
            else
                u = probes(k - numel(u_list)); probe = 1;
            end
            if branch == 0
                [X, Xuu, Xu] = forceSurgeDamping(0, u, p.m, p.S, p.L, p.T1, p.rho, p.u_max);
            else
                [X, Xuu, Xu] = forceSurgeDamping(0, u, p.m, p.S, p.L, p.T1, p.rho, p.u_max, p.thrust_max);
            end
            A11 = addedMassSurge(p.m, p.L, p.rho);
            fprintf(fid, '%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
                s, branch, probe, p.m, p.S, p.L, p.T1, p.rho, p.u_max, p.thrust_max, u, X, Xuu, Xu, A11);
        end
    end
end
fclose(fid);

%% ---------------- XuuITTC ----------------------------------------
fid = fopen(fullfile(scriptDir, 'xuu_ittc_mss_ac77394.csv'), 'w');
fprintf(fid, 'set,rho,L,B,T,C_B,u_r,Xuu\n');
for s = 1:numel(sets)
    p = sets(s);
    u_list = [u_fixed, -8 + 16 * rand(1, 20)];
    for k = 1:numel(u_list)
        Xuu = XuuITTC(u_list(k), p.rho, p.L, p.B, p.T, p.Cb);
        fprintf(fid, '%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
            s, p.rho, p.L, p.B, p.T, p.Cb, u_list(k), Xuu);
    end
end
fclose(fid);
% the default block coefficient (XuuITTC.m lines 28-30)
if XuuITTC(1.3, 1025, 10, 3, 1) ~= XuuITTC(1.3, 1025, 10, 3, 1, 0.65)
    error('XuuITTC default C_B is not 0.65');
end

%% ---------------- Dmtrx.m, surface-craft branch -------------------
fid = fopen(fullfile(scriptDir, 'floating_damping_mss_ac77394.csv'), 'w');
hdr = {'set'};
names = {'MRB', 'MA', 'G'};
for j = 1:3, for i = 1:36, hdr{end+1} = sprintf('%s_%02d', names{j}, i); end, end %#ok<AGROW>
hdr = [hdr, {'T1', 'T2', 'T6', 'zeta4', 'zeta5'}];
for i = 1:36, hdr{end+1} = sprintf('D_%02d', i); end %#ok<AGROW>
fprintf(fid, '%s\n', strjoin(hdr, ','));
for s = 1:20
    scale = 10^(1 + 3 * rand);                       % 10 kg .. 10 t
    A = randn(6); MRB = scale * (A * A' + 6 * eye(6));
    A = randn(6); MA = 0.3 * scale * (A * A' + 6 * eye(6));
    A = randn(6); G = scale * (A + A');               % off-diagonal coupling
    G(3,3) = scale * (5 + 50 * rand); G(4,4) = scale * (1 + 10 * rand);
    G(5,5) = scale * (5 + 50 * rand);
    T126 = 1 + 200 * rand(1, 3);
    z45 = 0.05 + 0.9 * rand(1, 2);
    D = Dmtrx(T126, z45, MRB, MA, G);
    row = [s, rowmajor(MRB), rowmajor(MA), rowmajor(G), T126, z45, rowmajor(D)];
    fprintf(fid, [strjoin(repmat({'%.17g'}, 1, numel(row)), ','), '\n'], row);
end
fclose(fid);
fprintf('Saved surge_damping, xuu_ittc, floating_damping CSVs to %s\n', scriptDir);

function v = rowmajor(M)
    M = M.';
    v = M(:)';
end
