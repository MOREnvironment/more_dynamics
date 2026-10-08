%% ================================================================
%  Wageningen B-series K_T, K_Q, computed by MSS wageningen.m
%
%  MSS LIBRARY/modeling/wageningen.m (T. I. Fossen, MIT) is called as
%  written, on a J grid -0.5 ... 2.0 (past both ends of the regression's
%  fitted range, because the block's default follows MSS at any J) plus
%  seeded random J, and on three propeller geometries (test constructions):
%    remus_like  P/D 1,    AE/AO 0.718, z 3  (remus100.m line 157 at ac77394)
%    four_blade  P/D 0.83, AE/AO 0.55,  z 4
%    outboard    P/D 0.83, AE/AO 0.718, z 3
%  First made at MSS ac77394 (2026-10-06); byte-identical at 72656d1 and
%  cc07579.
%  Output: <OUT_DIR>/force_producers/wageningen_kt_kq_mss_current.csv
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
% MSS is found only through MSS_DIR (no path relative to this repository).
% Run with, e.g.,
%   MSS_DIR=/path/to/MSS "$MATLAB_BIN" -batch "run('<this file>')"
mssDir = getenv('MSS_DIR');
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout root.');
end
addpath(genpath(mssDir));
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
for sub = {'force_producers'}
    if ~exist(fullfile(outDir, sub{1}), 'dir'), mkdir(fullfile(outDir, sub{1})); end
end
fprintf('MATLAB %s\nMSS: %s\n', version, which('wageningen'));

SEED = 20261006;
rng(SEED, 'twister');
J_grid = [linspace(-0.5, 2.0, 251), 0.6632, -0.5 + 2.5 * rand(1, 1000)];
sets = [1.0 0.718 3; 0.83 0.55 4; 0.83 0.718 3];

rows = zeros(size(sets, 1) * numel(J_grid), 6);
r = 0;
for s = 1:size(sets, 1)
    for J = J_grid
        [KT, KQ] = wageningen(J, sets(s, 1), sets(s, 2), sets(s, 3));
        r = r + 1;
        rows(r, :) = [sets(s, :) J KT KQ];
    end
end

T = array2table(rows, 'VariableNames', {'PD', 'AEAO', 'z', 'J', 'KT', 'KQ'});
outFile = fullfile(outDir, 'force_producers', 'wageningen_kt_kq_mss_current.csv');
writetable(T, outFile);
fprintf('Saved %d rows x %d columns to %s\n', height(T), width(T), outFile);
