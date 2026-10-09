%% ================================================================
%  Gravity reference, computed by MSS INS/functions/gravity.m
%
%  WGS-84 latitude gravity at the one latitude the submerged-hydrostatics
%  reference uses (remus100.m line 96, 63.446827 deg, Trondheim), so the
%  block's ``more_transformations`` gravity has a frozen MATLAB oracle
%  (Transforms rule, E-62: "a block's oracle is a frozen MSS/MATLAB output
%  ..., not a Python re-implementation of the transform") instead of a
%  Python transcription of gravity.m's coefficients.
%
%  Output, under <OUT_DIR>/restoring/gravity_mss_reference.csv:
%    mu_deg,g
%  Run (nothing relative to a machine): MSS_DIR (required) = the MSS
%  checkout, MATLAB_BIN = the MATLAB executable, OUT_DIR (optional) =
%  output folder, default <tempdir>/more_mss_references; never this folder,
%  so a run cannot overwrite a frozen file. Byte-compare with the frozen
%  file of the same name here:
%    MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this file>')"
%
%  Author:    Enio Krizman
%  Date:      2026-10-09
% ================================================================
clear functions;
format long g
format compact

mssDir = getenv('MSS_DIR');   % MSS is found only through MSS_DIR, never a sibling folder
if isempty(mssDir)
    error('MSS_DIR is not set: point it at the MSS checkout root.');
end
addpath(genpath(mssDir));
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
if ~exist(fullfile(outDir, 'restoring'), 'dir'), mkdir(fullfile(outDir, 'restoring')); end
fprintf('MATLAB %s\nMSS: %s\n', version, which('gravity'));

mu_deg = 63.446827;   % remus100.m line 96 (Trondheim, Norway)
g = gravity(deg2rad(mu_deg));

out = fullfile(outDir, 'restoring', 'gravity_mss_reference.csv');
fid = fopen(out, 'w');
fprintf(fid, 'mu_deg,g\n%.17g,%.17g\n', mu_deg, g);
fclose(fid);
fprintf('Wrote %s\n', out);
