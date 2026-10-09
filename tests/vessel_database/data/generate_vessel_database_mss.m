%% ================================================================
%  Vessel-database (BEM) reference, read from MSS's shipped LAUV_marie.mat
%
%  Loads the MSS `vessel` structure exactly as `wamit2vessel.m` 50-94
%  (HYDRO/wamit2vessel.m) and `computeManeuveringModel.m` document it, and
%  dumps the fields the U10 loader (ADR 0004 Sec.2, E-37, E-53) must read:
%  the added-mass matrix at the zero and infinite frequency limits (both
%  stored at the body CG, `hydrodynamic_reference`), the restoring matrix C
%  (constant over frequency for this submerged body), the rigid-body mass
%  matrix MRB, the scalar `main` fields (mass, water density, gravity,
%  displaced volume, CG, CB, radii of gyration k44/k55/k66, GM_T, GM_L,
%  length, beam, draught, block coefficient), the frequency and velocity
%  grids, and the three descriptive strings (hydrodynamic_axes,
%  hydrodynamic_reference, hydrodynamic_source). No computation beyond
%  indexing: the CO-shift (Fossen 2011 eq. 3.24, H(r)) and every refusal
%  are the loader's job (U10b), checked by the Python tests against this
%  frozen dump.
%
%  File read: MSS_DIR/HYDRO/vessels_capytaine/LAUV_marie/LAUV_marie.mat.
%  This is the file MSS 2.0.2 ships (byte-identical to the MSS-Capytaine
%  export, A-39 Sec.3.4): no MSS function is called to build it, only
%  MATLAB's own `load`, so field names and shapes come from the file as
%  committed, not from a transcription.
%
%  Sentinel convention (Fossen, MSS-Capytaine `capytaine_vessel.py` 183-186,
%  239-241, matched by `plotBv.m` 20-22): `vessel.freqs = [0, finite omegas,
%  10]`; the 10 rad/s entry holds Capytaine's omega=infinity solution, not a
%  literal 10 rad/s result.
%
%  Output, under <OUT_DIR>/vessel_database/:
%    vessel_database_lauv_marie_mss_reference.csv  (one row, named columns)
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
outDir = getenv('OUT_DIR');
if isempty(outDir), outDir = fullfile(tempdir, 'more_mss_references'); end
if ~exist(fullfile(outDir, 'vessel_database'), 'dir')
    mkdir(fullfile(outDir, 'vessel_database'));
end

matFile = fullfile(mssDir, 'HYDRO', 'vessels_capytaine', 'LAUV_marie', 'LAUV_marie.mat');
if ~isfile(matFile)
    error('LAUV_marie.mat not found at %s (MSS_DIR wrong or file moved).', matFile);
end
fprintf('MATLAB %s\nFile: %s\n', version, matFile);

loaded = load(matFile);
vessel = loaded.vessel;

% ------------------------------------------------------------------
% Frequency and velocity grids; zero/infinite-frequency column indices
% ------------------------------------------------------------------
freqs = vessel.freqs(:)';
nFreq = numel(freqs);
idxZero = find(freqs == 0, 1, 'first');
idxInf  = find(freqs == 10, 1, 'last');
if isempty(idxZero) || isempty(idxInf) || idxInf ~= nFreq
    error('LAUV_marie.mat does not carry the MSS 0/10 sentinel pair as its first/last freqs entry.');
end
velocities = vessel.velocities(:)';
idxVel0 = find(velocities == 0, 1, 'first');
if isempty(idxVel0)
    error('LAUV_marie.mat has no zero-speed column in vessel.velocities.');
end

A_zero = vessel.A(:,:,idxZero,idxVel0);
A_inf  = vessel.A(:,:,idxInf,idxVel0);
B_mid  = vessel.B(:,:,:,idxVel0);     % every frequency, zero speed (raw, unclipped)
C_const = vessel.C(:,:,idxZero,idxVel0);   % constant over frequency (checked below)
MRB = vessel.MRB;

% Sanity: C must not vary with frequency for this submerged body (A-35 Sec.3,
% item 5); if it ever does, the frozen dump must say so rather than silently
% picking one slice.
cVariesOverFreq = false;
for k = 1:nFreq
    if max(abs(vessel.C(:,:,k,idxVel0) - C_const), [], 'all') > 1e-12
        cVariesOverFreq = true;
        break
    end
end

m = vessel.main;

% ------------------------------------------------------------------
% Write one wide CSV row: header line, value line. Numbers at %.17g
% (round-trip double precision); strings quoted.
% ------------------------------------------------------------------
cols = {};
vals = {};

cols{end+1} = 'n_freq';            vals{end+1} = sprintf('%d', nFreq);
cols{end+1} = 'freq_zero_idx';     vals{end+1} = sprintf('%d', idxZero);
cols{end+1} = 'freq_inf_idx';      vals{end+1} = sprintf('%d', idxInf);
cols{end+1} = 'freq_inf_value';    vals{end+1} = sprintf('%.17g', freqs(idxInf));
cols{end+1} = 'velocity_zero_idx'; vals{end+1} = sprintf('%d', idxVel0);
cols{end+1} = 'c_varies_over_freq'; vals{end+1} = sprintf('%d', cVariesOverFreq);

cols{end+1} = 'hydrodynamic_axes';      vals{end+1} = sprintf('"%s"', vessel.hydrodynamic_axes);
cols{end+1} = 'hydrodynamic_reference'; vals{end+1} = sprintf('"%s"', vessel.hydrodynamic_reference);
cols{end+1} = 'hydrodynamic_source';    vals{end+1} = sprintf('"%s"', vessel.hydrodynamic_source);
cols{end+1} = 'name';                   vals{end+1} = sprintf('"%s"', m.name);

matrices = {A_zero, 'A_zero'; A_inf, 'A_inf'; C_const, 'C'; MRB, 'MRB'};
for r = 1:size(matrices,1)
    M = matrices{r,1};
    prefix = matrices{r,2};
    k = 0;
    for i = 1:6
        for j = 1:6
            k = k + 1;
            cols{end+1} = sprintf('%s_%02d', prefix, k); %#ok<AGROW>
            vals{end+1} = sprintf('%.17g', M(i,j));      %#ok<AGROW>
        end
    end
end

scalarFields = {'m','rho','g','nabla','k44','k55','k66','GM_T','GM_L', ...
    'Lpp','Lwl','T','B','C_B','submerged','submergenceDepth'};
for f = 1:numel(scalarFields)
    name = scalarFields{f};
    cols{end+1} = sprintf('main_%s', name); %#ok<AGROW>
    vals{end+1} = sprintf('%.17g', m.(name)); %#ok<AGROW>
end

vectorFields = {'CG','CB'};
for f = 1:numel(vectorFields)
    name = vectorFields{f};
    v = m.(name);
    for k = 1:3
        cols{end+1} = sprintf('main_%s_%d', name, k); %#ok<AGROW>
        vals{end+1} = sprintf('%.17g', v(k));          %#ok<AGROW>
    end
end

% Raw B(omega) at every frequency and the zero-speed column, diagonal only
% (off-diagonal radiation damping is small at zero speed for this
% axisymmetric-ish body and not needed to check "never clipped" below);
% kept exactly as Capytaine wrote it (Fossen X-5 reply 2: no manual clip).
for k = 1:nFreq
    for i = 1:6
        cols{end+1} = sprintf('B_diag_f%02d_dof%d', k, i); %#ok<AGROW>
        vals{end+1} = sprintf('%.17g', B_mid(i,i,k));       %#ok<AGROW>
    end
end
for k = 1:nFreq
    cols{end+1} = sprintf('freq_%02d', k); %#ok<AGROW>
    vals{end+1} = sprintf('%.17g', freqs(k)); %#ok<AGROW>
end

outFile = fullfile(outDir, 'vessel_database', 'vessel_database_lauv_marie_mss_reference.csv');
fid = fopen(outFile, 'w');
fprintf(fid, '%s\n', strjoin(cols, ','));
fprintf(fid, '%s\n', strjoin(vals, ','));
fclose(fid);
fprintf('Wrote %s\n', outFile);
