# Frozen references — Grethe's `Monohull` (U8a, no MSS model of Grethe)

**Frozen** files: a changed byte in a CSV or in `grethe_parameters.json` needs a decision record. The generator is pure Python + numpy (no CasADi, no `more_dynamics` import, no `MSS_DIR`/`MATLAB_BIN`): it is a second, independent implementation of the cited equations (G5 against the published anchors — Fossen 2011 and Radojčić et al. 2014 — and a G2 check against the library's own existing formulas where they already match, e.g. rigid body, Coriolis, metacentric restoring, Hoerner cross-flow, ITTC friction). There is no MSS model of Grethe (she is not a Fossen/MSS vehicle), so this is not a second-implementation-against-MATLAB check like the REMUS 100 and Otter references: it stands on the published equations alone.

Run: `python3 generate_grethe_references.py` in this folder; a re-run is checked with `cmp` against the committed files (deterministic arithmetic, no `MSS_DIR`/`MATLAB_BIN`; `skew` is taken from `more_transformations`, as every library generator does, AGENTS.md Transforms -- not run with `-I`, U8b).

Made 2026-10-09, library commit `5a527af` (branch `generic-models-port`, the commit the R9 porter handed off against).

## Parameter set

`grethe_parameters.json`: every value Grethe's `Monohull` composition needs, with its `kind` (`document` / `published` / `estimate` / `derived` / `gap` / `selector` / `owner decision` / `unidentified`) and `place` (the dossier row or inbox entry it comes from). Sources: `agents-more/60_working/2026-10-07_A43_vehicle_parameter_dossiers.md` Part 2 (M1–M50), `agents-more/20_sources/monohull_structure_for_grethe.md` Section 1.1 (drawing measurements) and Section 3.3 (option names and the G0 default set), `agents-more/80_owner/inbox/E-101_a64_grethe_monohull_choices_remus_density_at_repin.md` (Q1, Q2, Q4, Q5 — the owner's clicked decisions). `Y_r = N_v = 0`, marked `unidentified` (E-101 Q2 a). The CO is at the midpoint of the waterline length, on the waterline, `reference_point = [-2.24, 0, 0]` (E-101 Q4 a). `waterline_length` (4.63 m) and the overall `length` (5.2 m) are **both** declared because every RA14-keyed row (the residual coefficient, the RA14 wetted-surface option) is a **sensitivity case for the unmeasured projected chine length `L_P`** (A-64 Section 1.2 footnote 7, Codex pass A finding 6), not a validated substitution — the surge-resistance file below carries both, labelled `L463` and `L52`. No August data (safeguard rule 20): every number here is a document value, a published prior or an estimate already on record before this job: no log or sea-day value is used.

Several `kind: "owner decision"` / entries name a **library option that does not exist yet** (`added_mass_mass_basis`, `wetted_surface_method = "regression_table"`, `surge_resistance = "ittc_residual"`, `manoeuvring_damping = "linear_coupled"`): the `place` field says so explicitly. U8b (porter) builds them; this job only freezes what they should compute against.

## Files

| File | What it holds | Equations (cited in the generator) |
|---|---|---|
| `grethe_parameters.json` | the parameter set above | — |
| `grethe_rigid_body_added_mass.csv` | 4 rows (payload 0 / 300 kg x added-mass basis hull / displaced): `mass`, the combined CG `x_g,y_g,z_g`, the inertia diagonal `I11,I22,I33`, `M_A_11..M_A_66` (the vehicle's own `added_mass_matrix` diagonal, `M_A = -diag(derivatives)`, positive) | Fossen 2011 eq. 3.34 p. 50 (hull + point payload, MSS `otter.m` 94-98, 122-128); `addedMassSurge.m` 32-33, `otter.m` 152-159 (scaled derivatives, on the hull mass **or** the displaced mass, E-101 Q5 a) |
| `grethe_coriolis.csv` | 30 seeded 6-DOF states (seed 20261009, still water): `nu_01..06`, `C_RB_01..36`, `C_A_full_01..36`, `C_A_munk_removed_01..36` (column-major, as the payload's `matrix()` helper reads) | Fossen 2011 eq. 3.27 p. 50 ("co" form), eq. 6.43 p. 120 (`m2c`), eq. 6.52 p. 121 (Munk moment; MSS `remus100.m` 207-210 removes it) |
| `grethe_restoring.csv` | 2 rows (payload 0 / 300 kg): `draft`, `displaced_volume`, `wetted_surface_mumford`, `transverse_metacentric_height`, `longitudinal_metacentric_height`, `G_33`, `G_44`, `G_55` (at the waterline-midpoint CO, `reference_point = [-2.24, 0, 0]`, `hull_count = 1`) | Fossen 2011 eq. 4.24 p. 65, eqs. 4.32-4.35 pp. 65-66, eq. 7.250 p. 181; MSS `Gmtrx.m` 25-37, `exShipHydrostatics.m` 32-56, `otter.m` 121-128, 172-193; `XuuITTC.m` 38 (Mumford wetted surface) |
| `grethe_cross_flow.csv` | 30 seeded `nu_r` (seed 20261010), the bare-hull draft of `grethe_restoring.csv`: `nu_r_01..06`, `tau_overall_01..06` (today's basis — the overall length/beam 5.2/2.15 m, what `single_hull()` actually feeds the strip now), `tau_waterline_01..06` (the future basis — waterline length/beam 4.63/1.97 m, once `waterline_length`/`waterline_beam` are wired into the strip, A-64 Section 3.3 row R11) | Fossen 2011 eqs. 6.91-6.92 p. 127; MSS `crossFlowDrag.m` 54-69, `Hoerner.m` 25-51 (Hoerner's rectangular-section table, copied verbatim from the library's own `models/cross_flow/cross_flow.py` — published data, not re-digitised) |
| `grethe_surge_resistance.csv` | 33 rows, `u` 0.0–3.2 m/s step 0.1, both length cases (`L463` = waterline 4.63/1.97 m, `L52` = overall 5.2/2.15 m): `wetted_surface_mumford_<case>`, `C_R_RA14_<case>`, `S_RA14_<case>` (RA14's own wetted-surface option), `Fn_V_<case>`, `X_ittc_only_<case>` (today's `"ittc"` form, `C_R = 0`), `X_ittc_residual_<case>` (the not-yet-built `"ittc_residual"` form: `C_T = (1+k) C_F + C_R`, Fossen 2011 eq. 6.82's own decomposition, `C_R` from RA14) | Fossen 2011 eqs. 6.82-6.85 p. 125; RA14 Appendix 1 p. 24 (R_COEF, S_COEF, copied with citation from `agents-more/60_working/A64_scripts/a64_hull_regime_resistance.py`, itself transcribed and verified there — rule 4: not retyped from the paper in this job). **`C_R` is clamped to 0 below `Fn_V = 0.6`** (A-64 hidden assumption 3: the regression does not cover the hull there) |
| `grethe_damping.csv` | 20 seeded `(v_r, r_r)` pairs (seed 20261011), the bare-hull mass matrix: `Y_v`, `Y_r` (= 0), `N_v` (= 0), `N_r`, `tau_Y`, `tau_N` | `Y_v = -M22/T_sway`, `N_r = -M66/T_yaw` (the library's own `time_constant_damping_surface` formula, `otter.m` 203, 207, reproduced independently here); `Y_r = N_v = 0`, unidentified (E-101 Q2 a) |

## What this job does **not** freeze (rule 21, deferred to later jobs)

- **The outboard / leg force** (R13–R16): no force-producer composition exists yet (U12, the propulsor skeleton). No wrench, no turning-circle or coast-down trajectory is frozen here — those need thrust and steering.
- **Modulus (second-order) sway-yaw damping terms**: left at the placeholder `0` (A-64 Section 3.3), not identified.
- **The 2026 static-trim check** (−4.11°, −0.55°): an **owner design argument**, not a measurement (A-64 Section 1.1, E-101 addendum); the test contract carries it as a labelled, skipped-with-reason test until X-9 item G14.
- **RA14's own length `L_P`** (the projected chine length) is unmeasured; `L463`/`L52` above are sensitivity cases, never treated as that measurement.

## Readers

`tests/vehicle_models/test_monohull_grethe.py`.

Author: Enio Krizman · Date: 2026-10-09
