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


## Frozen references of the MATLAB model file (`grethe_m_*`, 2026-10-23)

A second implementation of the same craft in MSS style (`CRAFT/USV/models/grethe.m`, `SIMgrethe.m`, `exGrethe.m`, on a local branch of an MSS fork, revision `99fdf4d` plus local commits; not in this repository). **A self-made reference, not a validation:** it was written from the same cited equations and the same part values, so agreement shows two implementations of one set of equations, never that the boat behaves so (field level is G7a). Generator: `generate_grethe_mss.m`, run with `MSS_DIR=<MSS checkout holding grethe.m> OUT_DIR=<this folder>` (and `TMPDIR` unset); seeds 20261023 (derivative cases), 20261024 (payload cases). Closed-loop run: `closed_loop/` (plot and summary of `SIMgrethe(true, <folder>)`), a demonstration, not a gate.

| File | Holds |
|---|---|
| `grethe_m_derivative.csv` | 300 states and commands (thrust N, angle rad): outboard wrench `tau_out`, `xdot` (thrust vector, restoring at the body origin), `xdot_x_ref` (restoring at x = -2.24 m), leg-force variants |
| `grethe_m_payload_derivative.csv` | 12 cases, 25 kg payload at a seeded position (the payload trim: green since the surface restoring has `payload_trim_form`) |
| `grethe_m_mass_matrix.csv`, `grethe_m_outboard.json` | `M` and the outboard limits read from the model file |
| `grethe_m_trajectory_<run>[_x_ref|_leg].csv`, `grethe_m_trajectory_tolerances*.csv` | open-loop RK4 runs (h = 0.05 s, command held 0.5 s) and the per-state tolerance (largest effect of a 1e-9 derivative error or of reordering the RK4 sum) |

Runs (`_ramp` files: thrust step and coast-down with the smoothstep residual; `xdot_ramp` columns likewise): straight from rest (400 N), coast-down from 2.40 m/s, turning circle (20 deg, 450 N), steering square wave (5 deg, 20 s), roll free decay (5 deg), pitch free decay (3 deg). No heading zigzag was run: the first version of the model file (Y_r = N_v = 0) was yaw-unstable in open loop and no switching rule recovered it; with the chosen yaw damping a zigzag could now be added.

**Disagreements found (not tuned):** (1) the part's `reference_point` [-2.24, 0, 0] places the restoring matrix at the outboard, not at the body origin (the point P of `Gmtrx.m`); the model file defaults to the origin and both are gated (`gate`, `x_ref`); (2) payload: the model file includes the payload lever-arm trim moment, the library has no form for it (red test); (3) the outboard leg normal force exists only in the model file (red test); (4) the model file ramps the residual resistance below Fn_V 0.6 (smoothstep, default), the library steps; the gates use the file's `ramp = false` set and a red test names the missing ramp form. **Values chosen in the model file, not in the part:** waterplane coefficient 0.41 (fitted to the logged 1.85 s roll period), T_yaw 4 s and N_v 2500 N m s/m (chosen for a yaw-stable open loop and about 7 deg/s at 20 deg steering); the gate tests pass them to the library as overrides, the part is unchanged.

**Resolved 2026-10-26 (Y3):** the four disagreements are closed and the overrides are gone. The `Monohull`'s surface restoring has `payload_trim_form` (the payload trim, gated), the `Propulsor` has the leg force housing and the thrust-command forms, the surge resistance has the ramped residual (`ittc_residual_ramped`), the part's `reference_point` is the origin, and the part carries the model file's own values (waterplane coefficient 0.41, T_yaw 4 s, N_v 2500); the gates in `test_monohull_grethe.py` run on the part with form switches only.

Author: Enio Krizman · Date: 2026-10-09 (references of 2026-10-23 added)

## The plant of the part (`grethe_m_plant_*`, 2026-10-26, MATLAB R2026a Update 4, Y3)

The named part "Grethe" follows the model file (E-145 Q7) and carries its outboard with the actuator equations the file leaves
to its caller. The model file takes the actual thrust and steering angle as inputs; `SIMgrethe.m` 151-154 and 257-259 give
the thrust a first-order lag and the steering a lag with rate and angle limits. The generator builds that plant, in
continuous form (`SIMgrethe.m` takes an Euler step; the part integrates the states with the vehicle's own RK4):
`T' = (T_c - T) / T_thrust`, `delta' = satlim((satlim(delta_c, +-delta_max) - delta) / T_delta, +-ddelta_max)`, the model file
called with its defaults (leg force, residual ramp, restoring at the origin, no payload). The generator reads
`T_thrust`, `T_delta`, `ddelta_max`, `delta_max` from the file's `lim` (`grethe_m_plant.json`), it types none of them.

| File | What it holds |
|---|---|
| `grethe_m_plant.json` | the four actuator values read from `lim` |
| `grethe_m_plant_derivative.csv` | G1: 300 seeded cases, `xa_01..14` = [actual thrust, actual angle, x(12) in MSS order], `uc_01..02` the commands (some beyond the limits: thrust beyond +-600 N, angle beyond +-0.9 rad), `xadot_01..14` the plant's derivative |
| `grethe_m_plant_trajectory_<run>.csv` | G3: the six runs of `exGrethe.m` (its lines 62-71) on the plant; `t, ui_01..02` (commands), `x_01..12`, `act_01..02` (actual thrust, angle), the commands held 10 steps, RK4 at 0.05 s, actuators at rest at the start |
| `grethe_m_plant_trajectory_tolerances.csv` | per run and state (the 14), the largest effect of a 1e-9 derivative error or of reordering the RK4 sum; at most 9e-5 |

**The references do not move with the MSS re-pin.** The generator was run on MSS `ed2ce7b` (a clean `git archive` of it) with
`grethe.m` and `exGrethe.m` of the fork's `grethe-vehicle` branch at `679a6e8` added: the 24 files that existed before come out
byte-identical, the 9 new ones were reproduced byte for byte by a second run. (`grethe.m` itself is in neither MSS revision;
`MSS_DIR` must hold both.) The earlier model-file sets (`gate`, `leg`, `x_ref`, `ramp`) stay as the option-by-option references:
the part reproduces them with ideal actuators and the matching form switches (`housing`, the surge resistance form, the
restoring point). `grethe_parameters.json` and the U8a CSVs above are the parameter set and references of the independent
block gates, kept as they were; the named part's values are the model file's (see the part's own `SOURCE.md`).
