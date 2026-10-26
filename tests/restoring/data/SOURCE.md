# Frozen reference — hydrostatics (submerged restoring forces, surface restoring matrix)

**Frozen** files: a changed byte in a CSV needs a decision record; a new MSS revision that changes an output adds a file `*_mss_<revision>.csv` beside the old one. Every generator lives in this repository, finds MSS only through `MSS_DIR` (MATLAB through `MATLAB_BIN`) and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/<block>/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here.

## Files (rewritten 2026-10-07; MSS checked at `cc07579b63630342802a198c34c74ac7e096d5fb`)

| File | What it holds | Made by | Made at MSS | At MSS `cc07579` | sha256 |
|---|---|---|---|---|---|
| `inputs.csv` | the 50 frozen cases (`x1..x12`, `tau1..tau6`), as in `../rigid_body` | test input | — | — | `b653df012854130ebf7af4e916a78fda74abb0ffb0a2b3f9fe82fd593eb1cd48` |
| `spheroid_matlab_reference_mss_current.csv` | 50 × 7: `M_RB_01` (the mass, kg) and `g_01..g_06` (the restoring vector, N and N·m) of MSS `remus100.m` / `gRvect.m` at the attitudes of `inputs.csv`; W = m g(63.446827 deg), neutrally buoyant | `../hydrodynamics/generate_spheroid_mss.m` (its columns, cut as text) | `99bf0b3` | identical | `22a64de9587f77ddecc4fb74cc296d70f7f6f3c376fdfb11365757c19de21a7c` |
| `spheroid_matlab_reference.csv` | the same 7 columns from the template REMUS model (inline copy of `remus100.m` as of 2026-02) — byte-identical to the file above | `../hydrodynamics/generate_spheroid_template.m` | MSS of 2026-02 | identical | `22a64de9587f77ddecc4fb74cc296d70f7f6f3c376fdfb11365757c19de21a7c` |
| `surface_gmtrx_mss_current.csv` | 50 × 44: `LIBRARY/modeling/Gmtrx.m` on random small-craft argument sets: `nabla, A_wp, GMT, GML, x_F, r_bP_1..3`, `G_01..G_36`; `r_bP` off the CO in every case; `Gmtrx.m` fixes `rho = 1025`, `g = 9.81` (lines 25–26) | `generate_surface_hydrostatics_mss.m` | `ac77394` | identical | `ed0209d192ed211f9f9404cd99e6d1b2c6ffcef004760904155d0d5d9ce9a53f` |
| `surface_chain_mss_current.csv` | 22 × 61: the MSS hydrostatic coefficient chain — `kind` 1 = `CRAFT/USV/models/otter.m` (twin hull, 20 payloads; row 1 = `mp = 25` kg at `[0.05 0 -0.35]` m), `kind` 2 = `mssExamples/exShipHydrostatics.m`, `kind` 3 = `CRAFT/SHIP/models/osv.m`; inputs `hull_count, L, B_hull, T, nabla, Cw, r_bg_3, y_hull, I_L_factor, x_F, rho, g, mp, rp_1..3`, outputs `A_wp, I_T, I_L, KB, BM_T, BM_L, GM_T, GM_L, G_01..36` (G at the CO). Here `otter.m`'s KB uses `L·B_pont` | `generate_surface_hydrostatics_mss.m` (the chain file is named by the `otter.m` KB form it finds) | `ac77394` | **superseded for kind 1** by the file below; kinds 2, 3 byte-equal; identical when run on MSS ≤ `72656d1` | `304c0266eacda2ce237345872f94c6612743731c88edc54bc5e4e79983fa5849` |
| `surface_chain_mss_cc07579.csv` | the same chain at `cc07579`: `otter.m:178` computes KB with the pontoon waterplane area `Aw_pont = Cw L B_pont` (MSS commit `cf349d4`, 2026-10-07). Differs from the file above in `KB`, `GM_T`, `GM_L`, `G_22` (G44), `G_29` (G55) of the 20 kind-1 rows only; row 1: KB 0.136585 → 0.127913 m, GM_T 1.232500 → 1.223828 m, G44 967.266 → 960.460 N·m, G55 2743.370 → 2736.565 N·m | `generate_surface_hydrostatics_mss.m` | `cc07579` | new | `b5a6cc92ece53e56edc2065685fd6ce76c2ce209589d5f127ffafee06ef49d03` |
| `generate_surface_hydrostatics_mss.m` | the generator of the three `surface_*` files; `rng(20261006,'twister')`; numbers `%.17g` | — | — | re-run 2026-10-07 at `72656d1` (outputs equal the `_current` files) and `cc07579` (the `_cc07579` chain) | `acf1555114213c5ef8680e6ebf589fc83ea5f6a90b5a8eabe5510995bbcd2264` |
| `gravity_mss_cc07579.csv` | 1 row: `mu_deg,g` — `INS/functions/gravity.m` run by MATLAB at the REMUS latitude 63.446827 deg; the frozen oracle `test_gravity_matches_frozen_mss_output` checks `more_transformations` gravity against (T1, E-62: an oracle is a frozen MSS/MATLAB output, not a Python re-implementation of the transform) | `generate_gravity_mss.m` | `cc07579` | new | `b42deb206c7e53b4b3d879639db8ea2348074116559de1574c05a7cb4d31e41e` |

## Layout and conventions

- CO on the waterline, z down (`KG = T - z_g`, as `exShipHydrostatics.m` `r_bB = [.. T-KB]`); `x_F`, `r_bP` from the CO in BODY. `G` row-major `G_01..G_36`.
- Kind-1 intermediate terms are read from `otter.m`'s workspace through a copy instrumented at run time (in `tempdir`, never saved), checked equal to the unmodified function (max |diff| = 0). `I_L_factor`: 0.8 per pontoon ×2 (`otter.m`), 0.7 (`exShipHydrostatics.m`, `osv.m`); `rho`, `g` of kinds 2–3 are `Gmtrx.m`'s.
- Morrish KB: the block's default is the waterplane form `KB = (1/3)(5T/2 − ∇/A_wp)` per hull — the form of `exShipHydrostatics.m`, `osv.m` and, since `cf349d4`, `otter.m`; a wall-sided hull has KB = T/2 exactly in this form. The block's `center_of_buoyancy_area="length_times_beam"` reproduces `otter.m` up to `72656d1` (`L·B_pont`).
- Checks in the tests: chain row 1's `G` equals the catamaran reference's `G` of the same MSS revision (`../hydrodynamics/catamaran_matlab_reference_mss_{current,cc07579}.csv`) to 0.

## MSS line numbers

Written at the revision named beside them. From `72656d1` to `cc07579`: `otter.m` +1 from line 79 (KB 177 → 178, text changed); `osv.m` unchanged before line 187; `exShipHydrostatics.m`, `Gmtrx.m`, `gRvect.m`, `gravity.m`, `spheroid.m` unchanged.

## Readers

`tests/restoring/test_submerged_hydrostatics_block.py`, `tests/restoring/test_surface_hydrostatics_block.py`.

## Re-pin to MSS `99fdf4d` (2026-10-23, MATLAB R2026a Update 4; `99fdf4dbdfe0b3f124e26a689c6aabe0beb5654f`, 41 commits after `cc07579`)

`spheroid_matlab_reference_mss_current.csv` rewritten in place (`imlay61.m` `rho` 1025, `23c307e`; the hydrostatics columns the generator cuts out of the spheroid run). `surface_chain_mss_cc07579.csv`, `surface_gmtrx_mss_current.csv` and `gravity_mss_cc07579.csv` are byte-identical at `99fdf4d` (`Gmtrx.m` and `otter.m` read the same 1025 / 9.81 from `mssConstants.m`). `generate_surface_hydrostatics_mss.m` reads `ws.r_bG` (`otter.m` renamed `rg`).

## Payload trim at MSS `ed2ce7b` (2026-10-26, MATLAB R2026a Update 4; Y3, payload trim form)

| File | What it holds | Made by | Made at MSS |
|---|---|---|---|
| `payload_trim_mss_ed2ce7b.csv` | 48 cases (payload masses 0, 10, 25, 50 kg at two positions, six seeded states each, propeller speeds inside the limits, a current): `case_id, mp, rp_1..3, x_1..12, n_1..2, Vc, beta_c`, then for MSS's own payload trim (`otter.m` 259-269) `eta0_1..6` (`eta_0`), `g0_1..6`, `g_mss_1..6` (the restoring vector `G*eta` after the shift), `G_01..36` (row-major), `xdot_mss_1..12` (the unmodified function's derivative); and for the corrected form (lever-arm moment only in `g_0`, kinematics at the actual attitude) `eta0c_1..6`, `g_corr_1..6`, `xdot_corr_1..12`, built from the same call's instrumented workspace | `generate_payload_trim_mss.m` (seed 20261026, numbers `%.17g`; the instrumented copy of `otter.m` equals the unmodified function, max difference 0, checked in the run) | `ed2ce7b` |

Two runs on this machine gave byte-identical files (sha256 `9f1a905873a4fae7a1aac48c41e35b71175b0f3652a4ede7786106419975ee68`). Reader: `test_payload_trim.py`. Over the cases with a payload the two forms' derivatives differ by up to 3.2 (the heave acceleration, m/s^2), so the gate tells them apart.
