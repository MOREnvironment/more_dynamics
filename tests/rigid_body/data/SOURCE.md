# Frozen reference — rigid body and added mass

**Frozen** files: a changed byte in a CSV needs a decision record; a new MSS revision that changes an output adds a file `*_mss_<revision>.csv` beside the old one, never writes over it. Every generator below lives in this repository, finds MSS only through `MSS_DIR` (and MATLAB through `MATLAB_BIN`), and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/<block>/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here. Nothing here is relative to one machine.

## Files (rewritten 2026-10-07; MSS checked at `cc07579b63630342802a198c34c74ac7e096d5fb`)

| File | What it holds | Made by | Made at MSS | At MSS `cc07579` | sha256 |
|---|---|---|---|---|---|
| `inputs.csv` | 50 frozen cases: `x1..x12` = `[nu (6); eta (6)]`, `tau1..tau6` (shared with `../hydrostatics`, `../hydrodynamics`) | test input | — | — | `b653df012854130ebf7af4e916a78fda74abb0ffb0a2b3f9fe82fd593eb1cd48` |
| `matlab_reference_mss_current.csv` | 50 × 180: `M_RB`, `M_A`, `C_RB`, `C_A`, `C_total` of MSS `otter.m` (inertia about the combined CG) on `inputs.csv`, payload 25 kg at `[0.05 0 -0.35]` m | `../hydrodynamics/generate_catamaran_mss.m` (its columns, cut as text) | `99bf0b3` | identical (also at `ac77394`, `72656d1`) | `ce33d1401f586aaac6e3c09bae1b4b25157ec357311645ef8804140372229d06` |
| `rigid_body_rbody_mss_current.csv` | 600 × 123, 3 parameter sets: `rbody(m,R44,R55,R66,nu2,r_bG)` → `M_RB`, `C_RB`; `m2c(MRB,nu)` → `C_RB_m2c` | `generate_rigid_body_forms_mss.m` | `ac77394` | identical | `07c331d5d9e3a94a94952171a0315ad980e122d5f0b6ba0b259b52b01ef3bcc2` |
| `rigid_body_spheroid_mss_current.csv` | 400 × 194, 2 sets: `spheroid.m` → `M_RB`, `C_RB`; `imlay61.m` → `M_A`, `C_A`; `C_A_stab` = `C_A` with the eight `remus100.m` zeros | `generate_rigid_body_forms_mss.m` | `ac77394` | identical | `a3ab98e6559c98bf914dec76f5ea101d78e89b161b9fcb6be85397f84eac3abd` |
| `rigid_body_hull_mss_current.csv` | 400 × 207, 2 sets: `m = rho Cb L B T`, `A11` (`addedMassSurge.m`), `rbody.m` → `M_RB`, `C_RB`; `M_A = -diag(c .* [A11 m m Ig11 Ig22 Ig33])` (Ig about the CG); `m2c(MA,nu)` → `C_A`, `C_A_stab` | `generate_rigid_body_forms_mss.m` | `ac77394` | identical | `16e284ebcc97e9242651f45c0d6e68c51850e0e960777029d8bbca00801e4da6` |
| `rigid_body_m2c_3dof_mss_current.csv` | 400 × 27, 2 sets: `M3 = -[[Xu 0 0];[0 Yv Yr];[0 Yr Nr]]`, `m2c(M3,[u v r])` → `C3` | `generate_rigid_body_forms_mss.m` | `ac77394` | identical: these `M3` have `M(1,2) = M(1,3) = 0`, where `m2c.m`'s 3-DOF branch before and after its correction of 2026-10-07 (`9aef3ba`, `a3406cf`) agree | `d7bd59575b3150a3562ef69ff3d909307f28460f9c9390916780a7b05a46eacf` |
| `rigid_body_m2c_3dof_coupled_mss_cc07579.csv` | 402 × 23, 2 sets of fully coupled symmetric `M3` (set 1 `[[10 2 1];[2 20 3];[1 3 30]]`, set 2 `[[25 -1.5 4];[-1.5 60 7.5];[4 7.5 220]]`, test inputs): `m2c(M3,[u v r])` → `C3`; case 0 = `[1 2 3]` (set 1: `C3·nu3 = [-153, 51, 17]`), cases 1–200 = the states above | `generate_rigid_body_forms_mss.m` (section 5, added 2026-10-07) | `cc07579` | made here: the corrected 3-DOF branch (`m2c.m` 52–57, `p = M*nu`); before `9aef3ba` the branch dropped the surge coupling and gives a different `C3` on these matrices | `39cfa2db807ed127e29e8f50663b0a7c380f3e66cb6561904ede0b6c3442be75` |
| `spheroid_matlab_reference_mss_current.csv` | 50 × 150: `nu_r`, `M_RB`, `M_A`, `C_RB`, `C_A` of MSS `remus100.m` on `inputs.csv` | `../hydrodynamics/generate_spheroid_mss.m` (its columns) | `99bf0b3` | identical (these columns do not depend on cross-flow) | `8a1f0415412cd3cdedd2a9eab486438068dcf00d2abd1bc05485c46fd6b3dbf6` |
| `matlab_reference.csv` | 50 × 180, the same columns from the **template** Otter model (inline copy of `otter.m` as of 2026-02: inertia about the CO shifted a second time by `H`, latitude gravity) | `../hydrodynamics/generate_catamaran_template.m` (its columns) | MSS of 2026-02 | template record: at `cc07579` the generator gives `C_RB_*`, `C_total_*` different by ≤ 9.95e-14 (rounding through changed MSS helper functions); kept as frozen, no test reads it | `126106d05e78fb04dfcd2ec7844271f566fa016455e3322060298adc319a8ccb` |
| `generate_rigid_body_forms_mss.m` | the generator of the five `rigid_body_*` files | — | — | re-run 2026-10-07 at `72656d1` and `cc07579`: byte-identical outputs; section 5 (the coupled file) added the same day and re-run at `cc07579`: the four earlier files byte-identical | `e4852d394595cc05729e9efd933a534a83c69104fd26bf2deec4cc037670a47a` |

## Layout and conventions

- A matrix `X` is `X_01..X_36` **row-major** (3-DOF `X_01..X_09`); vectors one-based without padding (`nu_1..nu_6`, `r_bG_1..3`, `Rs_1..3`, `c_1..6`, `nu3_1..3`); `set`, `case` integers. Numbers written with `%.17g` (exact doubles) by `generate_rigid_body_forms_mss.m`; `writetable`'s 15 digits in the files cut from the two full-debug generators.
- `C_RB` is evaluated at `nu`; `C_A` at `nu_r = nu - nu_c`, `nu_c = [Vc cos(beta_c - psi), Vc sin(beta_c - psi), 0, 0, 0, 0]`, `Vc = 0.3` m/s, `beta_c = 30 deg`, `psi = x12`.
- Rigid-body forms: parameter sets are test inputs, echoed in every row; states: cases 1–50 = `x1..x6` of `inputs.csv`, cases 51–200 = `rng(20261006,'twister')`, uniform(−3, 3); the same `nu` feeds `C_RB` and `C_A`. The only transcribed MSS lines are `remus100.m` 207–210 (the zeroed `C_A` entries; 205–208 at `ac77394`) and the `otter.m` scaled-derivative pattern (153–158 at `cc07579`).
- The Otter inertia about the combined CG is MSS `otter.m` since its revision of 2026-04-20: `M_RB` rotational diagonal 15.52682 / 21.275 / 16.0125 kg·m²; the template file has 21.55470 / 28.26826 / 16.97788.

## MSS line numbers

Line numbers written in the generators are those of the MSS revision named beside them. From `72656d1` to `cc07579`: `otter.m` +1 from line 79 (a revision line), `m2c.m` +2 from line 32 and its 3-DOF branch 52–54 → 54–57 (rewritten), `remus100.m`, `rbody.m`, `spheroid.m`, `imlay61.m`, `addedMassSurge.m` unchanged.

## Readers

`tests/rigid_body/test_rigid_body_block.py`, `tests/rigid_body/test_rigid_body_forms.py`.
