# Frozen reference — propulsor parts and compositions (skeleton `Propulsor`)

**Frozen** files: a changed byte needs a decision record. The generator lives here and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/propulsor/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here. Nothing here needs MSS or MATLAB. The REMUS propeller rows, the two-thruster rows, the Wageningen rows and the outboard parameter set the propulsor gates also read are the existing frozen files of `tests/force_producers/data/` (not copied): `remus100_actuators_mss_current.csv`, `differential_thruster_matlab.csv`, `generate_differential_thruster_template.m`, `wageningen_kt_kq_mss_current.csv`, `parameter_sets.json` (`outboard_rpm`), and `tests/vehicle_models/data/remus100/remus100_parameters_consistent.json` (the propeller numbers, each with its `remus100.m` line).

## Files (made 2026-10-10)

| File | What it holds | Made by | sha256 |
|---|---|---|---|
| `jet_nozzle_reference.csv` | 600 rows (3 parameter sets × 200): the shared repository's jet-nozzle plugin's `generated_force` (6 outputs) and its dynamics (thrust and nozzle-angle derivatives) at structured states and inputs (zero, full ahead and astern, both angle limits, both signs of command) and seeded ones (seed 20261010, states and commands beyond the clamps), 17 significant digits | `generate_jet_nozzle_reference.py` | `1df0e56cc053837fe2b8bbc68f5494a25ebba1b454343e25616a5494bda505f7` |
| `jet_nozzle_test_parts.json` | the three parameter sets of the jet (this library's part names): `straight` (trim 0), `trimmed` (trim 0.2 rad), `biased` (trim −0.15 rad, yaw bias 0.05 rad, asymmetric angle limits, other thrust and times). **Arbitrary test numbers**, chosen to differ from the plugin's defaults and from each other; not any craft's | typed 2026-10-10 | `33b818476501514d2158c022958416e96416cf0adfd77dc4fd0e3680dc60a002` |
| `generate_jet_nozzle_reference.py` | builds his plugin through rpp's `ComponentContext` (`vehicle_models/vehicle_contract.py::context_of`, as the library's other plugin tests do), reads its payload with `RppCasadiGraph`, and writes the rows. Run from the library root; two runs were byte-identical | — | `f21f7603efd9e4f16c54e872ef378aa7c0b1e078d4630a8bf5408b8a259155ff` |
| `waterjet_published_values.json` | six numbers typed from the saved full texts, each with its table or page: Tajima et al. 2015 bollard coefficient (p. 378); Han et al. 2017 Table III `b11 k1`, `b11 k2` (p. 3354); Li & Kao 2023 nozzle outlet diameter (Table 1, p. 2280), rotor speed (Table 3, p. 2286), rotor torque (Table 6, p. 2288) | typed 2026-10-10 from `source-sim/0_literature/agents/` and `3_sim_to_real/` | `82290c276521af5a1c6572b66e59fb74da19f273102e12fad93e760639fc0ff1` |

## What was recorded from, and how

- The plugin is `more_dynamics/plugins/force_producers/jet_nozzle.py` of the shared repository, file sha256 `aa92025e4234e79927ddd04258a6cd7a671af6a6dde29a60ba62158cbc3df332`, library revision `8707610` (branch `generic-models-port`; the file's last change is `6b645e8`, 2026-10-02). His file is read, never edited. The `JetNozzle` parameter names are mapped `position → location`, `nozzle_velocity → nozzle_velocity_radpersec`.
- The rows are the only thing the propulsor gates read from his plugin. `test_propulsor_live_reference.py` runs the live plugin as an extra, only when it is present, and skips with a reason when it is not; if it ever disagrees with these rows his file changed, and the file's hash above says which.

## Typed numbers and their pages (rule: no number from memory)

| Value | Printed | Where |
|---|---|---|
| waterjet bollard coefficient | `F = 1.89×10⁻⁴ N²`, N in rpm, mean relative error 7.6 % over 10–90 % torque | Tajima, Hazuku, Shimizu, Ode, Takamasa (2015), J. JIME 50(3):376–381, eq. 7 and Fig. 3, p. 378 |
| thrust model of a waterjet USV | `T = k₁ ε^(2/3) + k₂ ε^(1/3) u`; `b11 k1 = 0.6006`, `b11 k2 = −0.0119` | Han, Xiong, He, Gu, Li (2017), IEEE TIE 64(4), eq. 7 p. 3350, Table III p. 3354 |
| jet-ski pump nozzle outlet diameter | `93.946` mm | Li & Kao (2023/2024), CMES 138(3), Table 1, p. 2280 |
| rotor speed at the trial top speed | `5950` rpm | Li & Kao, Table 3, p. 2286 |
| rotor torque | `391.550` N·m | Li & Kao, Table 6, p. 2288 |

Known from the memo (not used as a number): Li & Kao's "net thrust" in Table 6 is the rotor's, not the craft's (its 12.9 kN at 31.5 m/s would exceed the 244 kW shaft power); only the geometry, speed and torque are used.

## Layout and conventions

- BODY axes FRD; forces N, moments N·m about the CO; angles rad; relative velocity `nu_r` relative to the water; density fed by the vehicle.
- Command order [drive; steering] and state order [drive; steering], his: `[desired_thrust, desired_normalized_nozzle_angle]` and `[thrust, nozzle_angle]`.
- Rows: CSV with a header, one case per row, `set` names the parameter set.
