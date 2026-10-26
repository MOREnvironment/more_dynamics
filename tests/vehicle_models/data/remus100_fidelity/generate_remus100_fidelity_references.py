"""Generator of ``remus100_fidelity_turn_and_stability.csv``: for every case of
``fidelity_parts.cases()``, the straight-running equilibrium at the cruise
shaft speed and zero fin angles, the roots of the sway-yaw block of the
linearisation there (Fossen 2011, eqs. 12.17-12.21, p. 351; criterion C > 0,
Theorem 12.2, eq. 12.24) and the largest real part over all motion states,
and the steady yaw rate of a 4-degree rudder step (positive command, the
part's convention) held from the start of a 60 s run from 1.54 m/s
(RK4, h = 0.05 s, SIMremus100.m 60), averaged over the last 10 s.

Writes under ``OUT_DIR`` (default ``<tempdir>/more_references/vehicles/remus100_fidelity/``),
never into this folder; a re-run is checked against the frozen file with ``cmp``.

Run: ``python generate_remus100_fidelity_references.py`` from this folder.

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 12, eqs. 12.17-12.24, p. 351.
Fossen, T. I. MSS, MIT, @ ed2ce7b: CRAFT/AUV/SIMremus100.m 60;
    LIBRARY/numericalMethods/rk4.m 26-38.
Prestero, T. (2001). Verification of a six-degree of freedom simulation
    model for the REMUS autonomous underwater vehicle. MSc thesis, MIT/WHOI.
    Table 8.1, p. 60 (u = 1.54 m/s); p. 53 (a rudder of about 4 degrees).

Author:    Enio Krizman
Date:      2026-10-13
"""

import csv
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fidelity_parts as fp  # noqa: E402
from reference_torpedo_auv import TorpedoAuv, jacobian, rk4_step, straight_running, sway_yaw_roots  # noqa: E402

STEP = 0.05                       # SIMremus100.m 60
DURATION = 60.0
AVERAGE = 10.0
SPEED = 1.54                      # Prestero 2001, Table 8.1, p. 60
RUDDER = np.deg2rad(4.0)          # Prestero 2001, p. 53
FILE = "remus100_fidelity_turn_and_stability.csv"
COLUMNS = ("case", "pressure_centre_fraction", "cruise_u", "equilibrium_residual",
           "sway_yaw_root_1_re", "sway_yaw_root_1_im", "sway_yaw_root_2_re", "sway_yaw_root_2_im",
           "max_re_motion", "turn_rate_mean", "turn_rate_spread")


def start_state():
    x = np.zeros(12)
    x[6] = SPEED
    return x


def row(name, fraction, part):
    model = TorpedoAuv(part)
    zero = np.array(fp.rudder_command(part, 0.0), float)
    x_eq, residual = straight_running(model.xdot, zero, start_state())
    a = jacobian(model.xdot, x_eq, zero)
    roots = sorted(sway_yaw_roots(a), key=lambda z: (-z.real, z.imag))
    turn = np.array(fp.rudder_command(part, RUDDER), float)
    x, rates = start_state(), []
    n, keep = int(round(DURATION / STEP)), int(round(AVERAGE / STEP))
    for k in range(n):
        x = rk4_step(model.xdot, x, turn, STEP)
        if k >= n - keep:
            rates.append(x[11])
    rates = np.array(rates)
    return {"case": name, "pressure_centre_fraction": "" if fraction is None else f"{fraction:.2f}",
            "cruise_u": repr(float(x_eq[6])), "equilibrium_residual": repr(float(residual)),
            "sway_yaw_root_1_re": repr(float(roots[0].real)), "sway_yaw_root_1_im": repr(float(roots[0].imag)),
            "sway_yaw_root_2_re": repr(float(roots[1].real)), "sway_yaw_root_2_im": repr(float(roots[1].imag)),
            "max_re_motion": repr(float(np.linalg.eigvals(a).real.max())),
            "turn_rate_mean": repr(float(rates.mean())), "turn_rate_spread": repr(float(rates.max() - rates.min()))}


def main():
    out = Path(os.environ.get("OUT_DIR", Path(tempfile.gettempdir()) / "more_references" / "vehicles" / "remus100_fidelity"))
    out.mkdir(parents=True, exist_ok=True)
    with open(out / FILE, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        for name, fraction, part in fp.cases():
            writer.writerow(row(name, fraction, part))
    print(out / FILE)


if __name__ == "__main__":
    main()
