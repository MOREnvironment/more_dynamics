"""DUNE VSIM fin assembly: force ``K (delta / delta_max) U^2`` per fin (CasADi).

``vsim_fins_parameters(fin_count=...)`` declares the numbers (name, shape, SI
unit, meaning, admissible range); ``vsim_fins_casadi(fin_count=...)`` builds
``(delta, nu_r, <parameters by name>) -> tau``. ``fin_count`` (a selector)
sets the size of ``delta`` and of the two per-fin tables.

Equations (keys in References), per fin ``i`` with maximum force
``F_i = [F_x, F_y, F_z]`` at ``r_i``, deflection ``delta`` saturated at
``+-delta_max``:

* force ``f = U^2 (F_i / 2) (delta / delta_max)``, ``U^2 = |nu_r[0:3]|^2``
  (DUNE ``Fin.cpp`` 56-67, 86-93);
* torque ``t = U^2 [sqrt(y^2 + z^2) |F_x| / 2, x |F_y| / 2, x |F_z| / 2]
  (delta / delta_max)`` (``Fin.cpp`` 69-83), mapped to the wrench as
  ``[K, M, N] = [t_0, t_2, t_1]`` (``Fin.cpp`` 95-101).

``max_forces`` is a direct force gain of DUNE's simulator law: an identified
coefficient or a comparison form, not a physics value (a fin's physical
force follows from the local flow, its geometry and its lift and drag
curves, not from this gain).

This is DUNE's simulator law re-implemented in behaviour (not a rigid-body
moment ``r x f``). The speed ``U`` is ``|nu_r[0:3]|`` here (``Fin.cpp`` takes
it from its caller). Deviation from ``Fin.cpp``: below ``U = 1e-6`` m/s the
wrench is zero (a guard of this module).

References
----------
[DUNE] LSTS, Universidade do Porto. *DUNE: Unified Navigation Environment*,
    EUPL v1.1, ``src/Simulators/VSIM/VSIM/Fin.cpp`` 56-102 (B. Terra,
    J. Braga), revision ``555ef0b`` of the vehicle's DUNE source (read,
    re-implemented in behaviour, not copied).

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from more_dynamics.models.shared.force_producer_common import saturate

ZERO_SPEED = 1e-6  # m/s: below it the wrench is zero (guard of this module)


def _check_fin_count(fin_count):
    if not isinstance(fin_count, int) or isinstance(fin_count, bool) or fin_count < 1:
        raise ValueError(f"fin_count must be a positive integer, got {fin_count!r}")


def vsim_fins_parameters(*, fin_count):
    """The declared parameter set for ``fin_count`` fins: a tuple of
    ``Parameter`` in the order of the block's inputs."""
    _check_fin_count(fin_count)
    return (
        Parameter("max_forces", (fin_count, 3), "N",
                  "rows: maximum force [F_x, F_y, F_z] of each fin; direct force gain, identified "
                  "coefficient or comparison form, not a physics value"),
        Parameter("positions", (fin_count, 3), "m", "rows: CO -> fin, body axes (FRD)"),
        Parameter("max_deflection", (1, 1), "rad", "deflection limit delta_max", 0.0, minimum_exclusive=True),
    )


def vsim_fins_casadi(*, fin_count):
    """The block: ``(delta, nu_r, max_forces, positions, max_deflection) -> tau``.

    Contract
    --------
    ``fin_count`` a positive integer (else ``ValueError``). Inputs: ``delta``
    (``fin_count`` x 1, rad, saturated inside), ``nu_r`` (6x1, relative to
    the water) and the parameters by name. Output: ``tau`` (6x1, N and N m).
    """
    declared = vsim_fins_parameters(fin_count=fin_count)
    p = symbols(declared)
    delta = ca.SX.sym("delta", fin_count)
    nu_r = ca.SX.sym("nu_r", 6)
    limit = p["max_deflection"]
    normalized = saturate(delta, -limit, limit) / limit  # act / max_act (Fin.cpp 59-63, 66)
    forces, r = p["max_forces"], p["positions"]
    force_gain = (forces / 2.0).T  # F_i / 2, one column per fin (Fin.cpp 66)
    moment_gain = ca.vertcat(
        (ca.sqrt(r[:, 1] ** 2 + r[:, 2] ** 2) * ca.fabs(forces[:, 0]) / 2).T,  # (Fin.cpp 78-79)
        (r[:, 0] * ca.fabs(forces[:, 1]) / 2).T,  # (Fin.cpp 81)
        (r[:, 0] * ca.fabs(forces[:, 2]) / 2).T,  # (Fin.cpp 82)
    )
    speed_squared = ca.sumsqr(nu_r[0:3])  # U^2 = |nu_r[0:3]|^2
    force = speed_squared * (force_gain @ normalized)  # (Fin.cpp 92-93)
    torque = speed_squared * (moment_gain @ normalized)  # (Fin.cpp 95-101)
    tau = ca.vertcat(force, torque[0], torque[2], torque[1])  # [K, M, N] = [t_0, t_2, t_1] (Fin.cpp 97-101)
    tau = ca.if_else(ca.sqrt(speed_squared) <= ZERO_SPEED, ca.DM.zeros(6), tau)  # zero-speed guard of this module
    names = [d.name for d in declared]
    return ca.Function("vsim_fins", [delta, nu_r, *[p[name] for name in names]], [tau],
                       ["delta", "nu_r", *names], ["tau"])
