"""Cross-flow drag by strip theory (numpy constants, CasADi algebra).

Equations (keys in References):

* Strip theory, 20 strips of width ``dx = L / 20`` at ``x_i``::

      Y = -1/2 rho T Cd_2D sum |v_r + x_i r| (v_r + x_i r) dx
      Z = -1/2 rho T Cd_2D sum |w_r + x_i q| (w_r + x_i q) dx
      M = -1/2 rho T Cd_2D sum x_i |w_r + x_i q| (w_r + x_i q) dx
      N = -1/2 rho T Cd_2D sum x_i |v_r + x_i r| (v_r + x_i r) dx

  (Fossen 2011, eqs. 6.91-6.92, p. 127, sway and yaw as integrals; MSS
  ``crossFlowDrag.m`` 36-38, 54-69, which adds heave and pitch).
* Hoerner: ``Cd_2D`` interpolated in ``B / (2T)`` from the table of MSS
  ``Hoerner.m`` 25-45 (Fossen 2011, Fig. 6.5, p. 128), the last value beyond
  it (``Hoerner.m`` 47-51).
* Cylinder: ``Re = U_cf B / nu``, ``U_cf = sqrt(v_r^2 + w_r^2)``,
  ``nu = 1e-6 m^2/s`` (DNV-RP-C205 2017, §6.6.1.1, p. 113; MSS
  ``cylinderDrag.m`` 78-80); ``C_D(Re)`` interpolated in the smooth-cylinder
  table of ``cylinderDrag.m`` 25-53 (digitised by MSS from the DNV-RP-C205
  curves, 2017 edition Figure 6-6, p. 118), end values outside it
  (``cylinderDrag.m`` 83-89); ``Cd_2D = C_D kappa(L/B)`` with the finite-
  length factor ``kappa`` of DNV-RP-C205 2017, Table 6-2, p. 120 (sub-
  critical row A below ``Re = 2e5``, super-critical row B from there;
  ``cylinderDrag.m`` 59-76, 92-110).

The tables below are copied from those MSS files (the full-precision Re
column, not the numpy source's 4-digit one) and are unchanged at MSS
``72656d1``, whose line numbers the comments cite; the kappa rows equal
DNV-RP-C205 2017 Table 6-2 entry by entry.

Options and deviations:

* ``strip_grid="midpoint"`` (default): 20 strips evaluated at their midpoints
  ``x_i = -L/2 + (i - 1/2) dx`` (``crossFlowDrag.m`` 56, MSS since
  2026-08-26). ``"endpoint"``: the 21 points ``-L/2 : dx : L/2`` of the
  earlier MSS and of the numpy source.
* The cylinder Reynolds number is on the diameter, ``Re = U_cf B / nu``
  (``cylinderDrag.m`` 80, MSS since ``ac77394``, DNV-RP-C205 2017 §6.6.1.1).
  The numpy source and MSS before that revision put it on the length; that
  form is not kept (owner's decision of 2026-10-06).
* Inputs outside their domain (a non-positive length, beam, draft or density,
  a beam-to-draft ratio below Hoerner's data) raise ``ValueError`` naming the
  input, where MSS returns inf or NaN (``Hoerner.m`` 48 extrapolates with
  ``interp1`` to NaN below the table).

Conventions: ``nu_r = [u v w p q r]`` relative to the water, body axes z down;
the cylinder coefficient is one number per call, set by the cross-flow speed;
strip height = ``draft``; ``x`` along the body from the CO, strips centred on
it.

Ported from the numpy source [MGM] ``dynamics/plant/matrices/hydrodynamics.py``:
``HydroForces._hoerner`` (95-100), ``_cylinder_drag`` (103-120),
``cross_flow_drag`` (127-178) and ``_get_strips`` (181-190); the same tables
sit in ``drag_models.py::DragModels.hoerner`` (16-64) and ``cylinder_drag``
(68-162).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 6,
    §6.4.3, eqs. 6.91-6.92 and Fig. 6.5, pp. 127-128.
[DNV-RP-C205 2017] DNV GL (2017). *DNVGL-RP-C205: Environmental conditions and
    environmental loads*, Recommended practice, edition August 2017. §6.6.1.1,
    p. 113; §6.7.1, Figure 6-6, pp. 116-118; §6.8, Table 6-2, p. 120.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/crossFlowDrag.m`` 36-69; ``HYDRO/cylinderDrag.m``
    25-110 (M. Seidl); ``HYDRO/Hoerner.m`` 25-51.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/``
    files and lines listed above.
"""

from dataclasses import dataclass

import casadi as ca
import numpy as np

# HYDRO/cylinderDrag.m 25-53: [Re, C_D] of a smooth circular cylinder
# (digitised from DNV-RP-C205 2017, Figure 6-6, p. 118).
CYLINDER_DRAG_DATA = np.array(
    [
        [10211.0405297256, 1.20769],
        [15543.5423211490, 1.20369],
        [24434.8681540938, 1.20823],
        [35719.8807222993, 1.21060],
        [60880.7508493947, 1.21093],
        [86174.7436190610, 1.20901],
        [118086.431856986, 1.21134],
        [149262.255516420, 1.21148],
        [214708.876786124, 1.21171],
        [230940.445891927, 1.20109],
        [271581.385758833, 1.16918],
        [297098.185911873, 1.11163],
        [325252.018159382, 1.00926],
        [361942.848517458, 0.89411],
        [409798.156149477, 0.70855],
        [475302.392082697, 0.53155],
        [524734.960427919, 0.40785],
        [578921.928291835, 0.32683],
        [638292.879987135, 0.28422],
        [744046.366191920, 0.29711],
        [853237.628592033, 0.32280],
        [1068685.70450627, 0.38909],
        [1393381.37138896, 0.47033],
        [1831850.38455424, 0.53878],
        [2258395.50114555, 0.58372],
        [2899011.91287642, 0.62868],
        [3663279.03082383, 0.64803],
        [4937305.49552761, 0.67808],
    ]
)

# HYDRO/cylinderDrag.m 59-76: [L/B, kappa], finite-length reduction
# (DNV-RP-C205 2017, Table 6-2 rows A and B, p. 120).
KAPPA_SUBCRITICAL_DATA = np.array(
    [[2, 0.58], [5, 0.62], [10, 0.68], [20, 0.74], [40, 0.82], [50, 0.87], [100, 0.98]],
    dtype=float,
)
KAPPA_SUPERCRITICAL_DATA = np.array(
    [[2, 0.80], [5, 0.80], [10, 0.82], [20, 0.90], [40, 0.98], [50, 0.99], [100, 1.00]],
    dtype=float,
)

# HYDRO/cylinderDrag.m 92: sub-critical below, super-critical from here on.
CRITICAL_REYNOLDS_NUMBER = 2e5

# HYDRO/cylinderDrag.m 79: 1 / nu_water with nu_water = 1e-6 m^2/s.
INVERSE_KINEMATIC_VISCOSITY = 1e6

# HYDRO/Hoerner.m 25-45: [B/(2T), C_D] (Fossen 2011, Fig. 6.5, p. 128);
# the last value holds beyond (Hoerner.m 49-50).
HOERNER_DRAG_DATA = np.array(
    [
        [0.0108623, 1.96608],
        [0.176606, 1.96573],
        [0.353025, 1.89756],
        [0.451863, 1.78718],
        [0.472838, 1.58374],
        [0.492877, 1.27862],
        [0.493252, 1.21082],
        [0.558473, 1.08356],
        [0.646401, 0.998631],
        [0.833589, 0.87959],
        [0.988002, 0.828415],
        [1.30807, 0.759941],
        [1.63918, 0.691442],
        [1.85998, 0.657076],
        [2.31288, 0.630693],
        [2.59998, 0.596186],
        [3.00877, 0.586846],
        [3.45075, 0.585909],
        [3.7379, 0.559877],
        [4.00309, 0.559315],
    ]
)

# crossFlowDrag.m 37.
NUMBER_OF_STRIPS = 20

DRAG_MODELS = ("cylinder", "hoerner")
STRIP_GRIDS = ("midpoint", "endpoint")


@dataclass(frozen=True)
class CrossFlowDragConstants:
    """Strip grid and drag-coefficient data computed before the graph."""

    water_density: float
    draft: float
    strip_width: float
    strip_positions: np.ndarray
    drag_model: str
    drag_coefficient: float
    diameter: float
    kappa_subcritical: float
    kappa_supercritical: float


def _clamped_interp(x: float, table: np.ndarray) -> float:
    """``interp1`` inside the table, end values outside (``cylinderDrag.m`` 83-108)."""
    return float(np.interp(x, table[:, 0], table[:, 1]))


def _clamped_interp_casadi(x: ca.SX, table: np.ndarray) -> ca.SX:
    """Piecewise-linear ``C_D(x)``, end values outside the table.

    Per segment the ``np.interp`` form ``f_j + s_j (x - x_j)``, selected with
    ``if_else`` so the result stays a pure SX expression.
    """
    xs, fs = table[:, 0], table[:, 1]
    value = ca.SX(fs[-1])
    for j in range(len(xs) - 2, -1, -1):
        slope = (fs[j + 1] - fs[j]) / (xs[j + 1] - xs[j])
        value = ca.if_else(x < xs[j + 1], slope * (x - xs[j]) + fs[j], value)
    return ca.if_else(x < xs[0], fs[0], value)


def _hoerner_drag(beam: float, draft: float) -> float:
    """``Hoerner.m`` 47-51 (``_hoerner``): 2-D coefficient of B/(2T).
    Deviation from ``Hoerner.m`` 48: below the table this raises, where MSS
    returns NaN."""
    ratio = beam / (2.0 * draft)
    if ratio < HOERNER_DRAG_DATA[0, 0]:
        raise ValueError(
            f"beam / (2 draft) = {ratio} lies below Hoerner's data "
            f"({HOERNER_DRAG_DATA[0, 0]}); MSS returns NaN there"
        )
    return _clamped_interp(ratio, HOERNER_DRAG_DATA)


def _strip_positions(length: float, strip_grid: str) -> tuple:
    dx = length / NUMBER_OF_STRIPS  # (crossFlowDrag.m 38)
    if strip_grid == "midpoint":
        index = np.arange(1, NUMBER_OF_STRIPS + 1, dtype=float)
        return -length / 2 + (index - 0.5) * dx, dx           # crossFlowDrag.m 56
    return dx * (np.arange(NUMBER_OF_STRIPS + 1) - NUMBER_OF_STRIPS / 2), dx  # _get_strips


def preprocess_cross_flow_drag(
    length: float,
    beam: float,
    draft: float,
    water_density: float,
    drag_model: str,
    strip_grid: str = "midpoint",
) -> CrossFlowDragConstants:
    """Strip grid and coefficient data of ``crossFlowDrag(L, B, T, nu_r, model)``.

    Source names: ``L``, ``B``, ``T``, ``rho``, ``drag_model``. For
    ``"cylinder"`` the beam is the diameter (MSS calls ``crossFlowDrag(L, D, D,
    ...)``), the beam/diameter the Reynolds number is built on; the length sets
    the strip positions and the aspect ratio, and ``kappa(L/B)`` is fixed here
    for both flow regimes; for ``"hoerner"`` the whole coefficient is fixed
    here.
    """
    values = {"length": length, "beam": beam, "draft": draft, "water_density": water_density}
    for name, value in values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")
    model = str(drag_model).lower()
    if model not in DRAG_MODELS:
        raise ValueError(f"drag_model must be one of {DRAG_MODELS}, got {drag_model!r}")
    if strip_grid not in STRIP_GRIDS:
        raise ValueError(f"strip_grid must be one of {STRIP_GRIDS}")

    positions, dx = _strip_positions(float(length), strip_grid)
    aspect_ratio = length / beam
    return CrossFlowDragConstants(
        water_density=float(water_density),
        draft=float(draft),
        strip_width=dx,
        strip_positions=positions,
        drag_model=model,
        drag_coefficient=_hoerner_drag(beam, draft) if model == "hoerner" else float("nan"),
        diameter=float(beam),
        kappa_subcritical=_clamped_interp(aspect_ratio, KAPPA_SUBCRITICAL_DATA),  # (cylinderDrag.m 93-99)
        kappa_supercritical=_clamped_interp(aspect_ratio, KAPPA_SUPERCRITICAL_DATA),  # (cylinderDrag.m 101-107)
    )


def _cylinder_drag_casadi(constants: CrossFlowDragConstants, nu_r: ca.SX) -> ca.SX:
    """``C_D(Re) kappa`` (``cylinderDrag.m`` 78-110, ``_cylinder_drag``)."""
    cross_flow_speed = ca.sqrt(nu_r[1] ** 2 + nu_r[2] ** 2)  # (cylinderDrag.m 78)
    # Re = v D / nu (DNV-RP-C205 2017, §6.6.1.1, p. 113; cylinderDrag.m 79-80)
    reynolds = cross_flow_speed * constants.diameter * INVERSE_KINEMATIC_VISCOSITY
    kappa = ca.if_else(  # (DNV-RP-C205 2017, Table 6-2, p. 120; cylinderDrag.m 92)
        reynolds < CRITICAL_REYNOLDS_NUMBER,
        constants.kappa_subcritical,
        constants.kappa_supercritical,
    )
    return _clamped_interp_casadi(reynolds, CYLINDER_DRAG_DATA) * kappa  # (cylinderDrag.m 83-89, 110)


def cross_flow_drag_casadi(constants: CrossFlowDragConstants) -> ca.Function:
    """``nu_r -> tau = [0 Yh Zh 0 Mh Nh]`` (Fossen 2011, eqs. 6.91-6.92, p. 127;
    MSS ``crossFlowDrag.m`` 54-69)."""
    nu_r = ca.SX.sym("nu_r", 6)
    if constants.drag_model == "cylinder":
        drag_coefficient = _cylinder_drag_casadi(constants, nu_r)
    else:
        drag_coefficient = ca.SX(constants.drag_coefficient)
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    # -1/2 rho T Cd_2D dx (Fossen 2011, eqs. 6.91-6.92, p. 127; crossFlowDrag.m 63-66)
    factor = -0.5 * constants.water_density * constants.draft * drag_coefficient * constants.strip_width
    y_sum = z_sum = m_sum = n_sum = ca.SX(0.0)
    for x in constants.strip_positions:
        horizontal = v_r + x * r
        vertical = w_r + x * q
        u_h = ca.fabs(horizontal) * horizontal  # (crossFlowDrag.m 61)
        u_v = ca.fabs(vertical) * vertical  # (crossFlowDrag.m 62)
        y_sum += u_h  # (crossFlowDrag.m 63)
        z_sum += u_v  # (crossFlowDrag.m 64)
        m_sum += x * u_v  # (crossFlowDrag.m 65)
        n_sum += x * u_h  # (crossFlowDrag.m 66)
    tau = ca.vertcat(0.0, factor * y_sum, factor * z_sum, 0.0, factor * m_sum, factor * n_sum)  # (crossFlowDrag.m 69)
    return ca.Function("cross_flow_drag", [nu_r], [tau], ["nu_r"], ["tau"])
