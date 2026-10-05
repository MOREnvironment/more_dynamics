"""Cross-flow drag by strip theory (numpy constants, CasADi algebra).

Ported from ``more_generic_models`` ``dynamics/plant/matrices/hydrodynamics.py``:
``HydroForces._hoerner`` (lines 95-100), ``_cylinder_drag`` (103-120),
``cross_flow_drag`` (127-178) and ``_get_strips`` (181-190); the same tables sit
in ``drag_models.py::DragModels.hoerner`` (16-64) and ``cylinder_drag``
(68-162). Translated against MSS (MIT, T. I. Fossen)
``LIBRARY/modeling/crossFlowDrag.m`` lines 24-54, ``HYDRO/cylinderDrag.m``
lines 23-107 and ``HYDRO/Hoerner.m`` lines 25-51; the tables below are copied
from those two files at MSS ``99bf0b3`` (the full-precision Re column, not the
source's 4-digit one, owner decision E-10 a).

Options (owner decisions E-10 a, E-20 Q2 a and Q3 b), explicit arguments of
the preprocess:

* ``strip_grid="midpoint"`` (default): 20 strips evaluated at their midpoints
  (``crossFlowDrag.m`` since 2026-08-26). ``"endpoint"``: the 21 points
  ``-L/2 : dx : L/2`` of the earlier MSS and of the numpy source.
* ``cross_flow_reynolds_length="diameter"`` (default): ``Re = U_cf B / nu``.
  ``"length"``: ``Re = U_cf L / nu`` as ``cylinderDrag.m`` line 77.

Conventions: ``nu_r = [u v w p q r]`` relative to the water, body axes z down;
the cylinder coefficient is one number per call, set by the cross-flow speed
``sqrt(v_r^2 + w_r^2)``; kinematic viscosity 1e-6 m^2/s; strip height =
``draft``; ``x`` along the body from the CO, strips centred on it.
"""

from dataclasses import dataclass

import casadi as ca
import numpy as np

# HYDRO/cylinderDrag.m lines 23-51: [Re, C_D] of a circular cylinder.
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

# HYDRO/cylinderDrag.m lines 57-74: [L/B, kappa], finite-length reduction.
KAPPA_SUBCRITICAL_DATA = np.array(
    [[2, 0.58], [5, 0.62], [10, 0.68], [20, 0.74], [40, 0.82], [50, 0.87], [100, 0.98]],
    dtype=float,
)
KAPPA_SUPERCRITICAL_DATA = np.array(
    [[2, 0.80], [5, 0.80], [10, 0.82], [20, 0.90], [40, 0.98], [50, 0.99], [100, 1.00]],
    dtype=float,
)

# HYDRO/cylinderDrag.m line 89: sub-critical below, super-critical from here on.
CRITICAL_REYNOLDS_NUMBER = 2e5

# HYDRO/cylinderDrag.m line 77: 1 / nu with nu = 1e-6 m^2/s.
INVERSE_KINEMATIC_VISCOSITY = 1e6

# HYDRO/Hoerner.m lines 25-45: [B/(2T), C_D]; the last value holds beyond.
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

# crossFlowDrag.m line 25.
NUMBER_OF_STRIPS = 20

DRAG_MODELS = ("cylinder", "hoerner")
REYNOLDS_LENGTHS = ("diameter", "length")
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
    reynolds_length: float
    kappa_subcritical: float
    kappa_supercritical: float


def _clamped_interp(x: float, table: np.ndarray) -> float:
    """``interp1`` inside the table, end values outside (``cylinderDrag.m`` 90-104)."""
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
    """``Hoerner.m`` 47-51 (``_hoerner``): 2-D coefficient of B/(2T)."""
    ratio = beam / (2.0 * draft)
    if ratio < HOERNER_DRAG_DATA[0, 0]:
        raise ValueError(
            f"beam / (2 draft) = {ratio} lies below Hoerner's data "
            f"({HOERNER_DRAG_DATA[0, 0]}); MSS returns NaN there"
        )
    return _clamped_interp(ratio, HOERNER_DRAG_DATA)


def _strip_positions(length: float, strip_grid: str) -> tuple:
    dx = length / NUMBER_OF_STRIPS
    if strip_grid == "midpoint":
        index = np.arange(1, NUMBER_OF_STRIPS + 1, dtype=float)
        return -length / 2 + (index - 0.5) * dx, dx           # crossFlowDrag.m 41
    return dx * (np.arange(NUMBER_OF_STRIPS + 1) - NUMBER_OF_STRIPS / 2), dx  # _get_strips


def preprocess_cross_flow_drag(
    length: float,
    beam: float,
    draft: float,
    water_density: float,
    drag_model: str,
    cross_flow_reynolds_length: str = "diameter",
    strip_grid: str = "midpoint",
) -> CrossFlowDragConstants:
    """Strip grid and coefficient data of ``crossFlowDrag(L, B, T, nu_r, model)``.

    Source names: ``L``, ``B``, ``T``, ``rho``, ``drag_model``. For
    ``"cylinder"`` the beam is the diameter (MSS calls ``crossFlowDrag(L, D, D,
    ...)``) and ``kappa(L/B)`` is fixed here for both flow regimes; for
    ``"hoerner"`` the whole coefficient is fixed here.
    """
    values = {"length": length, "beam": beam, "draft": draft, "water_density": water_density}
    for name, value in values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")
    model = str(drag_model).lower()
    if model not in DRAG_MODELS:
        raise ValueError(f"drag_model must be one of {DRAG_MODELS}, got {drag_model!r}")
    if cross_flow_reynolds_length not in REYNOLDS_LENGTHS:
        raise ValueError(f"cross_flow_reynolds_length must be one of {REYNOLDS_LENGTHS}")
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
        reynolds_length=float(beam if cross_flow_reynolds_length == "diameter" else length),
        kappa_subcritical=_clamped_interp(aspect_ratio, KAPPA_SUBCRITICAL_DATA),
        kappa_supercritical=_clamped_interp(aspect_ratio, KAPPA_SUPERCRITICAL_DATA),
    )


def _cylinder_drag_casadi(constants: CrossFlowDragConstants, nu_r: ca.SX) -> ca.SX:
    """``C_D(Re) kappa`` (``cylinderDrag.m`` 76-107, ``_cylinder_drag``)."""
    cross_flow_speed = ca.sqrt(nu_r[1] ** 2 + nu_r[2] ** 2)
    reynolds = cross_flow_speed * constants.reynolds_length * INVERSE_KINEMATIC_VISCOSITY
    kappa = ca.if_else(
        reynolds < CRITICAL_REYNOLDS_NUMBER,
        constants.kappa_subcritical,
        constants.kappa_supercritical,
    )
    return _clamped_interp_casadi(reynolds, CYLINDER_DRAG_DATA) * kappa


def cross_flow_drag_casadi(constants: CrossFlowDragConstants) -> ca.Function:
    """``nu_r -> tau = [0 Yh Zh 0 Mh Nh]`` (``crossFlowDrag.m`` 39-54)."""
    nu_r = ca.SX.sym("nu_r", 6)
    if constants.drag_model == "cylinder":
        drag_coefficient = _cylinder_drag_casadi(constants, nu_r)
    else:
        drag_coefficient = ca.SX(constants.drag_coefficient)
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    factor = -0.5 * constants.water_density * constants.draft * drag_coefficient * constants.strip_width
    y_sum = z_sum = m_sum = n_sum = ca.SX(0.0)
    for x in constants.strip_positions:
        horizontal = v_r + x * r
        vertical = w_r + x * q
        u_h = ca.fabs(horizontal) * horizontal
        u_v = ca.fabs(vertical) * vertical
        y_sum += u_h
        z_sum += u_v
        m_sum += x * u_v
        n_sum += x * u_h
    tau = ca.vertcat(0.0, factor * y_sum, factor * z_sum, 0.0, factor * m_sum, factor * n_sum)
    return ca.Function("cross_flow_drag", [nu_r], [tau], ["nu_r"], ["tau"])
