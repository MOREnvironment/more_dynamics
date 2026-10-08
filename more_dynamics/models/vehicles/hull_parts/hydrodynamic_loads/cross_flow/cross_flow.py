"""Cross-flow drag by strip theory (CasADi).

``cross_flow_drag_parameters()`` declares the four numbers of the block
(``length``, ``beam``, ``draft``, ``water_density``: name, shape, SI unit,
meaning, admissible range); ``cross_flow_drag_casadi(drag_model=...,
strip_grid=...)`` builds ``(nu_r, <parameters by name>) -> tau``. The strip
positions, the Hoerner coefficient and the finite-length factor are built in
the graph from the parameters; the digitised tables are published data of
the method (module constants), not vehicle numbers.

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
* Interpolation inside a table is ``f_j + (x - x_j) (f_{j+1} - f_j) /
  (x_{j+1} - x_j)`` on the segment that holds ``x`` (MSS ``interp1``,
  linear), selected with ``if_else`` so the result stays one CasADi
  expression.

The tables below are copied from those MSS files (the full-precision Re
column) and are unchanged at MSS ``cc07579``, whose line numbers the comments
cite; the kappa rows equal DNV-RP-C205 2017 Table 6-2 entry by entry.

Options and deviations:

* ``strip_grid="midpoint"`` (default): 20 strips evaluated at their midpoints
  ``x_i = -L/2 + (i - 1/2) dx`` (``crossFlowDrag.m`` 56, MSS since
  2026-08-26). ``"endpoint"``: the 21 points ``-L/2 : dx : L/2`` of the
  earlier MSS.
* The cylinder Reynolds number is on the diameter, ``Re = U_cf B / nu``
  (``cylinderDrag.m`` 80, MSS since ``ac77394``, DNV-RP-C205 2017 §6.6.1.1).
  The length form of MSS before that revision is not kept.
* A beam-to-draft ratio below Hoerner's data is refused by
  ``check_cross_flow_drag_values`` naming the ratio, where MSS returns NaN
  (``Hoerner.m`` 48 extrapolates with ``interp1`` to NaN below the table);
  inside the graph the first value holds there. Non-positive dimensions are
  refused by the declared ranges.

Conventions: ``nu_r = [u v w p q r]`` relative to the water, body axes z down;
the cylinder coefficient is one number per call, set by the cross-flow speed;
strip height = ``draft``; ``x`` along the body from the CO, strips centred on
it.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 6,
    §6.4.3, eqs. 6.91-6.92 and Fig. 6.5, pp. 127-128.
[DNV-RP-C205 2017] DNV GL (2017). *DNVGL-RP-C205: Environmental conditions and
    environmental loads*, Recommended practice, edition August 2017. §6.6.1.1,
    p. 113; §6.7.1, Figure 6-6, pp. 116-118; §6.8, Table 6-2, p. 120.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``LIBRARY/modeling/crossFlowDrag.m``
    36-69; ``HYDRO/cylinderDrag.m`` 25-110 (M. Seidl); ``HYDRO/Hoerner.m``
    25-51.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, check_values, symbols

# HYDRO/cylinderDrag.m 25-53: [Re, C_D] of a smooth circular cylinder
# (digitised from DNV-RP-C205 2017, Figure 6-6, p. 118).
CYLINDER_DRAG_DATA = (
    (10211.0405297256, 1.20769),
    (15543.5423211490, 1.20369),
    (24434.8681540938, 1.20823),
    (35719.8807222993, 1.21060),
    (60880.7508493947, 1.21093),
    (86174.7436190610, 1.20901),
    (118086.431856986, 1.21134),
    (149262.255516420, 1.21148),
    (214708.876786124, 1.21171),
    (230940.445891927, 1.20109),
    (271581.385758833, 1.16918),
    (297098.185911873, 1.11163),
    (325252.018159382, 1.00926),
    (361942.848517458, 0.89411),
    (409798.156149477, 0.70855),
    (475302.392082697, 0.53155),
    (524734.960427919, 0.40785),
    (578921.928291835, 0.32683),
    (638292.879987135, 0.28422),
    (744046.366191920, 0.29711),
    (853237.628592033, 0.32280),
    (1068685.70450627, 0.38909),
    (1393381.37138896, 0.47033),
    (1831850.38455424, 0.53878),
    (2258395.50114555, 0.58372),
    (2899011.91287642, 0.62868),
    (3663279.03082383, 0.64803),
    (4937305.49552761, 0.67808),
)

# HYDRO/cylinderDrag.m 59-76: [L/B, kappa], finite-length reduction
# (DNV-RP-C205 2017, Table 6-2 rows A and B, p. 120).
KAPPA_SUBCRITICAL_DATA = ((2.0, 0.58), (5.0, 0.62), (10.0, 0.68), (20.0, 0.74), (40.0, 0.82), (50.0, 0.87), (100.0, 0.98),)
KAPPA_SUPERCRITICAL_DATA = ((2.0, 0.80), (5.0, 0.80), (10.0, 0.82), (20.0, 0.90), (40.0, 0.98), (50.0, 0.99), (100.0, 1.00),)

# HYDRO/cylinderDrag.m 92: sub-critical below, super-critical from here on.
CRITICAL_REYNOLDS_NUMBER = 2e5

# HYDRO/cylinderDrag.m 79: 1 / nu_water with nu_water = 1e-6 m^2/s.
INVERSE_KINEMATIC_VISCOSITY = 1e6

# HYDRO/Hoerner.m 25-45: [B/(2T), C_D] (Fossen 2011, Fig. 6.5, p. 128);
# the last value holds beyond (Hoerner.m 49-50).
HOERNER_DRAG_DATA = (
    (0.0108623, 1.96608),
    (0.176606, 1.96573),
    (0.353025, 1.89756),
    (0.451863, 1.78718),
    (0.472838, 1.58374),
    (0.492877, 1.27862),
    (0.493252, 1.21082),
    (0.558473, 1.08356),
    (0.646401, 0.998631),
    (0.833589, 0.87959),
    (0.988002, 0.828415),
    (1.30807, 0.759941),
    (1.63918, 0.691442),
    (1.85998, 0.657076),
    (2.31288, 0.630693),
    (2.59998, 0.596186),
    (3.00877, 0.586846),
    (3.45075, 0.585909),
    (3.7379, 0.559877),
    (4.00309, 0.559315),
)

# crossFlowDrag.m 37.
NUMBER_OF_STRIPS = 20

DRAG_MODELS = ("cylinder", "hoerner")
STRIP_GRIDS = ("midpoint", "endpoint")

_DECLARED = (
    Parameter("length", (1, 1), "m", "body length L (strip positions, aspect ratio L/B)", 0.0,
              minimum_exclusive=True),
    Parameter("beam", (1, 1), "m", "beam B, the diameter of a cylinder (Reynolds number)", 0.0,
              minimum_exclusive=True),
    Parameter("draft", (1, 1), "m", "draft T, the strip height", 0.0, minimum_exclusive=True),
    Parameter("water_density", (1, 1), "kg/m^3", "water density rho", 0.0, minimum_exclusive=True),
)


def cross_flow_drag_parameters():
    """The declared parameter set: a tuple of ``Parameter`` (name, shape, SI
    unit, meaning, admissible range), in the order of the block's inputs."""
    return _DECLARED


def _clamped_interp_casadi(x, table):
    """Piecewise-linear ``f(x)`` on a table, end values outside it
    (``cylinderDrag.m`` 83-108; ``Hoerner.m`` 47-51 above the table).

    Per segment ``f_j + s_j (x - x_j)``, the linear ``interp1`` of MSS,
    selected with ``if_else`` so the result stays a pure SX expression.
    """
    xs, fs = [row[0] for row in table], [row[1] for row in table]
    value = ca.SX(fs[-1])  # beyond the table: the last value (cylinderDrag.m 86-89; Hoerner.m 49-50)
    for j in range(len(xs) - 2, -1, -1):
        slope = (fs[j + 1] - fs[j]) / (xs[j + 1] - xs[j])  # linear interp1 on segment j (cylinderDrag.m 84)
        value = ca.if_else(x < xs[j + 1], slope * (x - xs[j]) + fs[j], value)  # (cylinderDrag.m 84)
    return ca.if_else(x < xs[0], fs[0], value)  # below the table: the first value (cylinderDrag.m 83)


def _strip_positions(length, strip_grid):
    """Strip positions ``x_i`` and width ``dx`` (symbols in ``length``)."""
    dx = length / NUMBER_OF_STRIPS  # (crossFlowDrag.m 38)
    if strip_grid == "midpoint":
        return [-length / 2 + (i - 0.5) * dx for i in range(1, NUMBER_OF_STRIPS + 1)], dx  # (crossFlowDrag.m 56)
    # the 21 points -L/2 : dx : L/2 of MSS before 2026-08-26 (crossFlowDrag.m 56 replaced them)
    return [dx * (i - NUMBER_OF_STRIPS / 2) for i in range(NUMBER_OF_STRIPS + 1)], dx


def _cylinder_drag_casadi(p, nu_r):
    """``C_D(Re) kappa(L/B)`` (``cylinderDrag.m`` 78-110)."""
    cross_flow_speed = ca.sqrt(nu_r[1] ** 2 + nu_r[2] ** 2)  # (cylinderDrag.m 78)
    # Re = v D / nu (DNV-RP-C205 2017, §6.6.1.1, p. 113; cylinderDrag.m 79-80)
    reynolds = cross_flow_speed * p["beam"] * INVERSE_KINEMATIC_VISCOSITY
    aspect_ratio = p["length"] / p["beam"]  # L/B (cylinderDrag.m 93, 101)
    kappa = ca.if_else(  # (DNV-RP-C205 2017, Table 6-2, p. 120; cylinderDrag.m 92)
        reynolds < CRITICAL_REYNOLDS_NUMBER,
        _clamped_interp_casadi(aspect_ratio, KAPPA_SUBCRITICAL_DATA),  # (cylinderDrag.m 93-99)
        _clamped_interp_casadi(aspect_ratio, KAPPA_SUPERCRITICAL_DATA),  # (cylinderDrag.m 101-107)
    )
    return _clamped_interp_casadi(reynolds, CYLINDER_DRAG_DATA) * kappa  # (cylinderDrag.m 83-89, 110)


def _check_selectors(drag_model, strip_grid):
    model = str(drag_model).lower()
    if model not in DRAG_MODELS:
        raise ValueError(f"drag_model must be one of {DRAG_MODELS}, got {drag_model!r}")
    if strip_grid not in STRIP_GRIDS:
        raise ValueError(f"strip_grid must be one of {STRIP_GRIDS}, got {strip_grid!r}")
    return model


def cross_flow_drag_casadi(*, drag_model, strip_grid="midpoint"):
    """The block: ``(nu_r, length, beam, draft, water_density) -> tau =
    [0 Yh Zh 0 Mh Nh]`` (Fossen 2011, eqs. 6.91-6.92, p. 127; MSS
    ``crossFlowDrag.m`` 54-69).

    Contract
    --------
    Keywords choose the graph: ``drag_model`` in ``DRAG_MODELS`` (required;
    for ``"cylinder"`` the beam is the diameter, as MSS calls
    ``crossFlowDrag(L, D, D, ...)``), ``strip_grid`` in ``STRIP_GRIDS``; an
    unknown value raises ``ValueError``. Inputs: ``nu_r`` (6x1, relative to
    the water, body axes, SI) and the parameters of
    ``cross_flow_drag_parameters()`` by name. Output: ``tau`` (6x1, the force
    on the vehicle). No value is checked inside the graph
    (``check_cross_flow_drag_values`` does that).
    """
    model = _check_selectors(drag_model, strip_grid)
    p = symbols(_DECLARED)
    nu_r = ca.SX.sym("nu_r", 6)
    if model == "cylinder":
        drag_coefficient = _cylinder_drag_casadi(p, nu_r)
    else:
        drag_coefficient = _clamped_interp_casadi(p["beam"] / (2.0 * p["draft"]), HOERNER_DRAG_DATA)  # Cd_2D(B/(2T)) (Hoerner.m 47-51)
    positions, dx = _strip_positions(p["length"], strip_grid)
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    # -1/2 rho T Cd_2D dx (Fossen 2011, eqs. 6.91-6.92, p. 127; crossFlowDrag.m 63-66)
    factor = -0.5 * p["water_density"] * p["draft"] * drag_coefficient * dx
    y_sum = z_sum = m_sum = n_sum = ca.SX(0.0)
    for x in positions:
        horizontal = v_r + x * r  # (crossFlowDrag.m 61)
        vertical = w_r + x * q  # (crossFlowDrag.m 62)
        u_h = ca.fabs(horizontal) * horizontal  # (crossFlowDrag.m 61)
        u_v = ca.fabs(vertical) * vertical  # (crossFlowDrag.m 62)
        y_sum += u_h  # (crossFlowDrag.m 63)
        z_sum += u_v  # (crossFlowDrag.m 64)
        m_sum += x * u_v  # (crossFlowDrag.m 65)
        n_sum += x * u_h  # (crossFlowDrag.m 66)
    tau = ca.vertcat(0.0, factor * y_sum, factor * z_sum, 0.0, factor * m_sum, factor * n_sum)  # (crossFlowDrag.m 69)
    names = [d.name for d in _DECLARED]
    return ca.Function("cross_flow_drag", [nu_r, *[p[name] for name in names]], [tau],
                       ["nu_r", *names], ["tau"])


def check_cross_flow_drag_values(values, *, drag_model):
    """Numbers checked: ``{name: ca.DM}``. ``check_values`` (names, shapes,
    finite, positive), and for ``"hoerner"`` the ratio ``B / (2T)`` must not
    lie below Hoerner's data (MSS returns NaN there, ``Hoerner.m`` 48): the
    refusal names the ratio."""
    model = _check_selectors(drag_model, "midpoint")
    numbers = check_values(_DECLARED, values)
    if model == "hoerner":
        ratio = float(numbers["beam"]) / (2.0 * float(numbers["draft"]))  # B/(2T) (Hoerner.m 47)
        if ratio < HOERNER_DRAG_DATA[0][0]:
            raise ValueError(
                f"beam / (2 draft) = {ratio} lies below Hoerner's data "
                f"({HOERNER_DRAG_DATA[0][0]}); MSS returns NaN there"
            )
    return numbers
