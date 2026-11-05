"""Reduction of frequency-dependent added mass and potential damping to constant matrices (the power-based method), as
one CasADi block.

The tables ``A(omega)`` and ``B(omega)`` are 6 x 6 matrices in BODY FRD about one reference point of the hull (the point
the tables are about, fixed in the body), degrees of freedom ``[surge sway heave roll pitch yaw]``; units kg, kg m,
kg m^2 for ``A`` and N s/m, N s, N m s for ``B``. Each table is given flat, row-major per frequency (``n x 36`` values,
the layout a hydrodynamic database tool writes).

Equations (keys in References)::

    A_eq = integral A(omega) S_N(omega) d omega,        B_eq = integral B(omega) S_N(omega) d omega,
    S_N = S / m_0,   m_0 = integral S d omega                                       (Fossen 2025, the power-based method)

    S(omega) ~ omega^-5 exp(-5/4 (omega_p / omega)^4) gamma^exp(-(omega - omega_p)^2 / (2 sigma^2 omega_p^2)),
    sigma = 0.07 for omega <= omega_p and 0.09 otherwise                  (DNVGL-RP-C205 2017, 3.5.5.1-3, pp. 64-65)

The scale ``(5/16) H_s^2 omega_p^4 A_gamma`` of the spectrum cancels in ``S_N`` and is left out. The entries are
interpolated onto ``FINE_GRID_POINTS`` equally spaced frequencies with the shape-preserving piecewise cubic Hermite
interpolant (Fritsch and Carlson 1980: interior slopes the weighted harmonic mean of the two secants when they have the
same sign, else zero; end slopes the one-sided three-point formula with the shape-preserving corrections of Moler 2004,
sec. 3.4) and integrated with the trapezoid rule. The table entry at ``omega >= SENTINEL_FREQUENCY`` stands for the
infinite frequency of the data format and is not used; a first frequency of zero starts the fine grid at
``MINIMUM_FREQUENCY``.

The frequencies are structural: they set the intervals of the interpolant, so they are an argument of the builder (as
a form selector is), and the block's inputs are the two tables, ``omega_p`` and ``spectrum_gamma`` by name. The integral
is written as ``sum_k (a_k y_k + b_k d_k)`` over the knots, with the weights ``a_k``, ``b_k`` of the interpolant and the
trapezoid rule on the spectrum (functions of ``omega_p`` and ``gamma``) and the slopes ``d_k`` of the table values ``y_k``:
the same sum as interpolating first and integrating after, in another order of the additions.

The minimum, sign and harmonic-mean logic of the slopes is written with ``if_else``; every division has a guarded
denominator, so the block and its derivatives are finite for any table.

References
----------
[Fossen 2025] Fossen, T. I. (2025). Maneuvering coefficient estimation from frequency-dependent added mass and damping:
    a power-based approach. Ocean Engineering 341:122494. (Not read: the paper is not on disk; the method is stated as
    the equations above.)
[DNVGL-RP-C205] DNV GL (2017). Environmental conditions and environmental loads, DNVGL-RP-C205. Sec. 3.5.5.1-3.5.5.3,
    pp. 64-65 (Pierson-Moskowitz and JONSWAP spectra, sigma_a = 0.07, sigma_b = 0.09).
[Fritsch 1980] Fritsch, F. N., Carlson, R. E. (1980). Monotone piecewise cubic interpolation. SIAM J. Numer. Anal.
    17(2):238-246.
[Moler 2004] Moler, C. (2004). Numerical Computing with MATLAB. SIAM. Sec. 3.4 (piecewise cubic Hermite interpolation).

Author:    Enio Krizman
Date:      2026-11-05
"""
import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.shared.wiring import function_from

SENTINEL_FREQUENCY = 10.0      # rad/s: the table entry at or above it stands for infinite frequency (data format)
FINE_GRID_POINTS = 1000        # equally spaced frequencies of the integration (numerical resolution)
MINIMUM_FREQUENCY = 1e-6       # rad/s: replaces a first frequency of zero (the spectrum is singular there)
SIGMA_BELOW_PEAK = 0.07        # JONSWAP width for omega <= omega_p (DNVGL-RP-C205 2017, 3.5.5.3, p. 65)
SIGMA_ABOVE_PEAK = 0.09        # JONSWAP width for omega > omega_p (DNVGL-RP-C205 2017, 3.5.5.3, p. 65)


def power_based_reduction_parameters(frequency_count):
    """The block's declared parameters for ``frequency_count`` table entries (the sentinel included)."""
    n = int(frequency_count)
    return (
        Parameter("added_mass", (36 * n, 1), "kg, kg m, kg m^2",
                  "added-mass table A(omega), 6 x 6 row-major per frequency, about the table reference point"),
        Parameter("potential_damping", (36 * n, 1), "N s/m, N s, N m s",
                  "potential-damping table B(omega), 6 x 6 row-major per frequency, about the table reference point"),
        Parameter("omega_p", (1, 1), "rad/s", "peak frequency of the spectrum that weights the tables", 0.0,
                  minimum_exclusive=True),
        Parameter("spectrum_gamma", (1, 1), "1", "peak factor gamma of the JONSWAP spectrum", 1.0),
    )


def check_frequencies(frequencies):
    """The used frequencies (below the sentinel) as floats; refuses a grid that is not strictly increasing or has fewer
    than two used entries."""
    w = [float(v) for v in frequencies]
    if any(b <= a for a, b in zip(w[:-1], w[1:])):
        raise ValueError("frequencies must be strictly increasing")
    used = [v for v in w if v < SENTINEL_FREQUENCY]
    if len(used) < 2:
        raise ValueError(f"fewer than two frequencies below the sentinel {SENTINEL_FREQUENCY} rad/s")
    return w, used


def _fine_grid(start, stop):
    """``FINE_GRID_POINTS`` equally spaced frequencies from ``start`` to ``stop`` (both included)."""
    step = (stop - start) / (FINE_GRID_POINTS - 1)
    grid = [i * step + start for i in range(FINE_GRID_POINTS)]
    grid[-1] = stop
    return grid


def jonswap_shape(omega, omega_p, gamma):
    """Unnormalised JONSWAP spectrum at the frequency ``omega`` (rad/s): ``omega^-5 exp(-5/4 (omega_p/omega)^4)
    gamma^exp(-1/2 ((omega - omega_p)/(sigma omega_p))^2)`` (DNVGL-RP-C205 2017, 3.5.5.1-3.5.5.2, p. 64), the scale left
    out (it cancels in the normalised spectrum)."""
    sigma = ca.if_else(omega <= omega_p, SIGMA_BELOW_PEAK, SIGMA_ABOVE_PEAK)  # sigma_a, sigma_b (3.5.5.2, p. 65)
    peak = ca.exp(-0.5 * ((omega - omega_p) / (sigma * omega_p)) ** 2)  # (DNVGL-RP-C205 2017, 3.5.5.2, p. 64)
    return omega ** -5 * ca.exp(-1.25 * (omega_p / omega) ** 4) * gamma ** peak  # S_PM shape times gamma^peak (p. 64)


def _sign(value):
    return ca.sign(value)


def _slopes(x, y):
    """Slopes of the shape-preserving cubic Hermite interpolant at the knots ``x`` (floats) of the values ``y`` (list of
    SX of one size), [Fritsch 1980] interior and [Moler 2004] sec. 3.4 ends."""
    n = len(x)
    h = [x[k + 1] - x[k] for k in range(n - 1)]
    delta = [(y[k + 1] - y[k]) / h[k] for k in range(n - 1)]
    if n == 2:
        return [delta[0], delta[0]]
    d = [None] * n
    for k in range(1, n - 1):
        same = _sign(delta[k - 1]) * _sign(delta[k]) > 0  # both secants of one sign (Fritsch 1980)
        w1, w2 = 2.0 * h[k] + h[k - 1], h[k] + 2.0 * h[k - 1]
        left = ca.if_else(same, delta[k - 1], 1.0)  # guarded denominators: used only when ``same``
        right = ca.if_else(same, delta[k], 1.0)
        d[k] = ca.if_else(same, (w1 + w2) / (w1 / left + w2 / right), 0.0)  # weighted harmonic mean (Fritsch 1980)

    def end(h0, h1, d0, d1):
        s = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)  # one-sided three-point slope (Moler 2004, sec. 3.4)
        s = ca.if_else(_sign(s) != _sign(d0), 0.0, s)  # not against the first secant (Moler 2004, sec. 3.4)
        limit = ca.logic_and(ca.logic_and(_sign(s) != 0, _sign(d0) != _sign(d1)), ca.fabs(s) > ca.fabs(3.0 * d0))
        return ca.if_else(limit, 3.0 * d0, s)  # at most three times the secant at a sign change (Moler 2004)

    d[0] = end(h[0], h[1], delta[0], delta[1])
    d[-1] = end(h[-1], h[-2], delta[-1], delta[-2])
    return d


def _hermite_weights(x, fine):
    """For each fine point: ``(k, h00, h10 h, h01, h11 h)``, the interval and the Hermite basis of the interpolant."""
    n = len(x)
    out = []
    for value in fine:
        k = 0
        while k < n - 2 and x[k + 1] <= value:  # the interval [x_k, x_k+1) holding the point, the last one closed
            k += 1
        h = x[k + 1] - x[k]
        t = (value - x[k]) / h
        out.append((k, 2 * t ** 3 - 3 * t ** 2 + 1, (t ** 3 - 2 * t ** 2 + t) * h, -2 * t ** 3 + 3 * t ** 2,
                    (t ** 3 - t ** 2) * h))
    return out


def power_based_reduction(frequencies):
    """The block ``(added_mass, potential_damping, omega_p, spectrum_gamma) -> (added_mass_at_table_reference,
    potential_damping_at_table_reference)`` for the table frequencies ``frequencies`` (rad/s, increasing; an entry at or
    above ``SENTINEL_FREQUENCY`` is not used): the 6 x 6 equivalent matrices about the table reference point (module
    docstring)."""
    w, used = check_frequencies(frequencies)
    n, m = len(w), len(used)
    tables = {name: ca.SX.sym(name, 36 * n) for name in ("added_mass", "potential_damping")}
    omega_p, gamma = ca.SX.sym("omega_p"), ca.SX.sym("spectrum_gamma")
    fine = _fine_grid(MINIMUM_FREQUENCY if used[0] == 0.0 else used[0], used[-1])
    spectrum = [jonswap_shape(value, omega_p, gamma) for value in fine]
    spacing = [b - a for a, b in zip(fine[:-1], fine[1:])]
    trapezoid = [0.0] * len(fine)  # trapezoid weights of the fine grid
    for i, dx in enumerate(spacing):
        trapezoid[i] += 0.5 * dx
        trapezoid[i + 1] += 0.5 * dx
    m0 = sum(0.5 * dx * (spectrum[i] + spectrum[i + 1]) for i, dx in enumerate(spacing))  # m_0 (Fossen 2025)
    # knot weights of the value (a) and of the slope (b): integral of the interpolant times S_N, by the trapezoid rule
    a, b = [ca.SX(0.0)] * m, [ca.SX(0.0)] * m
    for (k, h00, h10, h01, h11), weight, s in zip(_hermite_weights(used, fine), trapezoid, spectrum):
        ws = weight * s
        a[k], b[k], a[k + 1], b[k + 1] = a[k] + h00 * ws, b[k] + h10 * ws, a[k + 1] + h01 * ws, b[k + 1] + h11 * ws
    a = [value / m0 for value in a]  # S_N = S / m_0 (Fossen 2025)
    b = [value / m0 for value in b]
    outputs = {}
    for name, table in tables.items():
        rows = [table[36 * k:36 * (k + 1)] for k in range(m)]  # the used frequencies: the first m (increasing grid)
        slopes = _slopes(used, rows)
        flat = sum((a[k] * rows[k] + b[k] * slopes[k] for k in range(m)), ca.SX.zeros(36))
        outputs[f"{name}_at_table_reference"] = ca.reshape(flat, 6, 6).T  # row-major 6 x 6
    return function_from("power_based_reduction", {**tables, "omega_p": omega_p, "spectrum_gamma": gamma}, outputs)
