"""Pointwise (algebraic) step of the decomposition-coordination method for a general Young function.

Given w = eta/r + grad u, find nu minimising  G(x,|nu|) + (r/2)|nu - w|^2, i.e.

    nu = t w/|w|,    g(x,t) + r t = r |w|,   t >= 0                      (t = 0 if g(x,0+) >= r|w|).

phi(t) = g(t) + r t is strictly increasing (g is non-decreasing, r > 0), phi(0+) = g(0+) >= 0 and
phi(|w|) >= r|w|, so the root lies in [0, |w|].  It is found by bisection in log t (about 50 steps give a relative
accuracy ~1e-12 over 690 orders of magnitude, so even a tiny root, e.g. p -> 1, is resolved). Only monotonicity of g is used:
no smoothness, no Delta_2.
"""
import numpy as np

_LOGSPAN = 690.0          # smallest root resolved:  |w| * exp(-690) ~ 1e-300 |w|


def solve_scalar(gfun, s, r=1.0, iters=52):
    """Solve  gfun(t) + r t = s,  t >= 0  elementwise (s = r|w| >= 0).  ``gfun`` maps an array of t to g values
    (it may close over per-element parameters of the same shape as ``s``).  Returns t (0 where the root is < 1e-300 |w|,
    in particular where g(0+) >= s)."""
    s = np.asarray(s, float)
    t = np.zeros_like(s)
    pos = s > 0
    if not pos.any():
        return t
    idx = np.flatnonzero(pos.ravel())
    ss = s.ravel()[idx]
    hi_log = np.log(ss / r * (1 + 1e-12))
    lo_log = hi_log - _LOGSPAN
    gsub = lambda tt: gfun(tt, idx)
    # elements whose phi(e^lo) >= s have root ~ 0
    tiny = gsub(np.exp(lo_log)) + r * np.exp(lo_log) >= ss
    for _ in range(iters):
        mid = 0.5 * (lo_log + hi_log)
        tm = np.exp(mid)
        f = gsub(tm) + r * tm - ss
        neg = f < 0
        lo_log = np.where(neg, mid, lo_log)
        hi_log = np.where(neg, hi_log, mid)
    root = np.exp(0.5 * (lo_log + hi_log))
    root[tiny] = 0.0
    t.ravel()[idx] = root
    return t


def update_nu(gfun, wx, wy, r=1.0):
    """Pointwise nu = argmin G(|nu|) + (r/2)|nu - w|^2,  w = (wx, wy)."""
    mag = np.hypot(wx, wy)
    t = solve_scalar(gfun, r * mag, r)
    scale = np.divide(t, mag, out=np.zeros_like(mag), where=mag > 0)
    return wx * scale, wy * scale
