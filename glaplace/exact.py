"""Explicit / semi-explicit solutions for testing.

radial_torsion(G, c, R)
    On the disk B_R: -div( g(|grad u|) grad u/|grad u| ) = c, u = 0 on |x| = R.   The solution is radial and
    g(|u'(r)|) = c r/2  (2D), so  |u'(r)| = g^{-1}(c r/2)  (= 0 where c r/2 <= g(0+): dead core), and
    u(r) = int_r^R g^{-1}(c s/2) ds.   g may depend on r (radially symmetric weights a(r)); g^{-1} is computed by bisection
    and the integral by composite Gauss-Legendre (accuracy ~1e-10 unless g^{-1} is very singular).

manufactured(G_expr, u_expr)
    Given G(t, x, y) and u(x, y) as sympy expressions (strings), returns f = -div(g(x,y,|grad u|) grad u/|grad u|) and the
    exact u, grad u as numpy-callables.  u must have a non-vanishing gradient where g(t)/t is singular.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import CubicSpline


def _ginv(Gf, s, r):
    """Solve g(r, tau) = s, tau >= 0 (elementwise bisection in log; 0 if s <= g(r,0+))."""
    s = np.asarray(s, float)
    out = np.zeros_like(s)
    pos = s > 0
    if not pos.any():
        return out
    ss, rr = s[pos], r[pos]
    zero = np.zeros_like(rr)
    # bracket [lo, hi]: grow hi until g(hi) >= s
    hi = np.ones_like(ss)
    for _ in range(200):
        bad = Gf.g(hi, rr, zero) < ss
        if not bad.any():
            break
        hi = np.where(bad, hi * 2.0, hi)
    lo_log, hi_log = np.log(hi) - 690.0, np.log(hi)
    tiny = Gf.g(np.exp(lo_log), rr, zero) >= ss
    for _ in range(60):
        mid = 0.5 * (lo_log + hi_log)
        neg = Gf.g(np.exp(mid), rr, zero) < ss
        lo_log = np.where(neg, mid, lo_log)
        hi_log = np.where(neg, hi_log, mid)
    v = np.exp(0.5 * (lo_log + hi_log))
    v[tiny] = 0.0
    out[pos] = v
    return out


def radial_torsion(Gf, c=1.0, R=1.0, n=6000):
    """Returns (u, grad_u, u_r) as callables of (x, y); u_r(r) = |u'(r)|."""
    gl_x, gl_w = np.polynomial.legendre.leggauss(6)
    edges = np.linspace(0.0, R, n + 1)
    a, b = edges[:-1, None], edges[1:, None]
    xs = (0.5 * (b - a) * gl_x + 0.5 * (a + b)).ravel()
    ws = (0.5 * (b - a) * gl_w).ravel()
    tau = _ginv(Gf, c * xs / 2.0, xs).reshape(n, 6)
    cell = (tau * ws.reshape(n, 6)).sum(axis=1)             # int over each cell of |u'|
    cum = np.concatenate([[0.0], np.cumsum(cell)])          # int_0^{r_j} |u'|
    total = cum[-1]
    U = total - cum                                         # u(r_j) = int_{r_j}^R |u'|
    spl = CubicSpline(edges, U)

    def u_r(r):
        r = np.asarray(r, float)
        return _ginv(Gf, c * r / 2.0, r)

    def u(x, y):
        r = np.minimum(np.hypot(x, y), R)
        return spl(r)

    def grad_u(x, y):
        r = np.hypot(x, y)
        tr = u_r(r)
        sc = np.divide(-tr, r, out=np.zeros_like(r), where=r > 0)
        return sc * x, sc * y

    return u, grad_u, u_r


def manufactured(G_expr, u_expr, f_scale=1.0):
    """f = -div(g(x,y,|grad u|) grad u/|grad u|) for given sympy strings G(t,x,y), u(x,y)."""
    import sympy as sp
    x, y = sp.symbols("x y", real=True)
    t = sp.Symbol("t", positive=True)
    ns = {"x": x, "y": y, "t": t}
    G = sp.sympify(G_expr, locals=ns)
    u = sp.sympify(u_expr, locals=ns)
    g = sp.diff(G, t)
    ux, uy = sp.diff(u, x), sp.diff(u, y)
    mag = sp.sqrt(ux ** 2 + uy ** 2)
    coef = g.subs(t, mag) / mag
    f = -(sp.diff(coef * ux, x) + sp.diff(coef * uy, y))
    mod = "numpy"
    ff = sp.lambdify((x, y), f, mod)
    uf = sp.lambdify((x, y), u, mod)
    gxf, gyf = sp.lambdify((x, y), ux, mod), sp.lambdify((x, y), uy, mod)

    def bc(fun):
        return lambda X, Y: np.asarray(fun(X, Y), float) * np.ones_like(np.asarray(X, float))

    gx_, gy_ = bc(gxf), bc(gyf)
    return bc(ff), bc(uf), (lambda X, Y: (gx_(X, Y), gy_(X, Y)))
