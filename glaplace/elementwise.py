"""Element-wise integrands for the v4 solver (NEW module; existing modules are untouched).

On a triangulation with one-point sampling (or a resolved / minorant / majorant quadrature) the integrand is a family
of Young functions G_T(t), one per triangle T.  This module represents such a family by per-element arrays and provides,
fully vectorised:

    G(t)             G_T(t)
    g(t)             g_T(t) = G_T'(t)      (right derivative; g(0+) may be > 0)
    prox(s, r, t0)   the root t >= 0 of  g_T(t) + r t = s   (pointwise step of the ADMM, s = r|w|)
    conj(s)          a RIGOROUS UPPER BOUND for the convex conjugate  G_T^*(s) = sup_{t>=0} (s t - G_T(t))  (exact closed
                     form where available, otherwise s*t_hi - G(t_lo) with a verified bracket t_lo <= t* <= t_hi)

Classes
    PowerSum(c, p)          G_T(t) = sum_k c_{k,T} t^{p_{k,T}} / p_{k,T}  (c >= 0, p > 1): p(x), p-q, two-phase, double phase
    ExpMinus1MinusT(Nt,w,b) w (e^t - 1 - t) + b t^2/2                     (conjugate (1+s)log(1+s) - s if b = 0, w = 1)
    ExpMinus1(Nt)           e^t - 1        (g(0+) = 1)                    (conjugate 0 for s <= 1, s log s - s + 1 for s > 1)
    Generic(Young, xb, yb)  any glaplace.young.Young object: bisection (existing code) for prox, bisection for conj

Pointwise step.  In the variable tau = log t the function phi(tau) = g(e^tau) + r e^tau - s is convex and increasing for
PowerSum and the exponentials, so Newton's method in tau started at an upper bound tau_hi of the root (log(s/r), or a
sharper class-specific bound) decreases monotonically to the root.  Warm starts (previous |nu|) are clipped to
tau <= tau_hi; if they lie left of the root the first Newton step lands right of it (convexity) and the iteration is
monotone afterwards.  Elements not converged after ``maxit`` steps (never observed) fall back to the bisection of
glaplace.nonlinear.  Relative accuracy ~1e-15 in t.
"""
from __future__ import annotations
import numpy as np

from .young import Young, PowerLaw, Sum, Exp, ExpLiteral
from .nonlinear import solve_scalar


# ------------------------------------------------------------------------------------------------------------------
def _newton_log(F, s, r, tau_hi, tau0, maxit=200, tol=1e-12):
    """Newton in tau = log t for g(t) + r t = s on the integrand F (all elements of F active, s > g(0+)).
    tau_hi: upper bound of the root (phi(tau_hi) >= 0); tau0 <= tau_hi start.  Returns (t, indices not converged)."""
    n = s.size
    tau = tau0.copy()
    live = np.arange(n)
    Fl, sl, hil = F, s, tau_hi
    for _ in range(maxit):
        tl = tau[live]
        tt = np.exp(tl)
        g, tdg = Fl.g_tdg(tt)
        phi = g + r * tt - sl
        dphi = tdg + r * tt                               # d phi / d tau
        with np.errstate(divide="ignore", invalid="ignore"):
            step = np.where(dphi > 0, phi / dphi, 0.0)
        new = np.minimum(tl - step, hil)
        done = (np.abs(new - tl) <= tol * np.maximum(1.0, np.abs(new))) | (phi == 0) | ~np.isfinite(step)
        tau[live] = np.where(np.isfinite(new), new, tl)
        keep = ~done
        if not keep.any():
            live = live[:0]
            break
        if keep.sum() < 0.5 * keep.size:                  # compress the active set (slicing parameters once)
            live, Fl, sl, hil = live[keep], Fl.subset(np.flatnonzero(keep)), sl[keep], hil[keep]
        else:                                             # keep shapes; converged entries simply stay put
            pass
    return np.exp(tau), live


class ElementIntegrand:
    """Base class: subclasses implement G(t), g(t), g_tdg(t) -> (g, t g'), subset(idx), conj(s), t_upper(s)."""
    Nt: int
    g0 = None                                        # g(0+) per element (None = 0)

    def t_upper(self, s):                            # upper bound of the root of g(t) = s (default: none)
        return np.full_like(s, np.inf)

    def prox(self, s, r, t0=None):
        """t >= 0 with g(t) + r t = s (s >= 0 elementwise; t = 0 where s <= g(0+))."""
        s = np.asarray(s, float)
        t = np.zeros_like(s)
        act = s > (0.0 if self.g0 is None else self.g0)
        if not act.any():
            return t
        allact = act.all()
        idx = None if allact else np.flatnonzero(act)
        F = self if allact else self.subset(idx)
        ss = s if allact else s[idx]
        thi = np.log(np.minimum(ss / r, F.t_upper(ss)))
        if t0 is None:
            tau0 = thi.copy()
        else:
            tp = t0 if allact else t0[idx]
            tau0 = np.minimum(np.log(np.maximum(tp, 1e-300)), thi)
        tr, bad = _newton_log(F, ss, r, thi, tau0)
        if bad.size or not np.all(np.isfinite(tr)):     # robust fallback: bisection of glaplace.nonlinear
            bad = np.union1d(bad, np.flatnonzero(~np.isfinite(tr)))
            self.n_fallback = getattr(self, "n_fallback", 0) + bad.size
            Fb = F.subset(bad)
            tr[bad] = solve_scalar(lambda tt, ii: Fb.subset(ii).g(tt), ss[bad], r)
        if allact:
            return tr
        t[idx] = tr
        return t

    def update_nu(self, wx, wy, r, t0=None):
        mag = np.hypot(wx, wy)
        t = self.prox(r * mag, r, t0)
        sc = np.divide(t, mag, out=np.zeros_like(mag), where=mag > 0)
        return wx * sc, wy * sc, t

    def _conj_newton(self, s):
        """Upper bound of G^*(s) for elements with s > g(0+), via the root of g(t) = s and a verified bracket."""
        thi = self.t_upper(s)
        if not np.all(np.isfinite(thi)):
            raise ValueError("t_upper needed for the conjugate")
        tau_hi = np.log(thi)
        t, bad = _newton_log(self, s, 0.0, tau_hi, tau_hi.copy())
        return _bracket_bound(self, s, t)


def _bracket_bound(F, s, t):
    """Rigorous upper bound s*hi - G(lo) of G^*(s), with lo <= t* <= hi verified (g(lo) <= s <= g(hi))."""
    lo, hi = t * (1 - 1e-13), t * (1 + 1e-13)
    for _ in range(80):
        bad_lo = F.g(lo) > s
        bad_hi = F.g(hi) < s
        if not (bad_lo.any() or bad_hi.any()):
            break
        lo = np.where(bad_lo, lo * (1 - 1e-6), lo)
        hi = np.where(bad_hi, hi * (1 + 1e-6), hi)
    else:
        raise RuntimeError("conjugate bracket not found")
    return s * hi - F.G(lo)


class PowerSum(ElementIntegrand):
    """G_T(t) = sum_k c[k][T] t^{p[k][T]} / p[k][T];  c: (K, Nt) >= 0, p: (K, Nt) > 1.
    Zero coefficients are dropped (elements with a single active term are evaluated with one power only)."""

    def __init__(self, c, p, _check=True):
        c = np.atleast_2d(np.asarray(c, float)); p = np.atleast_2d(np.asarray(p, float))
        c, p = np.broadcast_arrays(c, p)
        self.c, self.p = np.ascontiguousarray(c), np.ascontiguousarray(p)
        if _check:
            if np.any(self.p <= 1) or np.any(self.c < 0):
                raise ValueError("need p > 1, c >= 0")
            if np.any(self.c.sum(axis=0) <= 0):
                raise ValueError("every element needs a positive coefficient")
        self.K, self.Nt = self.c.shape
        self.pm1 = self.p - 1.0
        act = self.c > 0
        self.single = act.sum(axis=0) == 1
        self.allsingle = bool(self.single.all())
        if self.allsingle:                                # compact one-term representation
            k = np.argmax(act, axis=0); el = np.arange(self.Nt)
            self.c1, self.p1 = self.c[k, el], self.p[k, el]
            self.unit = bool(np.all(self.c1 == 1.0))
        self._split = None

    def subset(self, idx):
        return PowerSum(self.c[:, idx], self.p[:, idx], _check=False)

    def _parts(self):
        """Split into (single-term elements, multi-term elements) for speed."""
        if self._split is None:
            si, mi = np.flatnonzero(self.single), np.flatnonzero(~self.single)
            S = self.subset(si) if si.size else None
            M = self.subset(mi) if mi.size else None
            if M is not None:                             # drop zero rows within the multi part (keep K)
                pass
            self._split = (si, S, mi, M)
        return self._split

    # evaluation ---------------------------------------------------------------------------------------------------
    def G(self, t):
        t = np.asarray(t, float)
        if self.allsingle:
            return self.c1 * t ** self.p1 / self.p1
        return np.sum(self.c * t[None] ** self.p / self.p, axis=0)

    def g(self, t):
        t = np.asarray(t, float)
        if self.allsingle:
            v = t ** (self.p1 - 1.0)
            return v if self.unit else self.c1 * v
        return np.sum(self.c * t[None] ** self.pm1, axis=0)

    def g_tdg(self, t):
        t = np.asarray(t, float)
        if self.allsingle:
            v = t ** (self.p1 - 1.0)
            if not self.unit:
                v = self.c1 * v
            return v, (self.p1 - 1.0) * v
        T = self.c * t[None] ** self.pm1
        return T.sum(axis=0), (self.pm1 * T).sum(axis=0)

    def t_upper(self, s):
        # g(t) >= c_k t^{p_k - 1} for every k  =>  t* <= min_k (s/c_k)^{1/(p_k-1)}
        if self.allsingle:
            return (s / self.c1) ** (1.0 / (self.p1 - 1.0))
        with np.errstate(divide="ignore"):
            cand = np.where(self.c > 0, (s[None] / np.where(self.c > 0, self.c, 1.0)) ** (1.0 / self.pm1), np.inf)
        return cand.min(axis=0)

    def prox(self, s, r, t0=None):
        if self.allsingle:
            return ElementIntegrand.prox(self, s, r, t0)
        si, S, mi, M = self._parts()
        t = np.zeros_like(np.asarray(s, float))
        if S is not None:
            t[si] = S.prox(s[si], r, None if t0 is None else t0[si])
        if M is not None:
            t[mi] = ElementIntegrand.prox(M, s[mi], r, None if t0 is None else t0[mi])
        return t

    def conj(self, s):
        """Upper bound for G_T^*(s) (exact closed form for single-term elements)."""
        s = np.asarray(s, float)
        if self.allsingle:
            q = self.p1 / (self.p1 - 1.0)
            return s * (s / self.c1) ** (1.0 / (self.p1 - 1.0)) / q
        si, S, mi, M = self._parts()
        out = np.zeros_like(s)
        if S is not None:
            out[si] = S.conj(s[si])
        if M is not None:
            sm = s[mi]
            pos = np.flatnonzero(sm > 0)
            vals = np.zeros_like(sm)
            if pos.size:
                vals[pos] = M.subset(pos)._conj_newton(sm[pos])
            out[mi] = vals
        return out


class ExpMinus1MinusT(ElementIntegrand):
    """G = w (e^t - 1 - t) + beta t^2/2 with per-element weights w >= 0, beta >= 0 (w + beta > 0)."""

    def __init__(self, Nt, w=1.0, beta=0.0):
        self.Nt = Nt
        self.w = np.broadcast_to(np.asarray(w, float), (Nt,)).copy()
        self.beta = np.broadcast_to(np.asarray(beta, float), (Nt,)).copy()

    def subset(self, idx):
        return ExpMinus1MinusT(len(idx), self.w[idx], self.beta[idx])

    def t_upper(self, s):
        with np.errstate(divide="ignore"):
            a = np.where(self.w > 0, np.log1p(s / np.where(self.w > 0, self.w, 1.0)), np.inf)
            b = np.where(self.beta > 0, s / np.where(self.beta > 0, self.beta, 1.0), np.inf)
        return np.minimum(a, b)

    def G(self, t):
        with np.errstate(over="ignore"):
            return self.w * (np.expm1(t) - t) + 0.5 * self.beta * t * t

    def g(self, t):
        with np.errstate(over="ignore"):
            return self.w * np.expm1(t) + self.beta * t

    def g_tdg(self, t):
        with np.errstate(over="ignore"):
            e = np.expm1(t)
            return self.w * e + self.beta * t, t * (self.w * (e + 1.0) + self.beta)

    def conj(self, s):
        s = np.asarray(s, float)
        out = np.zeros_like(s)
        pure = (self.beta == 0) & (s > 0)
        x = s[pure] / self.w[pure]
        out[pure] = self.w[pure] * ((1 + x) * np.log1p(x) - x)
        mix = np.flatnonzero((self.beta > 0) & (s > 0))
        if mix.size:
            out[mix] = self.subset(mix)._conj_newton(s[mix])
        return out


class ExpMinus1(ElementIntegrand):
    """G = e^t - 1 (g = e^t, g(0+) = 1: the ADMM step returns t = 0 when r|w| <= 1)."""

    def __init__(self, Nt):
        self.Nt = Nt
        self.g0 = np.ones(Nt)

    def subset(self, idx):
        return ExpMinus1(len(idx))

    def t_upper(self, s):
        return np.log(np.maximum(s, 1.0))

    def G(self, t):
        with np.errstate(over="ignore"):
            return np.expm1(t)

    def g(self, t):
        with np.errstate(over="ignore"):
            return np.exp(t)

    def g_tdg(self, t):
        with np.errstate(over="ignore"):
            e = np.exp(t)
            return e, t * e

    def conj(self, s):
        s = np.asarray(s, float)
        out = np.zeros_like(s)
        m = s > 1
        out[m] = s[m] * np.log(s[m]) - s[m] + 1.0
        return out


class Generic(ElementIntegrand):
    """Wrapper around a glaplace.young.Young object sampled at (xb, yb): bisection prox (existing code) and a bisection
    conjugate with a verified bracket (needs g unbounded)."""

    def __init__(self, Gf: Young, xb, yb):
        self.Gf, self.xb, self.yb = Gf, np.asarray(xb, float), np.asarray(yb, float)
        self.Nt = len(self.xb)
        self.g0 = np.asarray(Gf.g(np.zeros(self.Nt), self.xb, self.yb), float)

    def subset(self, idx):
        return Generic(self.Gf, self.xb[idx], self.yb[idx])

    def G(self, t):
        return self.Gf.G(t, self.xb, self.yb)

    def g(self, t):
        return self.Gf.g(t, self.xb, self.yb)

    def prox(self, s, r, t0=None):
        return solve_scalar(lambda tt, ii: self.Gf.g(tt, self.xb[ii], self.yb[ii]), s, r)

    def conj(self, s, iters=200):
        s = np.asarray(s, float)
        lo, hi = np.zeros_like(s), np.maximum(1.0, s)
        for _ in range(200):
            bad = self.g(hi) < s
            if not bad.any():
                break
            hi = np.where(bad, 2 * hi, hi)
        for _ in range(iters):
            m = 0.5 * (lo + hi)
            neg = self.g(m) < s
            lo, hi = np.where(neg, m, lo), np.where(neg, hi, m)
        out = s * hi - self.G(lo)
        return np.where(s <= self.g0, 0.0, out)


# ------------------------------------------------------------------------------------------------------------------
def from_young(Gf, xb, yb):
    """Convert an existing glaplace Young object (PowerLaw, VarPower, Sum of weighted PowerLaws, Exp, ExpLiteral) to the
    fastest element-wise representation; anything else -> Generic (bisection)."""
    from .cr import VarPower
    xb, yb = np.asarray(xb, float), np.asarray(yb, float)
    Nt = len(xb)
    if isinstance(Gf, PowerLaw):
        return PowerSum(np.ones((1, Nt)), np.full((1, Nt), Gf.p))
    if isinstance(Gf, VarPower):
        return PowerSum(np.ones((1, Nt)), Gf._p(xb, yb)[None])
    if isinstance(Gf, Sum) and all(isinstance(G, PowerLaw) for _, G in Gf.terms):
        c = np.vstack([np.broadcast_to(np.asarray(w(xb, yb), float) if callable(w) else float(w), (Nt,)) for w, _ in Gf.terms])
        p = np.vstack([np.full(Nt, G.p) for _, G in Gf.terms])
        return PowerSum(c, p)
    if type(Gf) is Exp:
        return ExpMinus1MinusT(Nt)
    if type(Gf) is ExpLiteral:
        return ExpMinus1(Nt)
    return Generic(Gf, xb, yb)
