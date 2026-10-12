"""Young functions  G(x, t)  and their derivatives  g(x, t) = d/dt G(x, t),  t >= 0.

The g-Laplacian is  -div( g(x,|grad u|) grad u/|grad u| ) = f.  A Young function object provides two
vectorised methods

    Gf.G(t, x, y)   -> G(x,y,t)
    Gf.g(t, x, y)   -> g(x,y,t) = dG/dt          (g(0) is understood as the limit g(0+))

where ``t`` is an array of |grad u| values (one per triangle) and ``x, y`` the barycentres of the
triangles (arrays of the same shape; they are only used by x-dependent functions).

Building blocks (all with G(0) = 0, G convex in t) and the algebra:

    PowerLaw(p)            t^p / p                          g = t^(p-1)
    PLog(p, a)             t^p log^a(1+t)                   (p-Laplacian with a logarithmic perturbation)
    Exp()                  e^t - 1 - t                      g = e^t - 1        (does NOT satisfy Delta_2)
    ExpLiteral()           e^t - 1                          g = e^t, g(0+) = 1 (flat where |flux| <= 1: "dead core")
    Custom(G, g)           any pair of callables G(t), g(t)
    from_expr("t**3/3")    G given as a sympy expression in t (needs sympy); g = dG/dt is computed symbolically

    G1 + G2, c*G, w(x,y)*G  (sums with constant or spatially varying non-negative weights), e.g.
    ``PowerLaw(2) + a * PowerLaw(4)`` with ``a = lambda x, y: np.where(x > 0, 1.0, 0.0)`` is a two-phase problem.

``check_young(Gf)`` tests numerically that G(0)=0, g >= 0 and g is non-decreasing (convexity of G).
"""
from __future__ import annotations
import numpy as np


def _arr(t):
    return np.asarray(t, float)


class Young:
    """Base class; subclasses implement ``_G(t)`` and ``_g(t)`` (x-independent) or override ``G``/``g``."""
    __array_ufunc__ = None          # so that  numpy_scalar * Young  calls Young.__rmul__

    def G(self, t, x=None, y=None):
        return self._G(_arr(t))

    def g(self, t, x=None, y=None):
        return self._g(_arr(t))

    def __add__(self, other):
        return Sum([(1.0, self)]).add(other)

    __radd__ = __add__

    def __rmul__(self, w):
        return Sum([(w, self)])

    __mul__ = __rmul__

    def __repr__(self):
        return self.__class__.__name__ + "()"


class PowerLaw(Young):
    """G(t) = t^p / p,  p > 1   (the p-Laplacian)."""

    def __init__(self, p):
        if p <= 1:
            raise ValueError("p must be > 1")
        self.p = float(p)

    def _G(self, t):
        return t ** self.p / self.p

    def _g(self, t):
        return t ** (self.p - 1.0)

    def __repr__(self):
        return f"PowerLaw(p={self.p})"


class PLog(Young):
    """G(t) = t^p log^a(1+t)  (p > 1; need p + a > 1 and G convex, see ``check_young``)."""

    def __init__(self, p, a):
        self.p, self.a = float(p), float(a)
        if self.p + self.a <= 1:
            raise ValueError("need p + a > 1")

    def _G(self, t):
        t = np.maximum(t, 1e-300)
        return t ** (self.p + self.a) * (np.log1p(t) / t) ** self.a

    def _g(self, t):
        # g = t^(p+a-1) [ p lt^a + a lt^(a-1)/(1+t) ],  lt = log(1+t)/t  (stable for tiny t; lt -> 1)
        t = np.maximum(t, 1e-300)
        lt = np.log1p(t) / t
        return t ** (self.p + self.a - 1.0) * (self.p * lt ** self.a + self.a * lt ** (self.a - 1.0) / (1.0 + t))

    def __repr__(self):
        return f"PLog(p={self.p}, a={self.a})"


class Exp(Young):
    """G(t) = e^t - 1 - t, g(t) = e^t - 1  (fails the Delta_2 condition)."""

    def _G(self, t):
        with np.errstate(over="ignore"):
            return np.expm1(t) - t

    def _g(self, t):
        with np.errstate(over="ignore"):
            return np.expm1(t)


class ExpLiteral(Young):
    """G(t) = e^t - 1,  g(t) = e^t, g(0+) = 1 > 0.  Convex, G(0)=0, no Delta_2.  Because g(0+) > 0 the 'flux'
    g(|grad u|) grad u/|grad u| has a jump at grad u = 0 (like the 1-Laplacian): the solution is flat where
    |flux| <= 1 (dead core, e.g. torsion with large load)."""

    def _G(self, t):
        with np.errstate(over="ignore"):
            return np.expm1(t)

    def _g(self, t):
        with np.errstate(over="ignore"):
            return np.exp(t)


class Custom(Young):
    """Young function from two callables:  Custom(G=lambda t: ..., g=lambda t: ...)."""

    def __init__(self, G, g):
        self._Gf, self._gf = G, g

    def _G(self, t):
        return _arr(self._Gf(t))

    def _g(self, t):
        return _arr(self._gf(t))


def from_expr(expr: str, var: str = "t") -> Custom:
    """Young function from a sympy expression for G(t), e.g. ``from_expr("t**3/3 + t**4/4")``."""
    import sympy as sp
    t = sp.Symbol(var, nonnegative=True)
    G = sp.sympify(expr, locals={var: t})
    g = sp.diff(G, t)
    return Custom(sp.lambdify(t, G, "numpy"), sp.lambdify(t, g, "numpy"))


class Sum(Young):
    """sum_i w_i(x,y) G_i(t),  weights w_i >= 0 constants or callables w(x, y)."""

    def __init__(self, terms):
        self.terms = list(terms)

    def add(self, other):
        if isinstance(other, Sum):
            return Sum(self.terms + other.terms)
        if isinstance(other, Young):
            return Sum(self.terms + [(1.0, other)])
        raise TypeError("can only add Young functions")

    def __add__(self, other):
        return self.add(other)

    __radd__ = __add__

    def __rmul__(self, c):
        return Sum([(_scale(c, w), Gf) for w, Gf in self.terms])

    __mul__ = __rmul__

    @staticmethod
    def _w(w, x, y, shape):
        if callable(w):
            return np.broadcast_to(np.asarray(w(x, y), float), shape)
        return float(w)

    def G(self, t, x=None, y=None):
        t = _arr(t)
        return sum(self._w(w, x, y, t.shape) * Gf.G(t, x, y) for w, Gf in self.terms)

    def g(self, t, x=None, y=None):
        t = _arr(t)
        return sum(self._w(w, x, y, t.shape) * Gf.g(t, x, y) for w, Gf in self.terms)

    def __repr__(self):
        return " + ".join(f"{'w(x)' if callable(w) else w}*{Gf!r}" for w, Gf in self.terms)


def _scale(c, w):
    if callable(c) or callable(w):
        cc = c if callable(c) else (lambda x, y, c=c: c)
        ww = w if callable(w) else (lambda x, y, w=w: w)
        return lambda x, y: cc(x, y) * ww(x, y)
    return c * w


def pq(p, q, a=1.0):
    """The p-q-Laplacian  t^p/p + a(x) t^q/q  (a >= 0 constant or callable a(x, y); a = 0 on a region gives the two-phase problem)."""
    return Sum([(1.0, PowerLaw(p)), (a, PowerLaw(q))])


def check_young(Gf, tmax=20.0, n=2001, x=0.0, y=0.0, verbose=True):
    """Numerical sanity check on a grid: G(0)=0, g >= 0, g non-decreasing (=> G convex), g = G' (finite differences).
    ``x, y`` : point at which an x-dependent G is tested.  Returns True if all checks pass."""
    t = np.linspace(0, tmax, n)
    X, Y = np.full_like(t, x), np.full_like(t, y)
    G, g = Gf.G(t, X, Y), Gf.g(t, X, Y)
    ok = {"G(0)=0": abs(G[0]) < 1e-12, "g>=0": bool(np.all(g >= -1e-12)),
          "g nondecreasing": bool(np.all(np.diff(g) >= -1e-9 * (1 + np.abs(g[1:])))),
          "g = G'": bool(np.allclose((np.diff(G) / np.diff(t))[2:], 0.5 * (g[1:] + g[:-1])[2:], rtol=1e-2, atol=1e-3))}
    if verbose:
        print(f"{Gf!r}: " + ", ".join(f"{k}: {'ok' if v else 'FAIL'}" for k, v in ok.items()))
    return all(ok.values())
