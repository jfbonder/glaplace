"""Tools for the v4 experiments (NEW module): disk meshes resolving circles, quadrature adapted to a point singularity at
the origin (loads f = |x|^-beta and true errors), CR evaluation at quadrature points, element-wise min / max of a
coefficient (minorant / majorant one-point quadrature), Lyapunov functional.
"""
from __future__ import annotations
import numpy as np
from scipy.special import roots_jacobi, roots_legendre

from .mesh import Mesh, boundary_nodes_from_triangles


# ------------------------------------------------------------------------------------------------------------------
def disk_mesh_resolved(h, circles=(0.5,), origin=True, min_angle=30.0):
    """Quality mesh (Shewchuk's Triangle) of the polygon inscribed in the unit circle, with vertices on the circles
    r = c for c in ``circles`` (segments are kept without Steiner points, so every triangle lies on one side of the
    polygonal interface) and a vertex at the origin.  h ~ target edge length."""
    import triangle as tr
    pts, segs = [], []
    off = 0
    for R in (1.0,) + tuple(circles):
        k = max(8, int(np.ceil(2 * np.pi * R / h)))
        th = 2 * np.pi * np.arange(k) / k
        pts.append(np.column_stack([R * np.cos(th), R * np.sin(th)]))
        segs.append(np.column_stack([off + np.arange(k), off + (np.arange(k) + 1) % k]))
        off += k
    if origin:
        pts.append(np.zeros((1, 2)))
    A = dict(vertices=np.vstack(pts), segments=np.vstack(segs))
    B = tr.triangulate(A, f"pq{min_angle:g}Ya{h * h * np.sqrt(3) / 4:.12f}")   # fixed-point: "e" would be read as a switch
    T = B["triangles"]
    return Mesh(B["vertices"], T, boundary_nodes_from_triangles(T))


def mesh_h(mesh):
    V, T = mesh.vertices, mesh.triangles
    e = np.concatenate([np.hypot(*(V[T[:, i]] - V[T[:, (i + 1) % 3]]).T) for i in range(3)])
    return float(e.max())


# ------------------------------------------------------------------------------------------------------------------
class SingularQuadrature:
    """Per-element quadrature  int_T F dx ~ sum_q W[T,q] F(X[T,q]).

    Every triangle uses a collapsed (Duffy) Gauss product rule with the collapse at vertex 0 (n x n points, exact for
    polynomials of degree 2n-1 in (s,t)).  On triangles having a vertex at ``origin`` the collapse is at that vertex and the
    radial direction uses Gauss-Jacobi with weight s^(1-gamma), so that F ~ |x|^-gamma (times smooth) is integrated
    accurately.  Also stores barycentric coordinates LAM[T,q,i] for evaluating P1 / CR functions."""

    def __init__(self, mesh, gamma=0.0, n=8, n_sing=16, origin=(0.0, 0.0)):
        V, T = mesh.vertices, mesh.triangles.copy()
        Nt = len(T)
        o = np.asarray(origin, float)
        at0 = np.linalg.norm(V[T] - o, axis=2) < 1e-14                   # (Nt,3)
        sing = at0.any(axis=1)
        # rotate local numbering so that the singular vertex is local vertex 0 (keeps orientation)
        perm = np.tile(np.arange(3), (Nt, 1))
        k = np.argmax(at0, axis=1)
        perm[sing] = (np.arange(3)[None, :] + k[sing, None]) % 3
        self.perm = perm
        Tp = np.take_along_axis(T, perm, axis=1)
        P = V[Tp]                                                         # (Nt,3,2)
        det = np.abs((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
        nq = max(n, n_sing) ** 2
        self.W = np.zeros((Nt, nq))
        self.LAMp = np.full((Nt, nq, 3), 1.0 / 3.0)              # unused slots: barycentre with weight 0
        self.X = np.repeat(P.mean(axis=1)[:, None, :], nq, axis=1)
        for flag, m in ((False, n), (True, n_sing)):
            idx = np.flatnonzero(sing == flag)
            if idx.size == 0:
                continue
            xt, wt = roots_legendre(m); t = 0.5 * (xt + 1); wt = 0.5 * wt
            if flag:
                xs, ws = roots_jacobi(m, 0.0, 1.0 - gamma)               # weight (1+x)^(1-gamma) on [-1,1]
                s = 0.5 * (xs + 1); ws = ws * 0.5 ** (2.0 - gamma)       # -> weight s^(1-gamma) on [0,1]
                wfac = s ** gamma                                        # s^gamma * F is smooth
            else:
                xs, ws = roots_legendre(m); s = 0.5 * (xs + 1); ws = 0.5 * ws * s; wfac = np.ones(m)
            S, Tt = np.meshgrid(s, t, indexing="ij")
            Wq = (ws * wfac)[:, None] * wt[None, :]
            l0 = (1 - S).ravel(); l1 = (S * (1 - Tt)).ravel(); l2 = (S * Tt).ravel()
            L = np.stack([l0, l1, l2], axis=1)                            # (m*m, 3)
            q = m * m
            self.LAMp[idx, :q] = L[None]
            self.X[idx, :q] = np.einsum("qi,tid->tqd", L, P[idx])
            self.W[idx, :q] = det[idx, None] * Wq.ravel()[None]
        # barycentric coordinates in the ORIGINAL local numbering of the mesh triangles
        self.LAM = np.zeros_like(self.LAMp)
        np.put_along_axis(self.LAM, np.repeat(perm[:, None, :], nq, axis=1), self.LAMp, axis=2)
        self.sing = sing

    def integrate(self, F):
        """F: array (Nt, nq) of values at X -> per-element integrals."""
        return np.sum(self.W * F, axis=1)

    def eval_P1(self, mesh, u):
        return np.einsum("tqi,ti->tq", self.LAM, u[mesh.triangles])

    def eval_CR(self, fem, u):
        return np.einsum("tqi,ti->tq", 1.0 - 2.0 * self.LAM, u[fem.t2e])

    def load_P1(self, mesh, f):
        """b_i = int f phi_i (phi_i = lambda_i)."""
        fv = f(self.X[..., 0], self.X[..., 1]) * self.W
        b = np.zeros(mesh.n_nodes)
        np.add.at(b, mesh.triangles.ravel(), np.einsum("tq,tqi->ti", fv, self.LAM).ravel())
        return b

    def load_CR(self, fem, f):
        """b_E = int f phi_E (phi_E = 1 - 2 lambda_i on T, i the vertex opposite E)."""
        fv = f(self.X[..., 0], self.X[..., 1]) * self.W
        b = np.zeros(fem.Ne)
        np.add.at(b, fem.t2e.ravel(), np.einsum("tq,tqi->ti", fv, 1.0 - 2.0 * self.LAM).ravel())
        return b


def barycentre_load(fem, f, scheme):
    """One-point barycentre rule: b = sum_T f(x_T) |T| / 3 for each of the three dofs of T (P1 nodes or CR edges)."""
    fb = f(fem.bary[:, 0], fem.bary[:, 1]) * fem.area / 3.0
    if scheme == "P1":
        b = np.zeros(fem.Nn); np.add.at(b, fem.mesh.triangles.ravel(), np.repeat(fb, 3)); return b
    b = np.zeros(fem.Ne); np.add.at(b, fem.t2e.ravel(), np.repeat(fb, 3)); return b


# ------------------------------------------------------------------------------------------------------------------
def element_minmax(fun, mesh, m=64, levels=4):
    """min_T fun and max_T fun for every triangle, assuming the extrema of fun over T lie on the boundary of T (true for
    a(x) = max(x2^2 - x1^2, 0)/|x|^(2-alpha) with alpha = 1: 1-homogeneous with nonvanishing gradient where a > 0).
    Each edge is sampled at m+1 points; the best sample is refined by resampling its two neighbouring sub-intervals
    (``levels`` times), i.e. the extremum is located to ~m^-(levels+1) of the edge length."""
    P = mesh.vertices[mesh.triangles]
    Nt = len(P)
    mins, maxs = np.full(Nt, np.inf), np.full(Nt, -np.inf)
    for i in range(3):
        a, b = P[:, i], P[:, (i + 1) % 3]
        for sign, out in ((1.0, mins), (-1.0, maxs)):
            lo, hi = np.zeros(Nt), np.ones(Nt)
            best = np.full(Nt, np.inf)
            for _ in range(levels + 1):
                s = lo[:, None] + (hi - lo)[:, None] * np.linspace(0, 1, m + 1)[None]
                X = a[:, None, :] + s[..., None] * (b - a)[:, None, :]
                F = sign * fun(X[..., 0], X[..., 1])
                k = np.argmin(F, axis=1)
                best = np.minimum(best, F[np.arange(Nt), k])
                dl = (hi - lo) / m
                lo, hi = np.maximum(s[np.arange(Nt), k] - dl, 0.0), np.minimum(s[np.arange(Nt), k] + dl, 1.0)
            if sign > 0:
                np.minimum(out, best, out=out)
            else:
                np.maximum(out, -best, out=out)
    return mins, maxs


def a_elm(alpha=1.0):
    """Esposito-Leonetti-Mingione weight a(x) = max(x2^2 - x1^2, 0) / |x|^(2-alpha)  (a(0) = 0)."""
    def a(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        r = np.hypot(x, y)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.maximum(y * y - x * x, 0.0) / r ** (2.0 - alpha)
        return np.where(r > 0, v, 0.0)
    return a


def lyapunov(area, r, nu, eta_next, nu_star, eta_star):
    """E_n = ||eta_{n+1} - eta*||^2/(2r) + r ||nu_n - nu*||^2 / 2 (Theorem A / C)."""
    de = (eta_next[0] - eta_star[0]) ** 2 + (eta_next[1] - eta_star[1]) ** 2
    dn = (nu[0] - nu_star[0]) ** 2 + (nu[1] - nu_star[1]) ** 2
    return float(np.sum(area * de) / (2 * r) + r * np.sum(area * dn) / 2)
