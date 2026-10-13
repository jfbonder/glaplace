"""Test problems with a Lavrentiev gap (Zhikov checkerboard) and a driver comparing P1 (conforming) with CR (nonconforming).
NEW module.  Setting (Zhikov 1986/1995; Balci-Ortner-Storn 2022, Sect. 3-4): Omega = (-1,1)^2, f = 0, Dirichlet data lambda*u0 with
    u0 = 1 (|x1| < x2),  -1 (|x1| < -x2),  x2/|x1| (otherwise)       (homogeneous of degree 0: singular at the origin),
and one of
    p(x)-Laplacian:   G = t^{p(x)}/p(x),  p = p_+ for |x1| < |x2|,  p_- otherwise  (1 < p_- < 2 < p_+),
    two-phase:        G = t^p/p + a(x) t^q/q,  a = 1 for |x1| < |x2|, 0 otherwise  (discontinuous, alpha = 0).
For lambda = 1 (p(x) case) the minimiser is u = x2 (its flux (0,1) is divergence free in both regions) and F(u) = 2/p_- + 2/p_+ = 2 for (3/2, 3).
For other lambda no closed form is used here; for lambda > 1 a Lavrentiev gap is reported in the literature (Zhikov; Balci-Ortner-Storn).
"""
from __future__ import annotations
import numpy as np
from .mesh import Mesh
from .young import pq
from .cr import VarPower, GLaplaceCRSolver
from .solver import GLaplaceSolver


def crisscross_mesh(n: int, L: float = 1.0) -> Mesh:
    """Uniform mesh of (-L,L)^2 with n x n squares (n even), each split into 4 triangles through its centre.
    Both lines |x1| = |x2| consist of mesh edges, so p(x), a(x) are constant on every triangle."""
    if n % 2:
        raise ValueError("n must be even")
    xs = np.linspace(-L, L, n + 1)
    X, Y = np.meshgrid(xs, xs)
    V = np.column_stack([X.ravel(), Y.ravel()])
    idx = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    xc = 0.5 * (xs[:-1] + xs[1:])
    XC, YC = np.meshgrid(xc, xc)
    C0 = len(V)
    V = np.vstack([V, np.column_stack([XC.ravel(), YC.ravel()])])
    cidx = C0 + np.arange(n * n).reshape(n, n)
    a, b = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel()
    c, d = idx[1:, 1:].ravel(), idx[1:, :-1].ravel()
    m = cidx.ravel()
    T = np.vstack([np.column_stack([a, b, m]), np.column_stack([b, c, m]),
                   np.column_stack([c, d, m]), np.column_stack([d, a, m])])
    bd = np.unique(np.concatenate([idx[0, :], idx[-1, :], idx[:, 0], idx[:, -1]]))
    return Mesh(V, T, bd)


def u0(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ax = np.abs(x)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(ax > 0, y / np.where(ax > 0, ax, 1.0), 0.0)
    return np.where(ax < y, 1.0, np.where(ax < -y, -1.0, r))


def psi(lam):
    return lambda x, y: lam * u0(x, y)


def p_zhikov(pm=1.5, pp=3.0):
    return lambda x, y: np.where(np.abs(x) < np.abs(y), pp, pm)


def a_checker(x, y):
    return np.where(np.abs(x) < np.abs(y), 1.0, 0.0)


def G_pxlap(pm=1.5, pp=3.0):
    return VarPower(p_zhikov(pm, pp))


def G_twophase(p=1.5, q=3.0):
    return pq(p, q, a_checker)


def energy_exact_pxlap_lambda1(pm=1.5, pp=3.0):
    """Exact minimum for lambda = 1 (u = x2):  F = 2/pm + 2/pp."""
    return 2 / pm + 2 / pp


def compare(G, lam, ns=(8, 16, 32, 64), r=3.0, tol=1e-8, max_iter=20000, verbose=True):
    """Solve with P1 (conforming) and CR (nonconforming) on criss-cross meshes; returns a list of dict rows.
    Energies are F_h(u_h) = sum_T |T| G(x_T, grad u_h|_T) (f = 0), the same formula for both schemes."""
    rows = []
    for n in ns:
        mesh = crisscross_mesh(n)
        for name, Solver in (("P1", GLaplaceSolver), ("CR", GLaplaceCRSolver)):
            S = Solver(mesh, G, 0.0, psi(lam), r=r)
            res = S.solve(tol=tol, max_iter=max_iter, stop="residual")
            row = dict(scheme=name, n=n, ndof=S.fem.Nn, iters=res.iterations, converged=res.converged,
                       energy=S.energy(res.u), residual=res.history["residual"][-1], time=res.time_iter)
            rows.append(row)
            if verbose:
                print(f"  lam={lam:g} n={n:3d} {name} ndof={row['ndof']:6d} it={row['iters']:5d} conv={row['converged']!s:5} "
                      f"F_h={row['energy']:.6f} res={row['residual']:.1e} ({row['time']:.1f}s)", flush=True)
    return rows
