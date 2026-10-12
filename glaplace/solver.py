"""Decomposition-coordination (augmented Lagrangian / ADMM) solver for the g-Laplacian

        -div( g(x,|grad u|) grad u / |grad u| ) = f  in Omega,      u = u_D on dOmega,          g = dG/dt,

with G(x,t) a Young function (convex, G(0) = 0), P1 finite elements on any 2D triangulation.
It minimises  J(u) = int G(x,|grad u|) - int f u.

Iteration (r > 0 = augmented-Lagrangian parameter, nu_0 = eta_0 = 0):
    1. r K u_n = F + int (r nu_{n-1} - eta_n) . grad(phi)                  linear, same matrix every step
    2. nu_n = t w/|w|,  w = eta_n/r + grad u_n,   g(x,t) + r t = r|w|       pointwise scalar equation
    3. eta_{n+1} = eta_n + r (grad u_n - nu_n)       (then eta_{n+1} = g(x,|nu_n|) nu_n/|nu_n| = the flux)
For G = t^p/p this is the method for the p-Laplacian (Glowinski); here G is arbitrary.
G(x,t) is sampled at triangle barycentres (piecewise-constant dependence on x).
"""
from __future__ import annotations
from dataclasses import dataclass, field
import time
import numpy as np
import scipy.sparse.linalg as spla

from .mesh import Mesh
from .fem import P1Space
from .nonlinear import update_nu
from .young import Young, PowerLaw


@dataclass
class SolveResult:
    u: np.ndarray                 # nodal values
    nu: tuple                     # (nu_x, nu_y) per element
    eta: tuple                    # (eta_x, eta_y) per element  (= flux g(|nu|) nu/|nu| at convergence)
    iterations: int
    converged: bool
    history: dict = field(default_factory=dict)   # per-iteration diagnostics
    time_setup: float = 0.0
    time_iter: float = 0.0


class GLaplaceSolver:
    """Set-up once (mesh, G, f, u_D, factorisation), then call :meth:`solve`.

    Parameters
    ----------
    mesh : Mesh  (any triangulation; boundary nodes in ``mesh.dirichlet``)
    G : a :class:`glaplace.young.Young` object (or a number p > 1, meaning G = t^p/p)
    f, u_D : constant, array, or callable ``fun(x, y)``.  f: F_i = int f phi_i.  u_D: Dirichlet data at the boundary nodes
             (or nodal values of length Nn)
    r : augmented-Lagrangian parameter
    f_degree : quadrature degree of the load vector (1 = barycentre rule; exact only for constant f)
    """

    def __init__(self, mesh: Mesh, G, f=0.0, u_D=0.0, r: float = 1.0, f_degree: int = 1):
        t0 = time.perf_counter()
        self.mesh, self.r = mesh, float(r)
        if not isinstance(G, Young):
            G = PowerLaw(float(G))
        self.G = G
        self.fem = P1Space(mesh)
        self.xb, self.yb = self.fem.bary[:, 0], self.fem.bary[:, 1]
        self.F = self.fem.load_vector(f, degree=f_degree)
        D = mesh.dirichlet
        if np.ndim(u_D) == 1 and not callable(u_D) and len(u_D) == mesh.n_nodes:
            self.gD = np.asarray(u_D, float)[D]
        else:
            self.gD = self.fem.node_values(u_D, D)
        mask = np.ones(mesh.n_nodes, bool)
        mask[D] = False
        self.free = np.flatnonzero(mask)
        K = self.r * self.fem.K
        self.Kff = K[self.free][:, self.free].tocsc()
        self.Kfd = K[self.free][:, D].tocsr()
        self._lu = spla.splu(self.Kff)            # factorised ONCE
        self.time_setup = time.perf_counter() - t0

    # element-wise functions closed over barycentres -----------------------------------------
    def _gfun(self, t, idx=None):
        x, y = (self.xb, self.yb) if idx is None else (self.xb[idx], self.yb[idx])
        return self.G.g(t, x, y)

    def energy(self, u):
        """J_h(u) = sum_T |T| G(x_T, |grad u|_T) - F.u"""
        gx, gy = self.fem.gradient(u)
        return float(np.sum(self.fem.area * self.G.G(np.hypot(gx, gy), self.xb, self.yb)) - self.F @ u)

    def _solve_linear(self, rhs):
        u = np.empty(self.mesh.n_nodes)
        u[self.mesh.dirichlet] = self.gD
        u[self.free] = self._lu.solve(rhs[self.free] - self.Kfd @ self.gD)
        return u

    def solve(self, tol: float = 1e-6, max_iter: int = 2000, stop: str = "rel_u",
              init: str = "zero", seed: int = 0, u0=None, callback=None, verbose: bool = False) -> SolveResult:
        """Run the iteration.

        stop : 'rel_u'    -> ||u_n - u_{n-1}||_2 / ||u_{n-1}||_2 < tol   (relative change; default)
               'residual' -> ||grad u_n - nu_n||_{L2} / (1 + ||grad u_n||_{L2}) < tol   (primal residual, scale-free)
        init : 'zero' (nu = eta = 0) or 'random' (0.15*rand)
        u0   : not used to start the iteration (the first step only needs nu, eta); kept for the first rel_u comparison.
        """
        if stop not in ("rel_u", "residual"):
            raise ValueError("stop must be 'rel_u' or 'residual'")
        t0 = time.perf_counter()
        fem, r, Nt = self.fem, self.r, self.fem.Nt
        if init == "random":
            rng = np.random.default_rng(seed)
            nux, nuy = 0.15 * rng.random(Nt), 0.15 * rng.random(Nt)
            etax, etay = 0.15 * rng.random(Nt), 0.15 * rng.random(Nt)
        else:
            nux = nuy = etax = etay = np.zeros(Nt)
        hist = {k: [] for k in ("rel_u", "residual", "energy")}
        u_old = np.zeros(self.mesh.n_nodes) if u0 is None else np.asarray(u0, float).copy()
        u_old[self.mesh.dirichlet] = self.gD
        converged, it = False, 0
        for it in range(1, max_iter + 1):
            rhs = self.F + fem.weak_div(r * nux - etax, r * nuy - etay)
            u = self._solve_linear(rhs)
            gx, gy = fem.gradient(u)
            nux, nuy = update_nu(self._gfun, etax / r + gx, etay / r + gy, r)
            etax = etax + r * (gx - nux)
            etay = etay + r * (gy - nuy)
            du = u - u_old
            nrm = np.linalg.norm(u_old)
            rel = np.linalg.norm(du) / nrm if nrm > 0 else np.inf
            res = float(np.sqrt(np.sum(fem.area * ((gx - nux) ** 2 + (gy - nuy) ** 2))))
            gn = float(np.sqrt(np.sum(fem.area * (gx ** 2 + gy ** 2))))
            hist["rel_u"].append(rel)
            hist["residual"].append(res / (1.0 + gn))
            hist["energy"].append(float(np.sum(fem.area * self.G.G(np.hypot(gx, gy), self.xb, self.yb)) - self.F @ u))
            if verbose and (it % 50 == 0 or it == 1):
                print(f"it {it:5d}  rel_u={rel:.3e}  residual={hist['residual'][-1]:.3e}  E={hist['energy'][-1]:.8f}")
            if callback is not None:
                callback(it, u, hist)
            u_old = u
            if (rel if stop == "rel_u" else hist["residual"][-1]) < tol:
                converged = True
                break
        return SolveResult(u=u, nu=(nux, nuy), eta=(etax, etay), iterations=it, converged=converged,
                           history={k: np.array(v) for k, v in hist.items()},
                           time_setup=self.time_setup, time_iter=time.perf_counter() - t0)


def solve(mesh: Mesh, G, f=0.0, u_D=0.0, r=1.0, **kw) -> SolveResult:
    """One-shot wrapper:  solve(mesh, G, f, u_D, tol=..., stop=...)."""
    ik = {k: kw.pop(k) for k in list(kw) if k in {"f_degree"}}
    return GLaplaceSolver(mesh, G, f, u_D, r=r, **ik).solve(**kw)
