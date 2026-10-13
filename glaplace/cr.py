"""Crouzeix-Raviart (CR) variant of the decomposition-coordination solver (NEW module; nothing else is modified).

Why: for integrands with a Lavrentiev gap (e.g. t^p/p + a(x) t^q/q, p < 2 < q, a discontinuous; or p(x) piecewise constant,
p_- < 2 < p_+) every CONFORMING space (P1, P2, ...) consists of Lipschitz functions and its discrete minima converge to the
infimum over Lipschitz (smooth) functions, which is larger than the true minimum over W^{1,1} (Carstensen-Ortner 2010).
The nonconforming CR element (Balci-Ortner-Storn 2022) converges to the W-minimiser for suitable one-point quadrature of x.

Discrete space.  CR_h = piecewise affine functions, continuous at the MIDPOINTS of edges; dofs = values at edge midpoints;
u = 0 (or the edge mean of u_D) at the midpoints of boundary edges.  The broken gradient is CONSTANT per triangle, exactly like
the P1 gradient, so the multiplier nu, eta (piecewise constants) and the pointwise nonlinear step are IDENTICAL to the P1 code.
Only the space of u changes (basis phi_i = 1 - 2 lambda_i on each triangle, i = edge opposite vertex i).
The energy is  F_h(v) = sum_T |T| G(x_T, grad_h v|_T) - int f v,  x_T = barycentre (one-point rule, as in the P1 code).

Because the iteration is the same finite-dimensional ADMM, convergence for every r > 0 and each h is Theorem C of the paper
(coercivity on CR holds: discrete Poincare inequality for CR functions with zero edge means on the boundary).

Usage:
    from glaplace.cr import GLaplaceCRSolver
    S = GLaplaceCRSolver(mesh, G, f, u_D, r=...);  res = S.solve(tol=1e-8, stop='residual')   # same API as GLaplaceSolver
    S.fem.element_values(res.u)  -> value at the barycentres;   S.energy(res.u) -> F_h(u)
"""
from __future__ import annotations
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .mesh import Mesh
from .fem import P1Space, _eval
from .solver import GLaplaceSolver
from .young import Young, PowerLaw


class CRSpace(P1Space):
    """Crouzeix-Raviart space on ``mesh`` (reuses the P1Space interface: area, bary, Gx, Gy, K, gradient, weak_div, errors)."""

    def __init__(self, mesh: Mesh):
        super().__init__(mesh)                                   # area, bary, gx, gy (= grad lambda_i), Nt
        T = mesh.triangles
        E = np.empty((self.Nt, 3, 2), dtype=np.int64)            # local edge i is opposite to vertex i
        for i in range(3):
            E[:, i, 0], E[:, i, 1] = T[:, (i + 1) % 3], T[:, (i + 2) % 3]
        E.sort(axis=2)
        uniq, inv, counts = np.unique(E.reshape(-1, 2), axis=0, return_inverse=True, return_counts=True)
        self.edges = uniq                                        # (Ne, 2) vertex indices
        self.t2e = inv.reshape(self.Nt, 3)                       # (Nt, 3) global edge of the local edge i
        self.Ne = len(uniq)
        self.Nn = self.Ne                                        # number of dofs (attribute reused by the solver)
        self.bdry_edges = np.flatnonzero(counts == 1)
        self.mid = mesh.vertices[uniq].mean(axis=1)              # (Ne, 2) edge midpoints
        # grad phi_i = -2 grad lambda_i  (constant per triangle)
        rows = np.repeat(np.arange(self.Nt), 3)
        cols = self.t2e.ravel()
        self.cgx, self.cgy = -2.0 * self.gx, -2.0 * self.gy
        self.Gx = sp.csr_matrix((self.cgx.ravel(), (rows, cols)), shape=(self.Nt, self.Ne))
        self.Gy = sp.csr_matrix((self.cgy.ravel(), (rows, cols)), shape=(self.Nt, self.Ne))
        W = sp.diags(self.area)
        self.K = (self.Gx.T @ W @ self.Gx + self.Gy.T @ W @ self.Gy).tocsr()

    def load_vector(self, f, degree: int = 2) -> np.ndarray:
        """b_e = int f phi_e with the edge-midpoint rule (exact for constant f; weights |T|/3, phi_i(mid_k) = delta_ik)."""
        P = self.mesh.vertices[self.mesh.triangles]
        mids = [0.5 * (P[:, 1] + P[:, 2]), 0.5 * (P[:, 2] + P[:, 0]), 0.5 * (P[:, 0] + P[:, 1])]
        fm = np.column_stack([_eval(f, m) for m in mids])
        b = np.zeros(self.Ne)
        np.add.at(b, self.t2e.ravel(), (self.area[:, None] / 3.0 * fm).ravel())
        return b

    def edge_values(self, fun_or_const, idx) -> np.ndarray:
        """Mean over the edges ``idx`` (Simpson rule) of a function / constant; used for Dirichlet data."""
        V = self.mesh.vertices
        a, b = V[self.edges[idx, 0]], V[self.edges[idx, 1]]
        m = 0.5 * (a + b)
        return (_eval(fun_or_const, a) + 4 * _eval(fun_or_const, m) + _eval(fun_or_const, b)) / 6.0

    def element_values(self, u) -> np.ndarray:
        """Value of the CR function at the barycentre of every triangle (= mean of the three edge values)."""
        return u[self.t2e].mean(axis=1)

    def nodal_average(self, u) -> np.ndarray:
        """Vertex values obtained by averaging the CR function over the triangles at each vertex (for plots only)."""
        T = self.mesh.triangles
        vals = np.zeros(self.mesh.n_nodes); cnt = np.zeros(self.mesh.n_nodes)
        for i in range(3):
            # value of the CR function at vertex i of T:  sum_j u_j (1 - 2 lambda_j(v_i)) = u_j(1-2 delta_ji)
            vi = u[self.t2e].sum(axis=1) - 2.0 * u[self.t2e[:, i]]
            np.add.at(vals, T[:, i], vi); np.add.at(cnt, T[:, i], 1.0)
        return vals / cnt

    def l2_error(self, uh, u_exact, relative=True) -> float:
        X, Y = self._qpoints()
        ue = np.asarray(u_exact(X, Y), float)
        uhq = (1.0 - 2.0 * self._QL) @ uh[self.t2e].T
        num = np.sqrt(np.sum(self.area * (self._QW @ (ue - uhq) ** 2)))
        den = np.sqrt(np.sum(self.area * (self._QW @ ue**2))) if relative else 1.0
        return float(num / den)


class _DofMesh:
    """Stand-in for the ``mesh`` attribute used by GLaplaceSolver.solve: only the dof count and the Dirichlet dofs."""
    def __init__(self, n, dirichlet):
        self.n_nodes, self.dirichlet = int(n), np.asarray(dirichlet)


class GLaplaceCRSolver(GLaplaceSolver):
    """Same iteration as :class:`glaplace.solver.GLaplaceSolver` (``solve`` is inherited), with the CR space for u.

    u_D : constant / callable(x, y) (its mean over every boundary edge is imposed) / nodal array of the P1 mesh.
    f_degree is accepted for API compatibility and ignored (midpoint rule)."""

    def __init__(self, mesh: Mesh, G, f=0.0, u_D=0.0, r: float = 1.0, f_degree: int = 2):
        t0 = time.perf_counter()
        if not isinstance(G, Young):
            G = PowerLaw(float(G))
        self.G, self.r, self.true_mesh = G, float(r), mesh
        self.fem = CRSpace(mesh)
        D = self.fem.bdry_edges
        self.mesh = _DofMesh(self.fem.Ne, D)
        self.xb, self.yb = self.fem.bary[:, 0], self.fem.bary[:, 1]
        self.F = self.fem.load_vector(f)
        if np.ndim(u_D) == 1 and not callable(u_D) and len(u_D) == mesh.n_nodes:
            ue = np.asarray(u_D, float)[self.fem.edges[D]].mean(axis=1)
            self.gD = ue
        else:
            self.gD = self.fem.edge_values(u_D, D)
        mask = np.ones(self.fem.Ne, bool); mask[D] = False
        self.free = np.flatnonzero(mask)
        K = self.r * self.fem.K
        self.Kff = K[self.free][:, self.free].tocsc()
        self.Kfd = K[self.free][:, D].tocsr()
        self._lu = spla.splu(self.Kff)
        self.time_setup = time.perf_counter() - t0


# ---------------------------------------------------------------------------------------------
class VarPower(Young):
    """Variable exponent  G(x,t) = t^{p(x)}/p(x),  p(x, y) a callable returning an array (>1), sampled at the barycentres."""

    def __init__(self, pfun):
        self.pfun = pfun

    def _p(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        return np.broadcast_to(np.asarray(self.pfun(x, y), float), x.shape)

    def G(self, t, x=None, y=None):
        p = self._p(x, y)
        return np.asarray(t, float) ** p / p

    def g(self, t, x=None, y=None):
        p = self._p(x, y)
        return np.asarray(t, float) ** (p - 1.0)

    def __repr__(self):
        return "VarPower(p(x))"
