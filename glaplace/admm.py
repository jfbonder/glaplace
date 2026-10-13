"""v4 ADMM driver for P1 and Crouzeix-Raviart with discrete duality gaps and guaranteed bounds (NEW module).

Same iteration as glaplace.solver.GLaplaceSolver.solve (Glowinski's ALG2 with rho = r, nu_0 = eta_1 = 0 by default):

    1. r (grad u_n, grad phi) = l_h(phi) + (r nu_{n-1} - eta_n, grad phi)        (one sparse LU, reused)
    2. nu_n = argmin_q  G_T(|q|) + r/2 |q - (grad u_n + eta_n/r)|^2              (element-wise, Newton in log t)
    3. eta_{n+1} = eta_n + r (grad u_n - nu_n)

What is new (all in new files; nothing in the existing package is modified):
  * the EQUILIBRATED multiplier  sigma_n = eta_n + r (grad_h u_n - nu_{n-1})  (NA_SECTIONS, Remark "equilibrated"),
    which satisfies (sigma_n, grad_h w)_h = l_h(w) for all w in V_h at every n (up to LU round-off);
  * the discrete dual value  D_h(sigma) = (sigma, grad_h L)_h - l_h(L) - sum_T |T| G_T^*(|sigma_T|),  L = lifting of the data
    (free dofs zero).  Fenchel-Young gives  J_h(v) >= D_h(sigma_n)  for every discrete v with the data, so
    J_h(u_n) - D_h(sigma_n) >= 0 is a discrete primal-dual gap -> stopping test (for P1 AND CR);
  * for CR, f = 0, data = exact edge means, G constant on every triangle (or a minorant G_T):  D_h(sigma_n) equals the
    bound D(sigma_n) of Theorem D and is a GUARANTEED LOWER BOUND for m_W at every iteration (the interior part of
    (sigma, grad_h L) vanishes because L has zero interior midpoint values; the boundary part is sum |E| (sigma.n) psibar_E);
    the max interior normal jump of sigma_n (zero in exact arithmetic) is recorded as a check;
  * for P1 with exact boundary trace and G constant on triangles (or a majorant): J_h(u_n) >= m_H >= m_W at every n;
  * dual residual  r ||nu_n - nu_{n-1}||_h,  warm starts (any nu_0, eta_1 allowed by Theorem C), mesh prolongation;
  * MMD_AT_PLUS_A ordering for the (symmetric) stiffness matrix (3-5x fewer LU nonzeros than the default COLAMD).
"""
from __future__ import annotations
import time
import numpy as np
import scipy.sparse.linalg as spla

from .fem import P1Space
from .cr import CRSpace
from .elementwise import ElementIntegrand, from_young


def _norm(area, x, y):
    return float(np.sqrt(np.sum(area * (x * x + y * y))))


class Discretisation:
    """P1 or CR space on a mesh with Dirichlet data.

    scheme : 'P1' or 'CR'
    gD     : values of the Dirichlet dofs (P1: nodal values at mesh.dirichlet; CR: edge means on fem.bdry_edges)
    F      : load vector (length = number of dofs), default 0
    """

    def __init__(self, mesh, scheme, gD=None, F=None):
        self.mesh, self.scheme = mesh, scheme
        if scheme == "P1":
            self.fem = P1Space(mesh)
            self.ndof = mesh.n_nodes
            self.D = np.asarray(mesh.dirichlet)
        elif scheme == "CR":
            self.fem = CRSpace(mesh)
            self.ndof = self.fem.Ne
            self.D = self.fem.bdry_edges
        else:
            raise ValueError(scheme)
        mask = np.ones(self.ndof, bool); mask[self.D] = False
        self.free = np.flatnonzero(mask)
        self.gD = np.zeros(len(self.D)) if gD is None else np.asarray(gD, float)
        self.F = np.zeros(self.ndof) if F is None else np.asarray(F, float)
        self.area = self.fem.area
        self.xb, self.yb = self.fem.bary[:, 0], self.fem.bary[:, 1]
        self.Nt = self.fem.Nt
        K = self.fem.K
        self.Kff = K[self.free][:, self.free].tocsc()
        self.Kfd = K[self.free][:, self.D].tocsr()
        self.Gxf = self.fem.Gx[:, self.free].T.tocsr()
        self.Gyf = self.fem.Gy[:, self.free].T.tocsr()
        L = np.zeros(self.ndof); L[self.D] = self.gD
        self.lift = L
        self.glx, self.gly = self.fem.gradient(L)            # gradient of the lifting (per element)
        self.lF = float(self.F @ L)

    def lift_full(self, ufree):
        u = self.lift.copy(); u[self.free] = ufree
        return u

    # --- CR edge geometry (for the RT0 normal-jump check and the boundary-flux formula of Theorem D) --------------
    def edge_geometry(self):
        if self.scheme != "CR":
            raise ValueError("CR only")
        if hasattr(self, "_N"):
            return self._N, self._L
        V, T = self.mesh.vertices, self.mesh.triangles
        N = np.zeros((self.Nt, 3, 2)); L = np.zeros((self.Nt, 3))
        for i in range(3):
            a, b, c = V[T[:, (i + 1) % 3]], V[T[:, (i + 2) % 3]], V[T[:, i]]
            t = b - a; L[:, i] = np.hypot(t[:, 0], t[:, 1])
            n = np.column_stack([t[:, 1], -t[:, 0]]) / L[:, i, None]
            sgn = np.sign(np.sum(n * (a - c), axis=1)); N[:, i] = n * sgn[:, None]
        self._N, self._L = N, L
        self._bd = np.zeros(self.ndof, bool); self._bd[self.D] = True
        return N, L

    def normal_jump(self, sx, sy):
        """max over interior edges of |sum_{T ∋ E} |E| sigma_T . n_T|  (0 for an RT0 field)."""
        N, L = self.edge_geometry()
        sn = sx[:, None] * N[:, :, 0] + sy[:, None] * N[:, :, 1]
        flux = np.zeros(self.ndof); np.add.at(flux, self.fem.t2e.ravel(), (sn * L).ravel())
        inner = ~self._bd
        return float(np.max(np.abs(flux[inner]))) if inner.any() else 0.0

    def boundary_flux_term(self, sx, sy, psibar=None):
        """sum_{E boundary} |E| (sigma_T(E).n_E) psibar_E  (Theorem D; psibar = exact edge means, default gD)."""
        N, L = self.edge_geometry()
        sn = sx[:, None] * N[:, :, 0] + sy[:, None] * N[:, :, 1]
        pb = np.zeros(self.ndof); pb[self.D] = self.gD if psibar is None else psibar
        e = self.fem.t2e
        return float(np.sum(np.where(self._bd[e], sn * L * pb[e], 0.0)))


class ADMM:
    """ADMM for min J_h(v) = sum_T |T| G_T(|grad_h v|) - F.v on a Discretisation.

    G : an ElementIntegrand, or a glaplace Young object (sampled at barycentres via elementwise.from_young)
    """

    def __init__(self, disc: Discretisation, G, r=1.0, ordering="MMD_AT_PLUS_A"):
        t0 = time.perf_counter()
        self.d, self.r = disc, float(r)
        self.G = G if isinstance(G, ElementIntegrand) else from_young(G, disc.xb, disc.yb)
        opts = dict(SymmetricMode=True) if ordering != "COLAMD" else {}
        self._lu = spla.splu((self.r * disc.Kff).tocsc(), permc_spec=ordering, options=opts)
        self._rhs0 = disc.F[disc.free] - self.r * (disc.Kfd @ disc.gD)
        self.time_setup = time.perf_counter() - t0

    # ------------------------------------------------------------------------------------------------------------
    def energy_from_grad(self, gx, gy, u):
        return float(np.sum(self.d.area * self.G.G(np.hypot(gx, gy))) - self.d.F @ u)

    def energy(self, u):
        gx, gy = self.d.fem.gradient(u)
        return self.energy_from_grad(gx, gy, u)

    def dual_value(self, sx, sy):
        """D_h(sigma) = (sigma, grad L)_h - l_h(L) - sum |T| G_T^*(|sigma_T|)."""
        d = self.d
        lin = float(np.sum(d.area * (sx * d.glx + sy * d.gly))) - d.lF
        return lin - float(np.sum(d.area * self.G.conj(np.hypot(sx, sy))))

    # ------------------------------------------------------------------------------------------------------------
    def solve(self, tol_res=1e-10, tol_gap=1e-9, tol_dual=None, max_iter=50000, min_iter=1, nu0=None, eta1=None,
              check_every=1, jumps=False, record=True, callback=None, verbose=0, time_limit=None, stop="both"):
        """stop: 'both' (residual < tol_res AND gap < tol_gap*max(1,|J_h|)), 'residual' (residual only; plus dual residual
        if tol_dual is given).  check_every: the gap (needs G^*) is evaluated every k iterations (and at the end).
        jumps: record the max interior normal jump of sigma_n (CR) at every check.
        Returns a dict with u, nu, eta (= eta_{n+1}), sigma (last equilibrated multiplier), histories and summary values."""
        d, r, area = self.d, self.r, self.d.area
        fem = d.fem
        Nt = d.Nt
        nux = np.zeros(Nt) if nu0 is None else np.array(nu0[0], float)
        nuy = np.zeros(Nt) if nu0 is None else np.array(nu0[1], float)
        etx = np.zeros(Nt) if eta1 is None else np.array(eta1[0], float)
        ety = np.zeros(Nt) if eta1 is None else np.array(eta1[1], float)
        tprev = np.hypot(nux, nuy)
        H = {k: [] for k in ("it", "res", "dres", "J", "D", "gap", "jump")} if record else None
        best_D, best_J = -np.inf, np.inf
        conv, it = False, 0
        t0 = time.perf_counter()
        J = Dv = gap = jump = np.nan
        for it in range(1, max_iter + 1):
            rhs = self._rhs0 + d.Gxf @ (area * (r * nux - etx)) + d.Gyf @ (area * (r * nuy - ety))
            u = d.lift_full(self._lu.solve(rhs))
            gx, gy = fem.gradient(u)
            sx, sy = etx + r * (gx - nux), ety + r * (gy - nuy)                   # sigma_n (equilibrated)
            nux_old, nuy_old = nux, nuy
            nux, nuy, tprev = self.G.update_nu(etx / r + gx, ety / r + gy, r, t0=tprev)
            etx = etx + r * (gx - nux); ety = ety + r * (gy - nuy)
            gn = _norm(area, gx, gy)
            res = _norm(area, gx - nux, gy - nuy) / (1.0 + gn)
            dres = r * _norm(area, nux - nux_old, nuy - nuy_old) / (1.0 + _norm(area, etx, ety))
            check = (it % check_every == 0) or it == 1 or (res < tol_res and it >= min_iter)
            if check:
                J = self.energy_from_grad(gx, gy, u)
                Dv = self.dual_value(sx, sy)
                gap = J - Dv
                best_D, best_J = max(best_D, Dv), min(best_J, J)
                jump = d.normal_jump(sx, sy) if (jumps and d.scheme == "CR") else np.nan
                if record:
                    for k, v in zip(("it", "res", "dres", "J", "D", "gap", "jump"), (it, res, dres, J, Dv, gap, jump)):
                        H[k].append(v)
                if callback is not None:
                    callback(it, dict(u=u, nu=(nux, nuy), eta=(etx, ety), sigma=(sx, sy), res=res, J=J, D=Dv))
                if verbose and (it % verbose == 0 or it == 1):
                    print(f"  it {it:6d} res={res:.2e} dres={dres:.2e} J={J:.10f} D={Dv:.10f} gap={gap:.2e}", flush=True)
                ok = res < tol_res and it >= min_iter
                if stop == "both":
                    ok = ok and gap < tol_gap * max(1.0, abs(J))
                if tol_dual is not None:
                    ok = ok and dres < tol_dual
                if ok:
                    conv = True
                    break
            if time_limit is not None and time.perf_counter() - t0 > time_limit:
                break
        # final values at the last iterate
        J = self.energy_from_grad(gx, gy, u); Dv = self.dual_value(sx, sy)
        best_D, best_J = max(best_D, Dv), min(best_J, J)
        out = dict(u=u, nu=(nux, nuy), eta=(etx, ety), sigma=(sx, sy), iterations=it, converged=conv,
                   residual=res, dual_residual=dres, J=J, D=Dv, gap=J - Dv, best_D=best_D, best_J=best_J,
                   time_iter=time.perf_counter() - t0, time_setup=self.time_setup, r=r,
                   max_jump=(max(np.nanmax(H["jump"]) if record and H["jump"] and np.isfinite(H["jump"]).any() else -1,
                                 d.normal_jump(sx, sy)) if d.scheme == "CR" else np.nan))
        if record:
            out["history"] = {k: np.array(v) for k, v in H.items()}
        out.update(self.diagnostics(u, (etx, ety)))
        return out

    def diagnostics(self, u, eta):
        """Quantities of Thm pde:thm:diag: N2 = ||grad_h u||, E2 = ||eta||, Einf = max|eta_T|, E0 = E2^2/(2r) + r N2^2/2."""
        gx, gy = self.d.fem.gradient(u)
        N2 = _norm(self.d.area, gx, gy)
        E2 = _norm(self.d.area, *eta)
        Einf = float(np.max(np.hypot(*eta)))
        Ginf = float(np.max(np.hypot(gx, gy)))
        return dict(N2=N2, E2=E2, Einf=Einf, Ginf=Ginf, E0=E2 ** 2 / (2 * self.r) + self.r * N2 ** 2 / 2)


# ----------------------------------------------------------------------------------------------------------------
def prolong(coarse_mesh, fine_mesh, fields):
    """Piecewise-constant transfer coarse -> fine: each fine triangle takes the value of the coarse triangle containing
    its barycentre (exact for nested meshes, e.g. crisscross_mesh(n) -> crisscross_mesh(2n)).  fields: list of arrays."""
    import matplotlib.tri as mtri
    tri = mtri.Triangulation(coarse_mesh.vertices[:, 0], coarse_mesh.vertices[:, 1], coarse_mesh.triangles)
    fb = fine_mesh.vertices[fine_mesh.triangles].mean(axis=1)
    idx = tri.get_trifinder()(fb[:, 0], fb[:, 1])
    if np.any(idx < 0):
        raise ValueError("fine barycentre outside the coarse mesh")
    return [np.asarray(f)[idx] for f in fields]


def source_hash(paths=None):
    """sha256 (first 12 hex digits) of the package source files: the 'code version' recorded in every result file."""
    import hashlib, glob, os
    here = os.path.dirname(__file__)
    paths = sorted(glob.glob(os.path.join(here, "*.py"))) if paths is None else paths
    h = hashlib.sha256()
    for p in paths:
        h.update(os.path.basename(p).encode()); h.update(open(p, "rb").read())
    return h.hexdigest()[:12]
