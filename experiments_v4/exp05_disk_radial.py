"""EXPERIMENT 5 (PDE requests E1 and E4(ii) + section 0): radial torsion on the unit disk with exact solutions,
variable exponent with p^- < 2 (new Theorem B / Cor. B) and radial two-phase (no gap, (S) holds).

-div(g(x,|grad u|) grad u/|grad u|) = 1 in B_1, u = 0 on the boundary, exact flux eta = -x/2, |grad u| = g^{-1}(|x|/2).
Cases (G sampled per triangle):
  E1a   t^{p(r)}/p(r), p = 1.4 + 0.4 r        (p+ = 1.8 < 2)          p sampled at the barycentre radius
  E1b   t^{p(r)}/p(r), p = 1.5 + r            (crosses 2)
  E1c   p = 1.5 (r < 1/2), 3 (r > 1/2)        (discontinuous)          triangle inside the polygon r=1/2 <=> all |vertices| <= 1/2
  p1.5, p3                                     (reference constant exponents)
  E4ii_24   t^2/2 + 1_{r>1/2} t^4/4            (old two-phase test)     same interface classification as E1c
  E4ii_153  t^1.5/1.5 + 1_{r>1/2} t^3/3
Meshes: Triangle meshes of the inscribed polygon with vertices on r = 1/2 and at 0 (glaplace.v4tools.disk_mesh_resolved),
target h = 0.2, 0.1, 0.05, 0.025, 0.0125.  The discrete interface is the inscribed polygon (O(h^2) from the circle) and
Omega_h is the inscribed polygon (u_h = 0 on its boundary).  P1 and CR; r = 1; nu_0 = eta_1 = 0;
stop: residual < 1e-9 and discrete duality gap < 1e-9 rel (f != 0: discrete gap, not a bound for the continuous minimum);
max_iter 50000.  Load: P1 int f phi_i exact (f = 1); CR midpoint rule exact.
Errors with a collapsed Gauss rule (64 points per triangle): relative ||u - u_h||_L2, ||grad(u - u_h)||_{L^s}, s = 1.2, 1.5,
1.9, 2 (relative), modular int |grad(u - u_h)|^{p(x)} (p at the quadrature point, true p for E1c).
Iteration histories on the finest mesh (E1a-c, P1): int |nu_n - grad u_h|^{p_T}, ||nu_n - grad u_h||_{L^1.5},
||eta_n - eta_h||_{L2}, with (u_h, eta_h) from a run to residual 1e-12.

    python experiments_v4/exp05_disk_radial.py [--hs ...] [--plot]
"""
import argparse, os, sys, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import write_rows, write_json, read_rows, style, savefig, C, MK, RESULTS
from glaplace.admm import Discretisation, ADMM
from glaplace.elementwise import PowerSum
from glaplace.v4tools import disk_mesh_resolved, mesh_h, SingularQuadrature
from glaplace.young import Young
from glaplace import exact


class RadialPow(Young):
    """sum_k w_k(r) t^{p_k(r)}/p_k(r) with r = hypot(x, y) (only for the exact radial solution)."""
    def __init__(self, terms):
        self.terms = terms                                   # list of (w(r), p(r))
    def G(self, t, x=None, y=None):
        r = np.hypot(x, y); t = np.asarray(t, float)
        return sum(w(r) * t ** p(r) / p(r) for w, p in self.terms)
    def g(self, t, x=None, y=None):
        r = np.hypot(x, y); t = np.asarray(t, float)
        return sum(w(r) * t ** (p(r) - 1.0) for w, p in self.terms)


one = lambda r: np.ones_like(np.asarray(r, float))
cst = lambda c: (lambda r: np.full_like(np.asarray(r, float), c))
ind = lambda r: (np.asarray(r, float) > 0.5).astype(float)
CASES = {
    "E1a": ([(one, lambda r: 1.4 + 0.4 * np.asarray(r))], "smooth"),
    "E1b": ([(one, lambda r: 1.5 + np.asarray(r))], "smooth"),
    "E1c": ([(one, lambda r: np.where(np.asarray(r) > 0.5, 3.0, 1.5))], "interface"),
    "p1.5": ([(one, cst(1.5))], "smooth"),
    "p3": ([(one, cst(3.0))], "smooth"),
    "E4ii_24": ([(one, cst(2.0)), (ind, cst(4.0))], "interface"),
    "E4ii_153": ([(one, cst(1.5)), (ind, cst(3.0))], "interface"),
}


def element_integrand(case, d):
    terms, kind = CASES[case]
    V, T = d.mesh.vertices, d.mesh.triangles
    if kind == "smooth":
        rT = np.hypot(d.xb, d.yb)
    else:   # outer phase iff some vertex outside the polygon r = 1/2
        rmax = np.max(np.hypot(V[T][..., 0], V[T][..., 1]), axis=1)
        rT = np.where(rmax > 0.5 + 1e-12, 0.75, 0.25)
    c = np.vstack([w(rT) for w, p in terms]); p = np.vstack([p(rT) for w, p in terms])
    return PowerSum(c, p), rT


def errors(case, d, u, Q, ue, gue, pfun):
    X, Y = Q.X[..., 0], Q.X[..., 1]
    if d.scheme == "P1":
        uq = Q.eval_P1(d.mesh, u)
    else:
        uq = Q.eval_CR(d.fem, u)
    ex = ue(X, Y)
    out = dict(L2_rel=float(np.sqrt(Q.integrate((ex - uq) ** 2).sum() / Q.integrate(ex ** 2).sum())))
    gx, gy = d.fem.gradient(u)
    ux, uy = gue(X, Y)
    e = np.hypot(ux - gx[:, None], uy - gy[:, None]); n = np.hypot(ux, uy)
    for s in (1.2, 1.5, 1.9, 2.0):
        out[f"grad_L{s:g}_rel"] = float((Q.integrate(e ** s).sum() / Q.integrate(n ** s).sum()) ** (1 / s))
    pq = pfun(np.hypot(X, Y))
    out["modular"] = float(Q.integrate(e ** pq).sum())
    return out


def run_case(case, h, scheme, r=1.0, tol=1e-9, history=False):
    terms, kind = CASES[case]
    m = disk_mesh_resolved(h)
    d = Discretisation(m, scheme)
    G, rT = element_integrand(case, d)
    if scheme == "P1":
        F = d.fem.load_vector(1.0, degree=1)          # exact for f = 1
    else:
        F = d.fem.load_vector(1.0)                    # midpoint rule, exact for f = 1
    d = Discretisation(m, scheme, F=F)
    A = ADMM(d, G, r=r)
    o = A.solve(tol_res=tol, tol_gap=1e-9, max_iter=50000, check_every=5, record=False)
    Gex = RadialPow(terms)
    ue, gue, _ = exact.radial_torsion(Gex, 1.0)
    Q = SingularQuadrature(m, gamma=0.0, n=8, n_sing=8)
    pfun = terms[0][1] if case.startswith("E1") or case.startswith("p") else (lambda rr: np.full_like(rr, 2.0))
    row = dict(case=case, scheme=scheme, h_target=h, h=mesh_h(m), nodes=m.n_nodes, ndof=len(d.free), r=r,
               iterations=o["iterations"], converged=o["converged"], residual=o["residual"], gap=o["gap"], J=o["J"],
               N2=o["N2"], E2=o["E2"], Einf=o["Einf"], Ginf=o["Ginf"], E0=o["E0"], time=o["time_iter"])
    row.update(errors(case, d, o["u"], Q, ue, gue, pfun))
    # flux error ||eta_h - eta||_L2 with eta = -x/2 (exact flux), quadrature
    X, Y = Q.X[..., 0], Q.X[..., 1]
    ex, ey = o["eta"]
    row["eta_L2_err"] = float(np.sqrt(Q.integrate((ex[:, None] + X / 2) ** 2 + (ey[:, None] + Y / 2) ** 2).sum()))
    hist = None
    if history:
        ref = A.solve(tol_res=1e-12, tol_gap=1e-12, max_iter=50000, check_every=10, record=False)
        ux, uy = d.fem.gradient(ref["u"]); hx, hy = ref["eta"]
        pT = G.p[0]
        rec = dict(mod=[], L15=[], eta=[])
        def cb(it, st):
            nx, ny = st["nu"]; e = np.hypot(nx - ux, ny - uy)
            rec["mod"].append(float(np.sum(d.area * e ** pT)))
            rec["L15"].append(float(np.sum(d.area * e ** 1.5) ** (1 / 1.5)))
            # eta_n (multiplier used at step n) = st['eta'] is eta_{n+1}; record that
            ex_, ey_ = st["eta"]; rec["eta"].append(float(np.sqrt(np.sum(d.area * ((ex_ - hx) ** 2 + (ey_ - hy) ** 2)))))
        A.solve(tol_res=1e-12, tol_gap=1e-12, max_iter=ref["iterations"], check_every=1, record=False, callback=cb)
        hist = dict(case=case, h=h, ref_iterations=ref["iterations"], ref_residual=ref["residual"], **rec)
    return row, hist


def compute(hs):
    rows, hists = [], []
    for case in CASES:
        for h in hs:
            for s in ("P1", "CR"):
                want = (h == min(hs)) and case in ("E1a", "E1b", "E1c") and s == "P1"
                t0 = time.time()
                row, hist = run_case(case, h, s, history=want)
                rows.append(row)
                if hist:
                    hists.append(hist)
                print(f"{case:9s} h={h:<7g} {s} ndof={row['ndof']:7d} it={row['iterations']:6d} conv={row['converged']} "
                      f"L2={row['L2_rel']:.3e} gL2={row['grad_L2_rel']:.3e} N2={row['N2']:.4f} E2={row['E2']:.4f} "
                      f"Einf={row['Einf']:.4f} ({time.time() - t0:.0f}s)", flush=True)
    write_rows("exp05_disk_radial", rows, __file__, setting="unit disk, f=1, r=1, tol_res=1e-9")
    write_json("exp05_histories", hists, __file__)


def fit(rows, key, k=3):
    rows = sorted(rows, key=lambda x: -x["h"])[-k:]
    return float(np.polyfit(np.log([x["h"] for x in rows]), np.log([x[key] for x in rows]), 1)[0])


def plot():
    plt = style()
    rows = read_rows("exp05_disk_radial")
    cases = list(CASES)
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6))
    for j, case in enumerate(cases):
        R = sorted([x for x in rows if x["case"] == case and x["scheme"] == "P1"], key=lambda x: -x["h"])
        hh = [x["h"] for x in R]
        axs[0].loglog(hh, [x["L2_rel"] for x in R], color=C[j % 7], marker=MK[j % 7], label=case)
        axs[1].loglog(hh, [x["grad_L2_rel"] for x in R], color=C[j % 7], marker=MK[j % 7])
        axs[2].semilogx(hh, [x["iterations"] for x in R], color=C[j % 7], marker=MK[j % 7])
    for ax, t in zip(axs, (r"rel. $L^2$ error of $u_h$ (P1)", r"rel. $L^2$ error of $\nabla u_h$ (P1)", "iterations ($r=1$)")):
        ax.set_xlabel("mesh size $h$ (max edge)"); ax.set_title(t); ax.invert_xaxis()
    axs[0].legend(fontsize=6.5)
    fig.tight_layout(); savefig(fig, "F9_radial_errors.pdf")
    # diagnostics N2, E2, Einf vs log(1/h)
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6), sharex=True)
    for j, case in enumerate(cases):
        for s, ls in (("P1", "-"), ("CR", "--")):
            R = sorted([x for x in rows if x["case"] == case and x["scheme"] == s], key=lambda x: -x["h"])
            lh = [np.log(1 / x["h"]) for x in R]
            for ax, k in zip(axs, ("N2", "E2", "Einf")):
                ax.plot(lh, [x[k] for x in R], ls, color=C[j % 7], marker=MK[j % 7], mfc="white" if s == "CR" else C[j % 7],
                        label=f"{case}" if s == "P1" else None)
    for ax, t in zip(axs, (r"$\|\nabla_h u_h\|_{L^2}$", r"$\|\eta_h\|_{L^2}$", r"$\max_T|\eta_h|$")):
        ax.set_title(t); ax.set_xlabel(r"$\log(1/h)$")
    axs[0].legend(fontsize=6.5); fig.suptitle("radial torsion: discrete saddle point size (solid P1, dashed CR)", fontsize=9)
    fig.tight_layout(); savefig(fig, "F10_radial_diagnostics.pdf")
    H = json.load(open(os.path.join(RESULTS, "exp05_histories.json")))["data"]
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.5))
    for j, h in enumerate(H):
        for ax, k in zip(axs, ("mod", "L15", "eta")):
            v = np.array(h[k]); ax.semilogy(np.arange(1, len(v) + 1), np.maximum(v, 1e-18), color=C[j], ls=["-", "--", "-."][j], label=h["case"])
    for ax, t in zip(axs, (r"$\sum_T|T|\,|\nu_n-\nabla u_h|^{p_T}$", r"$\|\nu_n-\nabla u_h\|_{L^{1.5}}$", r"$\|\eta_{n+1}-\eta_h\|_{L^2}$")):
        ax.set_title(t); ax.set_xlabel("iteration $n$")
    axs[0].legend(); fig.tight_layout(); savefig(fig, "F11_radial_history.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--hs", type=float, nargs="+", default=[0.2, 0.1, 0.05, 0.025, 0.0125])
    ap.add_argument("--plot", action="store_true"); a = ap.parse_args()
    plot() if a.plot else compute(a.hs)
