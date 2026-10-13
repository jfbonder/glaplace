"""EXPERIMENT 6 (PDE request E2, Example pde:ex:noS): (S) fails but there is no gap.

p = 3/2 on the unit disk, f = |x|^-beta, u = 0 on the boundary.  Radial solution: |eta| = r^{1-beta}/(2-beta),
|grad u| = |eta|^2 = r^{2(1-beta)}/(2-beta)^2.
  beta = 1.5: u = 4 log(1/r), J(u) = -32 pi/3, grad u not in L^2 ((S) fails), ||eta||^2 = 8 pi.
  beta = 1.4: companion with (S): u = (1 - r^0.2)/0.072.
  beta = 1.6: u = 31.25 (r^-0.2 - 1) (unbounded, in W^{1,3/2}_0, not in H^1).
Meshes: inscribed polygon with a vertex at 0 (v4tools.disk_mesh_resolved(h, circles=())), h = 0.2 ... 0.0125.
Load: 'accurate' = collapsed Gauss rule with Gauss-Jacobi weight s^{1-beta} on the triangles at 0 (exact up to ~1e-10), or
'bary' = barycentre rule f(x_T)|T|/3.  P1 and CR, r = 1, stop residual < 1e-9 and discrete gap < 1e-9 rel, max_iter 50000.
J(u_h) for P1 is evaluated with the accurate load: since u_h (extended by 0) is in W^{1,3/2}_0(B_1), J(u_h) >= J(u) rigorously.

    python experiments_v4/exp06_rough_load.py [--plot]
"""
import argparse, os, sys, time
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import write_rows, read_rows, style, savefig, C, MK
from glaplace.admm import Discretisation, ADMM
from glaplace.elementwise import PowerSum
from glaplace.v4tools import disk_mesh_resolved, mesh_h, SingularQuadrature, barycentre_load


def exact(beta):
    A = 1.0 / (2 - beta) ** 2; a = 2 * (1 - beta)
    ur = lambda r: A * r ** a
    if abs(beta - 1.5) < 1e-14:
        u = lambda r: 4.0 * np.log(1.0 / r)
    else:
        u = lambda r: A * (1.0 - r ** (a + 1)) / (a + 1)
    Jint = lambda r: (ur(r) ** 1.5 / 1.5 - r ** (-beta) * u(r)) * r
    J = 2 * np.pi * (quad(lambda r: ur(r) ** 1.5 / 1.5 * r, 0, 1, limit=200)[0] - quad(lambda r: r ** (1 - beta) * u(r), 0, 1, limit=200)[0])
    eta2 = 2 * np.pi * quad(lambda r: (r ** (1 - beta) / (2 - beta)) ** 2 * r, 0, 1)[0]
    return u, ur, J, eta2


def run(beta, h, scheme, load, r=1.0):
    m = disk_mesh_resolved(h, circles=())
    f = lambda x, y: np.hypot(x, y) ** (-beta)
    d0 = Discretisation(m, scheme)
    Qf = SingularQuadrature(m, gamma=beta, n=8, n_sing=10)
    Facc = Qf.load_P1(m, f) if scheme == "P1" else Qf.load_CR(d0.fem, f)
    F = Facc if load == "accurate" else barycentre_load(d0.fem, f, scheme)
    d = Discretisation(m, scheme, F=F)
    A = ADMM(d, PowerSum(np.ones((1, d.Nt)), np.full((1, d.Nt), 1.5)), r=r)
    o = A.solve(tol_res=1e-9, tol_gap=1e-9, max_iter=50000, check_every=5, record=False)
    u_ex, ur, Jex, eta2 = exact(beta)
    u = o["u"]
    gx, gy = d.fem.gradient(u)
    Jacc = float(np.sum(d.area * np.hypot(gx, gy) ** 1.5 / 1.5) - Facc @ u)
    # true errors
    Qg = SingularQuadrature(m, gamma=1.5 * (2 * beta - 2), n=8, n_sing=10)
    X, Y = Qg.X[..., 0], Qg.X[..., 1]; R = np.hypot(X, Y)
    ux, uy = -ur(R) * X / R, -ur(R) * Y / R
    eg = np.hypot(ux - gx[:, None], uy - gy[:, None])
    gerr = float(Qg.integrate(eg ** 1.5).sum() ** (1 / 1.5))
    gnorm = float(Qg.integrate(np.hypot(ux, uy) ** 1.5).sum() ** (1 / 1.5))
    Qu = SingularQuadrature(m, gamma=max(0.0, 2 * (2 * beta - 3)), n=8, n_sing=10)
    uq = Qu.eval_P1(m, u) if scheme == "P1" else Qu.eval_CR(d.fem, u)
    ue = u_ex(np.hypot(Qu.X[..., 0], Qu.X[..., 1]))
    l2 = float(np.sqrt(Qu.integrate((ue - uq) ** 2).sum() / Qu.integrate(ue ** 2).sum()))
    return dict(beta=beta, scheme=scheme, load=load, h_target=h, h=mesh_h(m), ndof=len(d.free), r=r,
                iterations=o["iterations"], converged=o["converged"], residual=o["residual"], gap=o["gap"],
                J_discrete=o["J"], J_acc=Jacc, J_exact=Jex, J_minus_exact=Jacc - Jex, N2=o["N2"], N2sq=o["N2"] ** 2,
                E2=o["E2"], E2sq_exact=eta2, Einf=o["Einf"], Ginf=o["Ginf"], E0=o["E0"], grad_L15_err=gerr,
                grad_L15_rel=gerr / gnorm, L2_rel=l2, time=o["time_iter"])


def compute():
    rows = []
    for beta in (1.5, 1.4, 1.6):
        for h in (0.2, 0.1, 0.05, 0.025, 0.0125):
            for s in ("P1", "CR"):
                for load in (("accurate", "bary") if beta == 1.5 else ("accurate",)):
                    t0 = time.time()
                    row = run(beta, h, s, load)
                    rows.append(row)
                    print(f"beta={beta} h={h:<7g} {s} {load:8s} it={row['iterations']:6d} conv={row['converged']} "
                          f"J_acc={row['J_acc']:.6f} (ex {row['J_exact']:.6f}) N2^2={row['N2sq']:.3f} E2^2={row['E2']**2:.4f} "
                          f"Einf={row['Einf']:.3f} gradL1.5rel={row['grad_L15_rel']:.3e} ({time.time()-t0:.0f}s)", flush=True)
    write_rows("exp06_rough_load", rows, __file__)


def plot():
    plt = style()
    rows = read_rows("exp06_rough_load")
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6))
    for j, beta in enumerate((1.4, 1.5, 1.6)):
        for s, ls in (("P1", "-"), ("CR", "--")):
            R = sorted([x for x in rows if x["beta"] == beta and x["scheme"] == s and x["load"] == "accurate"], key=lambda x: -x["h"])
            lh = [np.log(1 / x["h"]) for x in R]
            kw = dict(color=C[j], marker=MK[j], mfc="white" if s == "CR" else C[j], label=rf"{s}, $\beta={beta:g}$")
            axs[0].plot(lh, [x["N2sq"] for x in R], ls, **kw)
            axs[1].plot(lh, [x["E2"] ** 2 for x in R], ls, **kw)
            axs[2].plot(lh, [x["iterations"] for x in R], ls, **kw)
    R = sorted([x for x in rows if x["beta"] == 1.5 and x["scheme"] == "P1" and x["load"] == "accurate"], key=lambda x: -x["h"])
    lh = np.array([np.log(1 / x["h"]) for x in R])
    axs[0].plot(lh, 32 * np.pi * lh + (R[-1]["N2sq"] - 32 * np.pi * lh[-1]), "k:", lw=0.8, label=r"slope $32\pi$")
    axs[1].axhline(8 * np.pi, color="k", ls=":", lw=0.8, label=r"$8\pi$")
    axs[0].set_title(r"$\|\nabla_h u_h\|_{L^2}^2$"); axs[1].set_title(r"$\|\eta_h\|_{L^2}^2$"); axs[2].set_title("iterations ($r=1$)")
    for ax in axs:
        ax.set_xlabel(r"$\log(1/h)$")
    axs[0].legend(fontsize=6); axs[1].legend(fontsize=6)
    fig.suptitle(r"$p=3/2$, $f=|x|^{-\beta}$ on $B_1$ ($\beta=1.5$: (S) fails, no gap)", fontsize=9)
    fig.tight_layout(); savefig(fig, "F12_rough_load.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--plot", action="store_true"); a = ap.parse_args()
    plot() if a.plot else compute()
