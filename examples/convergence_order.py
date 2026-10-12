"""Order of convergence in h for radial g-torsion problems with (semi-)explicit solutions.

    -div( g(x,|grad u|) grad u/|grad u| ) = c  in B_1,   u = 0 on dB_1,     g(|u'|) = c r/2.

For each Young function G the solution u(r) = int_r^1 g^{-1}(c s/2) ds is computed to ~1e-10 (glaplace.exact.radial_torsion).
We report TRUE errors (quadrature, not nodal differences): relative L2 error of u, relative L2 error of grad u, nodal sup error,
and fit  error ~ h^t  on the 3 finest meshes.

Run:  python examples/convergence_order.py [--cases pq24 ...] [--h 0.2 0.1 0.05 0.025 0.0125] [--plot] [--jobs 2]
"""
import argparse, csv, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, __file__.rsplit("/examples", 1)[0])
import glaplace as gl

ring = lambda x, y: np.where(np.hypot(x, y) > 0.5, 1.0, 0.0)      # a = 0 for r < 0.5 (two-phase, radial)

CASES = {   # name: (G, c, r_aug, description)
    "p1.5":      (lambda: gl.PowerLaw(1.5), 1.0, 1.0, "t^1.5/1.5"),
    "p2":        (lambda: gl.PowerLaw(2.0), 1.0, 1.0, "t^2/2 (Laplacian)"),
    "p4":        (lambda: gl.PowerLaw(4.0), 1.0, 1.0, "t^4/4"),
    "pq24":      (lambda: gl.pq(2, 4), 1.0, 1.0, "t^2/2 + t^4/4"),
    "pq1.5_3":   (lambda: gl.pq(1.5, 3), 1.0, 1.0, "t^1.5/1.5 + t^3/3"),
    "plog2_1":   (lambda: gl.PLog(2, 1), 1.0, 1.0, "t^2 log(1+t)"),
    "plog3_-1":  (lambda: gl.PLog(3, -1), 1.0, 1.0, "t^3 / log(1+t)"),
    "exp":       (lambda: gl.Exp(), 1.0, 1.0, "e^t - 1 - t"),
    "twophase":  (lambda: gl.pq(2, 4, ring), 1.0, 1.0, "t^2/2 + a(r) t^4/4, a = 1_{r>1/2}"),
    "expdead":   (lambda: gl.ExpLiteral(), 8.0, 30.0, "e^t - 1, c = 8 (dead core r < 1/4)"),
}


def run_one(args):
    name, h, tol, max_iter = args
    mk, c, r, _ = CASES[name]
    G = mk()
    u, gu, _ = gl.exact.radial_torsion(G, c)
    mesh = gl.disk_mesh(h)
    t0 = time.perf_counter()
    res = gl.GLaplaceSolver(mesh, G, c, 0.0, r=r).solve(tol=tol, max_iter=max_iter, stop="residual")
    F = gl.P1Space(mesh)
    ue = u(*mesh.vertices.T)
    return dict(case=name, h=h, nodes=mesh.n_nodes, iterations=res.iterations, converged=res.converged,
                residual=float(res.history["residual"][-1]), L2_rel=F.l2_error(res.u, u),
                gradL2_rel=F.grad_error_Lp(res.u, gu, 2.0), sup_nodes=float(np.abs(res.u - ue).max()),
                time=time.perf_counter() - t0)


def fit(rows, key, k=3):
    rows = sorted(rows, key=lambda r: -r["h"])[-k:]
    return np.polyfit(np.log([r["h"] for r in rows]), np.log([r[key] for r in rows]), 1)[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="+", default=list(CASES))
    ap.add_argument("--h", type=float, nargs="+", default=[0.2, 0.1, 0.05, 0.025, 0.0125])
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--max-iter", type=int, default=3000)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--csv", default="results/convergence_order.csv")
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    jobs = [(n, h, a.tol, a.max_iter) for n in a.cases for h in a.h]
    if a.jobs > 1:
        with ProcessPoolExecutor(a.jobs) as ex:
            out = list(ex.map(run_one, jobs))
    else:
        out = [run_one(j) for j in jobs]
    keys = list(out[0])
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, keys); w.writeheader(); w.writerows(out)
    print(f"{'case':10s} {'G':34s} {'t(L2)':>6} {'t(grad)':>8} {'t(sup)':>7}  converged(finest)")
    summary = {}
    for n in a.cases:
        rows = [r for r in out if r["case"] == n]
        summary[n] = (fit(rows, "L2_rel"), fit(rows, "gradL2_rel"), fit(rows, "sup_nodes"))
        fin = sorted(rows, key=lambda r: r["h"])[0]
        print(f"{n:10s} {CASES[n][3]:34s} {summary[n][0]:6.2f} {summary[n][1]:8.2f} {summary[n][2]:7.2f}  "
              f"{fin['converged']} ({fin['iterations']} it)")
    if a.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for n in a.cases:
            rows = sorted([r for r in out if r["case"] == n], key=lambda r: -r["h"])
            h = [r["h"] for r in rows]
            ax[0].loglog(h, [r["L2_rel"] for r in rows], "o-", label=n)
            ax[1].loglog(h, [r["gradL2_rel"] for r in rows], "o-", label=n)
        for k, (e, ttl) in enumerate([(2, "relative $L^2$ error of $u$"), (1, "relative $L^2$ error of $\\nabla u$")]):
            hh = np.array([a.h[0], a.h[-1]]); ax[k].loglog(hh, 0.5 * hh ** e * ax[k].get_ylim()[1] / hh[0] ** e * 0.3, "k--", lw=0.8, label=f"$h^{e}$")
            ax[k].set_xlabel("h"); ax[k].set_title(ttl)
        ax[0].legend(fontsize=7)
        plt.tight_layout(); plt.savefig("docs/convergence_order.png", dpi=110)
