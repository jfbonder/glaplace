"""EXPERIMENT 8: re-run, with the v4 code (glaplace.admm, P1), the experiments of the v3 article kept in the new story,
so that all numbers come from one code version.

(a) radial torsion orders (v3 Table 'orders', Fig. 'orders'): 10 Young functions on the DistMesh disk meshes
    gl.disk_mesh(h), h = 0.2 ... 0.0125, f = c, u = 0; r = 1 (expdead: c = 8, r = 30); tol: residual < 1e-9 AND discrete
    duality gap < 1e-9 rel (new; v3 used the residual only), max_iter 20000 (v3: 3000).  True errors (degree-4 Dunavant);
    orders fitted over the 3 finest meshes.
(b) p -> 1: t^p/p, p = 1.2, 1.1, h = 0.2 ... 0.025, max_iter 50000.
(c) manufactured x-dependent G = t^2/2 + (1+x1) t^4/4, u = exp(x1+x2) on (0,1)^2, n = 10 ... 80, r = 100, tol 1e-10
    (P1 3-point load).  r-study at n = 40: r = 1, 10, 30, 100, 300.
(d) dead core G = e^t - 1, c = 8, h = 0.05: r = 1, 3, 10, 30, 100, tol 1e-8 (v3 'r' section), error and flatness of the core.
(e) two-phase Hoelder t^2/2 + dist(x,S) t^3/3 on the disk, f = 4, r = 3, tol 1e-9, h = 0.2 ... 0.025 (+ indicator).
(f) Lyapunov check (v3 Fig. 'lyap') on the disk h = 0.1 for 5 integrands; E_n with the proxy saddle point of a long run.

    python experiments_v4/exp08_old_experiments.py [--parts a b c d e f] [--plot]
"""
import argparse, os, sys, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import write_rows, write_json, read_rows, style, savefig, C, MK, RESULTS
import glaplace as gl
from glaplace.admm import Discretisation, ADMM
from glaplace.v4tools import lyapunov

ring = lambda x, y: np.where(np.hypot(x, y) > 0.5, 1.0, 0.0)
CASES = {
    "p1.5": (lambda: gl.PowerLaw(1.5), 1.0, 1.0, "t^1.5/1.5"),
    "p2": (lambda: gl.PowerLaw(2.0), 1.0, 1.0, "t^2/2"),
    "p4": (lambda: gl.PowerLaw(4.0), 1.0, 1.0, "t^4/4"),
    "pq24": (lambda: gl.pq(2, 4), 1.0, 1.0, "t^2/2+t^4/4"),
    "pq1.5_3": (lambda: gl.pq(1.5, 3), 1.0, 1.0, "t^1.5/1.5+t^3/3"),
    "plog2_1": (lambda: gl.PLog(2, 1), 1.0, 1.0, "t^2 log(1+t)"),
    "plog3_-1": (lambda: gl.PLog(3, -1), 1.0, 1.0, "t^3/log(1+t)"),
    "exp": (lambda: gl.Exp(), 1.0, 1.0, "e^t-1-t"),
    "twophase": (lambda: gl.pq(2, 4, ring), 1.0, 1.0, "t^2/2+1_{|x|>1/2} t^4/4"),
    "expdead": (lambda: gl.ExpLiteral(), 8.0, 30.0, "e^t-1, c=8"),
    "p1.2": (lambda: gl.PowerLaw(1.2), 1.0, 1.0, "t^1.2/1.2"),
    "p1.1": (lambda: gl.PowerLaw(1.1), 1.0, 1.0, "t^1.1/1.1"),
}
HS = [0.2, 0.1, 0.05, 0.025, 0.0125]


def p1_disc(mesh, f, uD=0.0, degree=1):
    fem = gl.P1Space(mesh)
    gD = fem.node_values(uD, mesh.dirichlet)
    return Discretisation(mesh, "P1", gD=gD, F=fem.load_vector(f, degree=degree)), fem


def part_a(hs=HS, names=None, max_iter=20000):
    rows = []
    names = names or [k for k in CASES if k not in ("p1.2", "p1.1")]
    for name in names:
        mk, c, r, desc = CASES[name]
        G = mk(); u, gu, _ = gl.exact.radial_torsion(G, c)
        for h in hs:
            mesh = gl.disk_mesh(h)
            d, fem = p1_disc(mesh, c)
            A = ADMM(d, G, r=r)
            o = A.solve(tol_res=1e-9, tol_gap=1e-9, max_iter=max_iter, check_every=5, record=False)
            row = dict(case=name, G=desc, c=c, r=r, h=h, nodes=mesh.n_nodes, iterations=o["iterations"], converged=o["converged"],
                       residual=o["residual"], gap=o["gap"], J=o["J"], L2_rel=fem.l2_error(o["u"], u),
                       gradL2_rel=fem.grad_error_Lp(o["u"], gu, 2.0), sup_nodes=float(np.abs(o["u"] - u(*mesh.vertices.T)).max()),
                       time=o["time_iter"])
            if name == "expdead":
                gx, gy = fem.gradient(o["u"]); core = np.hypot(*fem.bary.T) < 0.2
                row["core_max_grad"] = float(np.hypot(gx, gy)[core].max())
            rows.append(row)
            print(f"{name:9s} h={h:<7g} it={row['iterations']:6d} conv={row['converged']} L2={row['L2_rel']:.3e} "
                  f"g={row['gradL2_rel']:.3e} sup={row['sup_nodes']:.3e} ({row['time']:.0f}s)", flush=True)
    return rows


def fit(rows, key, k=3):
    rows = sorted(rows, key=lambda x: -x["h"])[-k:]
    return float(np.polyfit(np.log([x["h"] for x in rows]), np.log([x[key] for x in rows]), 1)[0])


def part_c():
    f, u, gu = gl.exact.manufactured("t**2/2 + (1+x)*t**4/4", "exp(x+y)")
    G = gl.pq(2, 4, lambda x, y: 1 + x)
    rows = []
    for n in (10, 20, 40, 80):
        m = gl.rectangle_mesh(n, n); d, fem = p1_disc(m, f, u, degree=2)
        o = ADMM(d, G, r=100.0).solve(tol_res=1e-10, tol_gap=1e-10, max_iter=20000, record=False)
        rows.append(dict(part="manufactured", n=n, r=100.0, iterations=o["iterations"], converged=o["converged"],
                         L2_rel=fem.l2_error(o["u"], u), gradL2_rel=fem.grad_error_Lp(o["u"], gu, 2.0), time=o["time_iter"]))
        print(rows[-1], flush=True)
    m = gl.rectangle_mesh(40, 40); d, fem = p1_disc(m, f, u, degree=2)
    for r in (1.0, 10.0, 30.0, 100.0, 300.0):
        o = ADMM(d, G, r=r).solve(tol_res=1e-10, tol_gap=1e-10, max_iter=20000, record=False)
        rows.append(dict(part="manufactured_r", n=40, r=r, iterations=o["iterations"], converged=o["converged"],
                         residual=o["residual"], dual_residual=o["dual_residual"], gap=o["gap"],
                         L2_rel=fem.l2_error(o["u"], u), gradL2_rel=fem.grad_error_Lp(o["u"], gu, 2.0), time=o["time_iter"]))
        print(rows[-1], flush=True)
    return rows


def part_d():
    G = gl.ExpLiteral(); u, gu, _ = gl.exact.radial_torsion(G, 8.0)
    mesh = gl.disk_mesh(0.05); d, fem = p1_disc(mesh, 8.0)
    rows = []
    for r in (1.0, 3.0, 10.0, 30.0, 100.0):
        for stop, tg in (("residual", None), ("both", 1e-8)):
            A = ADMM(d, G, r=r)
            o = A.solve(tol_res=1e-8, tol_gap=1e-8, max_iter=20000, record=False, stop=stop)
            gx, gy = fem.gradient(o["u"]); core = np.hypot(*fem.bary.T) < 2 / 8 - 0.05
            rows.append(dict(part="deadcore_r", h=0.05, r=r, stop=stop, iterations=o["iterations"], converged=o["converged"],
                             residual=o["residual"], dual_residual=o["dual_residual"], gap=o["gap"],
                             L2_rel=fem.l2_error(o["u"], u), core_max_grad=float(np.hypot(gx, gy)[core].max()), time=o["time_iter"]))
            print(rows[-1], flush=True)
            if r == 30.0 and stop == "both":
                np.save(os.path.join(RESULTS, "exp08_deadcore_u.npy"), np.column_stack([np.hypot(*mesh.vertices.T), o["u"]]))
    return rows


def part_e():
    dist = lambda x, y: np.hypot(np.minimum(x, 0.0), np.minimum(y, 0.0))
    aI = lambda x, y: np.where((x >= 0) & (y >= 0), 0.0, 1.0)
    rows = []
    for name, aa in (("Holder", dist), ("indicator", aI)):
        for h in (0.2, 0.1, 0.05, 0.025):
            mesh = gl.disk_mesh(h); d, fem = p1_disc(mesh, 4.0)
            o = ADMM(d, gl.pq(2, 3, aa), r=3.0).solve(tol_res=1e-9, tol_gap=1e-9, max_iter=20000, record=False)
            gx, gy = fem.gradient(o["u"]); gn = np.hypot(gx, gy)
            rows.append(dict(part="twophase_holder", a=name, h=h, nodes=mesh.n_nodes, iterations=o["iterations"],
                             converged=o["converged"], J=o["J"], max_u=float(o["u"].max()), max_grad=float(gn.max()),
                             gradL2=o["N2"], E2=o["E2"], Einf=o["Einf"], time=o["time_iter"]))
            print(rows[-1], flush=True)
    return rows


def part_f():
    out = {}
    mesh = gl.disk_mesh(0.1)
    for name, Gf, c, r, N in (("pq24", gl.pq(2, 4), 1.0, 1.0, 1500), ("p1.5", gl.PowerLaw(1.5), 1.0, 1.0, 2500),
                              ("p4", gl.PowerLaw(4.0), 1.0, 1.0, 2500), ("exp", gl.Exp(), 1.0, 1.0, 1500),
                              ("expdead", gl.ExpLiteral(), 8.0, 30.0, 2500)):
        d, fem = p1_disc(mesh, c)
        A = ADMM(d, Gf, r=r)
        ref = A.solve(tol_res=0, tol_gap=0, max_iter=N + 1, record=False)
        # proxy saddle point: nu* = grad u_L, eta* = sigma_L (equilibrated)
        gx, gy = fem.gradient(ref["u"]); nus = (gx, gy); etas = ref["sigma"]
        E = []
        A.solve(tol_res=0, tol_gap=0, max_iter=N, record=False,
                callback=lambda it, st: E.append(lyapunov(d.area, r, st["nu"], st["eta"], nus, etas)))
        E = np.array(E)
        out[name] = dict(E=E.tolist(), r=r, monotone_until_1e_20=bool(np.all(np.diff(E[E / E[0] > 1e-20]) < 0)))
        print(name, out[name]["monotone_until_1e_20"], E[-1] / E[0], flush=True)
    return out


def plot():
    plt = style()
    rows = read_rows("exp08_orders")
    names = [k for k in CASES if k not in ("p1.2", "p1.1")]
    hl = {"p4", "plog2_1", "twophase", "expdead"}
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6))
    j = 0
    for name in names:
        R = sorted([x for x in rows if x["case"] == name], key=lambda x: -x["h"])
        hh = [x["h"] for x in R]
        if name in hl:
            kw = dict(color=C[j], marker=MK[j], label=name, zorder=3); j += 1
        else:
            kw = dict(color="#999999", lw=0.8, marker=".", zorder=1)
        for ax, k in zip(axs, ("L2_rel", "gradL2_rel", "sup_nodes")):
            ax.loglog(hh, [x[k] for x in R], **kw)
    hh = np.array([0.2, 0.0125])
    for ax, e, t in zip(axs, (2, 1, 2), (r"rel. $L^2$ error of $u$", r"rel. $L^2$ error of $\nabla u$", "nodal max error")):
        ax.loglog(hh, (hh / 0.2) ** e * (2e-2 if e == 2 else 2e-1), "k:", lw=0.8, label=f"$h^{e}$")
        ax.set_xlabel("target mesh size $h$"); ax.set_title(t)
    axs[0].legend(fontsize=6.5)
    fig.tight_layout(); savefig(fig, "F15_convergence_order.pdf")
    # dead core profile
    D = np.load(os.path.join(RESULTS, "exp08_deadcore_u.npy"))
    u, _, _ = gl.exact.radial_torsion(gl.ExpLiteral(), 8.0)
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    rr = np.linspace(0, 1, 400); ax.plot(rr, u(rr, 0 * rr), color=C[6], lw=1.2, label="exact")
    ax.plot(D[:, 0], D[:, 1], ".", color=C[0], ms=1.5, label="$u_h$ (nodes)")
    ax.axvline(0.25, color="gray", ls=":", lw=0.8); ax.set_xlabel(r"$|x|$"); ax.set_ylabel("$u$"); ax.legend()
    ax.set_title(r"$G=e^t-1$, $c=8$: dead core $|x|\leq1/4$")
    fig.tight_layout(); savefig(fig, "F16_dead_core.pdf")
    L = json.load(open(os.path.join(RESULTS, "exp08_lyapunov.json")))["data"]
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    for j, (k, v) in enumerate(L.items()):
        E = np.array(v["E"]); m = E / E[0] > 1e-24
        ax.semilogy(np.arange(1, len(E) + 1)[m], (E / E[0])[m], color=C[j], ls=["-", "--", "-.", ":", "-"][j], label=k)
    ax.set_xlabel("iteration $n$"); ax.set_ylabel(r"$E_n/E_1$"); ax.legend(fontsize=6.5)
    fig.tight_layout(); savefig(fig, "F17_lyapunov_disk.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["a", "b", "c", "d", "e", "f"])
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    if a.plot:
        plot(); sys.exit()
    if "a" in a.parts:
        rows = part_a()
        write_rows("exp08_orders", rows, __file__)
        print(f"{'case':10s} {'t(L2)':>6} {'t(grad)':>8} {'t(sup)':>7} it(coarsest/finest)")
        for name in sorted(set(x["case"] for x in rows), key=list(CASES).index):
            R = [x for x in rows if x["case"] == name]; Rs = sorted(R, key=lambda x: -x["h"])
            print(f"{name:10s} {fit(R,'L2_rel'):6.2f} {fit(R,'gradL2_rel'):8.2f} {fit(R,'sup_nodes'):7.2f} "
                  f"{Rs[0]['iterations']}/{Rs[-1]['iterations']} conv={[x['converged'] for x in Rs]}")
    if "b" in a.parts:
        rows = part_a(hs=[0.2, 0.1, 0.05, 0.025], names=["p1.2", "p1.1"], max_iter=50000)
        write_rows("exp08_p_to_1", rows, __file__)
    if "c" in a.parts:
        write_rows("exp08_manufactured", part_c(), __file__)
    if "d" in a.parts:
        write_rows("exp08_deadcore_r", part_d(), __file__)
    if "e" in a.parts:
        write_rows("exp08_twophase_holder", part_e(), __file__)
    if "f" in a.parts:
        write_json("exp08_lyapunov", part_f(), __file__)
