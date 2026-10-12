"""Two-phase problem with a HÖLDER coefficient:  G(x,t) = t^p/p + a(x) t^q/q,  a(x) = dist(x, S)^alpha,
S = closed first quadrant {x >= 0, y >= 0} (so a = 0 on the quarter-disk sector, a in C^{0,alpha}, alpha <= 1).

Why this example: with p <= q <= p + alpha (here 2 <= 3 <= 2 + 1) the coefficient condition of the double-phase regularity theory
(Colombo-Mingione; Baroni-Colombo-Mingione; De Filippis-Mingione, interior estimates) and of the no-Lavrentiev-gap theory
(Balci-Gwiazda-Skrzeczkowski) holds, in contrast with the indicator coefficient of ``two_phase.py'' (a is not Hölder).
Nothing is claimed up to the boundary (see the paper).  No exact solution is available: we monitor mesh stability of
max|grad u_h|, the energy, max u_h.

Run:  python examples/two_phase_holder.py [--p 2] [--q 3] [--alpha 1] [--c 4] [--plot]
"""
import argparse, sys
import numpy as np
sys.path.insert(0, __file__.rsplit("/examples", 1)[0])
import glaplace as gl

ap = argparse.ArgumentParser()
ap.add_argument("--p", type=float, default=2.0)
ap.add_argument("--q", type=float, default=3.0)
ap.add_argument("--alpha", type=float, default=1.0)
ap.add_argument("--c", type=float, default=4.0, help="load f = c")
ap.add_argument("--r", type=float, default=3.0)
ap.add_argument("--hs", type=float, nargs="+", default=[0.2, 0.1, 0.05, 0.025])
ap.add_argument("--plot", action="store_true")
a = ap.parse_args()

dist = lambda x, y: np.hypot(np.minimum(x, 0.0), np.minimum(y, 0.0))      # distance to the closed first quadrant
aH = lambda x, y: dist(x, y) ** a.alpha
aI = lambda x, y: np.where((x >= 0) & (y >= 0), 0.0, 1.0)                  # indicator (not Hölder), for comparison

print(f"p={a.p:g} q={a.q:g} alpha={a.alpha:g}  (q <= p+alpha: {a.q <= a.p + a.alpha})  f={a.c:g}  r={a.r:g}")
print(f"{'a':9s} {'h':>6s} {'nodes':>6s} {'it':>5s} {'conv':>5s} {'max u':>9s} {'max|grad u|':>12s} {'||grad u||_L2':>14s} {'J_h(u_h)':>11s}")
out = {}
for name, aa in (("Holder", aH), ("indicator", aI)):
    for h in a.hs:
        mesh = gl.disk_mesh(h)
        S = gl.GLaplaceSolver(mesh, gl.pq(a.p, a.q, aa), a.c, 0.0, r=a.r)
        res = S.solve(tol=1e-9, max_iter=6000, stop="residual")
        fem = S.fem
        gx, gy = fem.gradient(res.u)
        gn = np.hypot(gx, gy)
        l2 = np.sqrt(np.sum(fem.area * gn**2))
        print(f"{name:9s} {h:6.3f} {mesh.n_nodes:6d} {res.iterations:5d} {str(res.converged):>5s} {res.u.max():9.5f} {gn.max():12.5f} {l2:14.6f} {S.energy(res.u):11.6f}", flush=True)
        out[(name, h)] = (mesh, res, gn, S)

if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    mesh, res, gn, S = out[("Holder", min(a.hs))]
    fem = S.fem
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    mesh.plot(res.u, ax=ax[0]); ax[0].set_title(f"u  (a = dist^{a.alpha:g}, p={a.p:g}, q={a.q:g})")
    c = ax[1].tripcolor(mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.triangles, facecolors=gn)
    plt.colorbar(c, ax=ax[1]); ax[1].set_title("|grad u|"); ax[1].set_aspect("equal")
    c = ax[2].tripcolor(mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.triangles, facecolors=aH(fem.bary[:, 0], fem.bary[:, 1]))
    plt.colorbar(c, ax=ax[2]); ax[2].set_title("a(x) = dist(x, S)^alpha"); ax[2].set_aspect("equal")
    plt.tight_layout(); plt.savefig("docs/two_phase_holder.png", dpi=110)
