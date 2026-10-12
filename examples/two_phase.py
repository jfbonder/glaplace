"""Two-phase problem  -div( (|grad u| + a(x)|grad u|^{q-1}) grad u/|grad u| ) = f  i.e.  G(x,t) = t^p/p + a(x) t^q/q,
with a(x) >= 0 vanishing on a sector (p-Laplacian there, p-q-Laplacian elsewhere).

Domain: unit disk, a(x,y) = 0 if the polar angle is in (0, pi/2) (a sector), 1 elsewhere.  f = 1, u = 0 on the boundary.
No explicit solution; we check the discrete minimality and the effect of a.  A radial two-phase problem with exact solution is in
``examples/convergence_order.py`` (case ``twophase``).

Run:  python examples/two_phase.py [--p 2] [--q 6] [--h 0.03] [--plot]
"""
import argparse, sys
import numpy as np
sys.path.insert(0, __file__.rsplit("/examples", 1)[0])
import glaplace as gl

ap = argparse.ArgumentParser()
ap.add_argument("--p", type=float, default=2.0)
ap.add_argument("--q", type=float, default=6.0)
ap.add_argument("--c", type=float, default=4.0, help="load f = c")
ap.add_argument("--h", type=float, default=0.03)
ap.add_argument("--r", type=float, default=3.0, help="augmented-Lagrangian parameter")
ap.add_argument("--plot", action="store_true")
a = ap.parse_args()

ang = lambda x, y: np.mod(np.arctan2(y, x), 2 * np.pi)
aw = lambda x, y: np.where((ang(x, y) > 0) & (ang(x, y) < np.pi / 2), 0.0, 1.0)
G = gl.pq(a.p, a.q, aw)
mesh = gl.disk_mesh(a.h)
S = gl.GLaplaceSolver(mesh, G, a.c, 0.0, r=a.r)
res = S.solve(tol=1e-9, max_iter=4000, stop="residual")
print(f"{mesh.n_nodes} nodes: {res.iterations} it, converged={res.converged}, {res.time_iter:.1f}s, max u={res.u.max():.5f}")
# reference: the same problem with a = 1 everywhere and with the pure p-Laplacian (a = 0 everywhere)
for name, GG in (("a = 1 everywhere", gl.pq(a.p, a.q)), ("a = 0 everywhere", gl.PowerLaw(a.p))):
    r2 = gl.GLaplaceSolver(mesh, GG, a.c, 0.0, r=a.r).solve(tol=1e-9, max_iter=4000, stop="residual")
    print(f"  {name}: max u = {r2.u.max():.5f} ({r2.iterations} it)")

if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fem = gl.P1Space(mesh)
    gx, gy = fem.gradient(res.u)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    mesh.plot(res.u, ax=ax[0]); ax[0].set_title(f"u  (two-phase, p={a.p:g}, q={a.q:g}, f={a.c:g})")
    c = ax[1].tripcolor(fem.bary[:, 0], fem.bary[:, 1], np.hypot(gx, gy), shading="gouraud") if False else \
        ax[1].tripcolor(mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.triangles, facecolors=np.hypot(gx, gy))
    plt.colorbar(c, ax=ax[1]); ax[1].set_title("|grad u|"); ax[1].set_aspect("equal")
    c = ax[2].tripcolor(mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.triangles, facecolors=aw(fem.bary[:, 0], fem.bary[:, 1]))
    ax[2].set_title("a(x)  (0 on the sector)"); ax[2].set_aspect("equal")
    plt.tight_layout(); plt.savefig("docs/two_phase.png", dpi=110)
