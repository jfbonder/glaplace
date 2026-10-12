"""G(t) = e^t - 1  (g = e^t, g(0+) = 1: no Delta_2 condition and g(0+) > 0).

Torsion on the unit disk with load c:  g(|u'|) = c r/2 has a positive solution only where c r/2 > 1, so u is flat ("dead core")
for r <= 2/c.  With c = 8 the core is r < 1/4.  We compare with the exact solution and study the augmented-Lagrangian parameter r
(the iteration is slow when g(0+) > 0, like for the 1-Laplacian, and r matters).

Run:  python examples/exp_dead_core.py [--plot]
"""
import argparse, sys
import numpy as np
sys.path.insert(0, __file__.rsplit("/examples", 1)[0])
import glaplace as gl

ap = argparse.ArgumentParser()
ap.add_argument("--c", type=float, default=8.0)
ap.add_argument("--h", type=float, default=0.05)
ap.add_argument("--plot", action="store_true")
a = ap.parse_args()

G = gl.ExpLiteral()
u, gu, ur = gl.exact.radial_torsion(G, a.c)
mesh = gl.disk_mesh(a.h); fem = gl.P1Space(mesh)
print(f"c = {a.c}: dead core radius 2/c = {2 / a.c:.3f},  u(0) = {u(0., 0.):.6f}")
print(f"{'r':>6} {'iterations':>11} {'converged':>10} {'rel L2 err':>11}")
best = None
for r in (1.0, 3.0, 10.0, 30.0, 100.0):
    res = gl.GLaplaceSolver(mesh, G, a.c, 0.0, r=r).solve(tol=1e-8, max_iter=6000, stop="residual")
    print(f"{r:6.1f} {res.iterations:11d} {str(res.converged):>10} {fem.l2_error(res.u, u):11.3e}")
    if best is None or res.iterations < best[1].iterations:
        best = (r, res)
r, res = best
gx, gy = fem.gradient(res.u)
core = np.hypot(*fem.bary.T) < 2 / a.c - 0.05
print(f"best r = {r}: max |grad u| in the core (r < 2/c - 0.05) = {np.hypot(gx, gy)[core].max():.2e}  (exact: 0)")
if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rr = np.linspace(0, 1, 400)
    plt.figure(figsize=(6, 4))
    plt.plot(rr, u(rr, 0 * rr), "k", label="exact")
    order = np.argsort(np.hypot(*mesh.vertices.T))
    plt.plot(np.hypot(*mesh.vertices.T)[order], res.u[order], ".", ms=2, label="computed")
    plt.axvline(2 / a.c, ls=":", c="gray"); plt.xlabel("r"); plt.legend(); plt.title(f"G = e^t - 1, c = {a.c:g}: dead core")
    plt.tight_layout(); plt.savefig("docs/exp_dead_core.png", dpi=110)
