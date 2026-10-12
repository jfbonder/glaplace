"""Manufactured solution with an x-dependent Young function (needs sympy):

    G(x,y,t) = t^2/2 + (1+x) t^4/4,    u = exp(x+y)    on (0,1)^2,    f = -div(g(x,|grad u|) grad u/|grad u|)  (symbolic).

Reports true L2 / gradient-L2 errors and the fitted order, and the iterations needed for several values of the
augmented-Lagrangian parameter r (here |grad u| is up to ~ e^2 * sqrt(2), so g(t)/t ~ t^2 is large and r = 1 is slow).

Run:  python examples/manufactured.py [--n 10 20 40 80]
"""
import argparse, sys
import numpy as np
sys.path.insert(0, __file__.rsplit("/examples", 1)[0])
import glaplace as gl

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, nargs="+", default=[10, 20, 40, 80])
ap.add_argument("--r", type=float, default=100.0)
a = ap.parse_args()

f, u, gu = gl.exact.manufactured("t**2/2 + (1+x)*t**4/4", "exp(x+y)")
G = gl.pq(2, 4, lambda x, y: 1 + x)
print(f"r = {a.r}")
print(f"{'n':>4} {'nodes':>7} {'it':>5} {'L2 rel':>10} {'order':>6} {'grad L2 rel':>12} {'order':>6}")
prev = None
for n in a.n:
    m = gl.rectangle_mesh(n, n)
    res = gl.GLaplaceSolver(m, G, f, u, r=a.r, f_degree=2).solve(tol=1e-10, max_iter=5000, stop="residual")
    F = gl.P1Space(m)
    e2, eg = F.l2_error(res.u, u), F.grad_error_Lp(res.u, gu, 2.0)
    o = ("", "") if prev is None else (f"{np.log(prev[0] / e2) / np.log(2):.2f}", f"{np.log(prev[1] / eg) / np.log(2):.2f}")
    print(f"{n:4d} {m.n_nodes:7d} {res.iterations:5d} {e2:10.3e} {o[0]:>6} {eg:12.3e} {o[1]:>6}")
    prev = (e2, eg)

print("\nEffect of r (n = 40):")
m = gl.rectangle_mesh(40, 40)
for r in (1, 10, 30, 100, 300):
    res = gl.GLaplaceSolver(m, G, f, u, r=float(r), f_degree=2).solve(tol=1e-10, max_iter=5000, stop="residual")
    print(f"  r = {r:4d}: {res.iterations:5d} iterations, converged={res.converged}")
