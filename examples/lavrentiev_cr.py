"""P1 (conforming) vs Crouzeix-Raviart (nonconforming) on Zhikov-type problems with a Lavrentiev gap.

    python examples/lavrentiev_cr.py [--problem pxlap|twophase] [--lam 3] [--ns 8 16 32 64] [--r 3] [--plot]

Omega = (-1,1)^2, f = 0, boundary data lam*u0 (see glaplace/lavrentiev.py).  Prints F_h(u_h) for both schemes.
For lam = 1 and the p(x) problem the exact minimum is 2; for lam = 3 the P1 energy stays above the CR energy (gap)."""
import argparse, sys
import numpy as np
sys.path.insert(0, __file__.rsplit("/examples", 1)[0])
from glaplace.lavrentiev import compare, G_pxlap, G_twophase

ap = argparse.ArgumentParser()
ap.add_argument("--problem", choices=["pxlap", "twophase"], default="pxlap")
ap.add_argument("--lam", type=float, default=3.0)
ap.add_argument("--ns", type=int, nargs="+", default=[8, 16, 32, 64])
ap.add_argument("--r", type=float, default=3.0)
ap.add_argument("--plot", action="store_true")
a = ap.parse_args()
G = G_pxlap() if a.problem == "pxlap" else G_twophase()
rows = compare(G, a.lam, ns=a.ns, r=a.r)
if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    for s, c in (("P1", "C3"), ("CR", "C0")):
        R = [x for x in rows if x["scheme"] == s]
        plt.plot([1 / np.sqrt(x["ndof"]) for x in R], [x["energy"] for x in R], "o-", c=c, label=s)
    plt.xlabel("ndof$^{-1/2}$"); plt.ylabel("$F_h(u_h)$"); plt.legend(); plt.title(f"{a.problem}, lambda={a.lam:g}")
    plt.savefig(f"docs/lavrentiev_{a.problem}_lam{a.lam:g}.png", dpi=110)
