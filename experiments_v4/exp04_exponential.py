"""EXPERIMENT 4 (NA request 6): exponential integrands on the square (no BOS theory; Theorem D still applies).

(a) autonomous G = e^t - 1 - t, data lam*u0, lam in {0.5, 3}: the minimiser is lam*x2 (affine), so P1 = CR = D exactly
    (code test only).
(b) EXPLORATORY (no theory): G(x,t) = a(x)(e^t - 1 - t) + t^2/2, a = 1 on |x1| < |x2| (checkerboard), lam in {0.5, 1, 2, 3}.
Criss-cross meshes n = 8 ... 256, r = 3, same stopping test as experiment 1.

    python experiments_v4/exp04_exponential.py [--plot]
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import run_lav, write_rows, read_rows, style, savefig, C, MK
from glaplace.elementwise import ExpMinus1MinusT
from glaplace.lavrentiev import a_checker


def G_auto(d):
    return ExpMinus1MinusT(d.Nt)


def G_x(d):
    return ExpMinus1MinusT(d.Nt, w=a_checker(d.xb, d.yb), beta=1.0)


def compute():
    rows = []
    for case, G, lams, ns in (("auto", G_auto, (0.5, 3.0), (8, 16, 32, 64, 128)),
                              ("xdep", G_x, (0.5, 1.0, 2.0, 3.0), (8, 16, 32, 64, 128, 256))):
        for lam in lams:
            for n in ns:
                for s in ("P1", "CR"):
                    row = run_lav("px", lam, n, s, r=3.0, G=G)
                    row["problem"] = "exp_" + case
                    rows.append(row)
                    print(f"{case} lam={lam} n={n} {s} it={row['iterations']} conv={row['converged']} J={row['J']:.10f} D={row['D']:.10f} ({row['time']:.0f}s)", flush=True)
    write_rows("exp04_exponential", rows, __file__)


def plot():
    plt = style()
    rows = read_rows("exp04_exponential")
    fig, ax = plt.subplots(figsize=(3.8, 2.8))
    for j, lam in enumerate((0.5, 1.0, 2.0, 3.0)):
        P = {x["n"]: x for x in rows if x["problem"] == "exp_xdep" and x["lam"] == lam and x["scheme"] == "P1"}
        Q = {x["n"]: x for x in rows if x["problem"] == "exp_xdep" and x["lam"] == lam and x["scheme"] == "CR"}
        ns = sorted(set(P) & set(Q))
        ax.loglog([Q[n]["ndof"] for n in ns], [P[n]["J"] - Q[n]["best_D"] for n in ns], color=C[j], marker=MK[j],
                  label=rf"$\lambda={lam:g}$")
    x = np.array([1e2, 1e6]); ax.loglog(x, 30 / x, "k:", lw=0.8, label=r"$\propto\mathrm{ndof}^{-1}$")
    ax.set_xlabel("CR free dofs"); ax.set_ylabel(r"width $\Delta_h=J_h^{P1}-D(\sigma_h)$")
    ax.set_title(r"$a(x)(e^t-1-t)+t^2/2$ (exploratory)"); ax.legend()
    fig.tight_layout(); savefig(fig, "F8_exponential_width.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--plot", action="store_true"); a = ap.parse_args()
    plot() if a.plot else compute()
