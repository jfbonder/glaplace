"""EXPERIMENT 3 (NA request 4c; Balci-Ortner-Storn 2022, Sect. 4.2 'bad quadrature'): p(x) checkerboard, lam = 5.

The exponent on the triangles having the origin as a vertex (the nodal patch of 0; 8 triangles of crisscross_mesh(n)) is
replaced by p_max = 3 (variant 'pmax') or p_min = 3/2 ('pmin'); 'exact' = unperturbed.  P1 and CR, n = 8 ... 256.
CAUTION (NA): the perturbed integrands differ from the original one on the patch, so their energies are NOT bounds for the
original m_W or m_H; D(sigma) of each run is a lower bound for the m_W of ITS OWN (perturbed) integrand only.

    python experiments_v4/exp03_bos_quadrature.py [--plot]
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import run_lav, write_rows, read_rows, style, savefig, C, MK
from glaplace.elementwise import PowerSum
from glaplace.lavrentiev import p_zhikov


def G_variant(kind):
    def make(d):
        p = p_zhikov()(d.xb, d.yb).astype(float)
        V, T = d.mesh.vertices, d.mesh.triangles
        patch = np.any(np.linalg.norm(V[T], axis=2) < 1e-14, axis=1)
        if kind == "pmax":
            p[patch] = 3.0
        elif kind == "pmin":
            p[patch] = 1.5
        return PowerSum(np.ones((1, d.Nt)), p[None])
    return make


def compute():
    rows = []
    for kind in ("exact", "pmax", "pmin"):
        for n in (8, 16, 32, 64, 128, 256):
            for s in ("P1", "CR"):
                row = run_lav("px", 5.0, n, s, r=3.0, G=G_variant(kind))
                row["variant"] = kind
                rows.append(row)
                print(f"{kind} n={n} {s} it={row['iterations']} conv={row['converged']} J={row['J']:.8f} D={row['D']:.8f} ({row['time']:.0f}s)", flush=True)
    write_rows("exp03_bos_quadrature", rows, __file__, lam=5.0)


def plot():
    plt = style()
    rows = read_rows("exp03_bos_quadrature")
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    for j, kind in enumerate(("exact", "pmax", "pmin")):
        for s, ls in (("P1", "--"), ("CR", "-")):
            R = sorted([x for x in rows if x["variant"] == kind and x["scheme"] == s], key=lambda x: x["n"])
            lab = {"exact": "unperturbed", "pmax": r"$p_{\max}$ on patch", "pmin": r"$p_{\min}$ on patch"}[kind]
            ax.semilogx([x["ndof"] for x in R], [x["J"] for x in R], ls, color=C[j], marker=MK[j],
                        mfc="white" if s == "P1" else C[j], label=f"{s}, {lab}")
    ax.set_xlabel("number of free dofs"); ax.set_ylabel(r"discrete energy $J_h(u_h)$")
    ax.set_title(r"$p(x)$ checkerboard, $\lambda=5$: perturbed quadrature at the origin")
    ax.legend(fontsize=6.5, ncol=1)
    fig.tight_layout(); savefig(fig, "F7_bos_quadrature.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--plot", action="store_true"); a = ap.parse_args()
    plot() if a.plot else compute()
