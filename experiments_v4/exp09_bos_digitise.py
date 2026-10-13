"""EXPERIMENT 9 (NA request 4d): digitise Balci-Ortner-Storn (arXiv 2106.06837v1) Figure 1 (vector graphics, p. 16):
F(u_c) - F_h(u_nc) vs ndof (adaptive meshes, Newton) for Experiment 1 (p(x), lam = 0.99 ... 1.08, left) and
Experiment 3 (two-phase, lam = 0.3 ... 0.6, right), and overlay our widths Delta_h (uniform criss-cross, ADMM).

The polylines are read with pdfplumber; axes are calibrated from the tick marks (left: x = 164.2 + 25.35 log10(ndof/10),
1e-3 at y = 132.26, 14.61 pt/decade; right: 10^1 at x = 336.68 (label centres), 25.33 pt/decade, 1e-2 at y = 129.92,
15.97 pt/decade).  Colour -> lambda assignment follows the legend order (black, red, green, blue, [yellow]); this is an
ASSUMPTION (the legend entries are drawn in that order), consistent with the decay of the lowest (black) curve for
lam = 0.99.  Digitisation accuracy ~0.05 decades.

    python experiments_v4/exp09_bos_digitise.py        -> results_v4/exp09_bos_fig1.csv, figures_v4/F6c_bos_overlay.pdf
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import write_rows, read_rows, style, savefig, C, MK

import os
PDF = os.environ.get("GL_BOS_PDF", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "biblio", "2106.06837v1.pdf"))  # local copy of Balci-Ortner-Storn (not in the repo)


def digitise():
    import pdfplumber
    pg = pdfplumber.open(PDF).pages[15]
    curves = [c for c in pg.curves if len(c["pts"]) >= 25 and c["top"] < 240]
    left = [c for c in curves if c["x0"] < 300]; right = [c for c in curves if c["x0"] > 300]
    rows = []
    for side, cs, lams, cal in (("Exp1_px", left, [0.99, 1.02, 1.04, 1.06, 1.08], (164.2, 25.35, 132.26, 14.61, -3)),
                                ("Exp3_tp", right, [0.3, 0.4, 0.5, 0.6], (336.68, 25.33, 129.92, 15.97, -2))):
        x0, dx, y0, dy, e0 = cal
        for c, lam in zip(cs, lams):
            for (x, y) in c["pts"]:
                rows.append(dict(experiment=side, lam=lam, colour=str(c["stroking_color"]),
                                 ndof=10 ** (1 + (x - x0) / dx), value=10 ** (e0 - (y - y0) / dy)))
    return rows


if __name__ == "__main__":
    rows = digitise()
    write_rows("exp09_bos_fig1", rows, __file__, source="BOS arXiv 2106.06837v1 Fig. 1, digitised")
    ours = read_rows("exp01_lavrentiev")
    plt = style()
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=False)
    for ax, side, p, lams in ((axs[0], "Exp1_px", "px", [0.99, 1.02, 1.04, 1.06, 1.08]), (axs[1], "Exp3_tp", "tp", [0.3, 0.4, 0.5, 0.6])):
        for j, lam in enumerate(lams):
            B = [r for r in rows if r["experiment"] == side and r["lam"] == lam]
            ax.loglog([r["ndof"] for r in B], [r["value"] for r in B], "-", color=C[j], lw=1.0, alpha=0.8)
            P = {x["n"]: x for x in ours if x["problem"] == p and abs(x["lam"] - lam) < 1e-12 and x["scheme"] == "P1"}
            Q = {x["n"]: x for x in ours if x["problem"] == p and abs(x["lam"] - lam) < 1e-12 and x["scheme"] == "CR"}
            ns = sorted(set(P) & set(Q))
            ax.loglog([Q[n]["ndof"] for n in ns], [max(P[n]["J"] - Q[n]["best_D"], 1e-16) for n in ns], "--", color=C[j],
                      marker=MK[j], label=rf"$\lambda={lam:g}$")
        ax.set_xlabel("ndof"); ax.set_ylabel(r"$J_h^{P1}-J_h^{CR}$")
        ax.set_title(("BOS Exp. 1 ($p(x)$)" if p == "px" else "BOS Exp. 3 (two-phase)") + ": solid BOS (adaptive), dashed ours", fontsize=8)
        ax.legend(fontsize=6.5)
    fig.tight_layout(); savefig(fig, "F6c_bos_overlay.pdf")
