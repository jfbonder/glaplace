"""Figures for the article from the raw data (convergence_order.csv, lyapunov_curves.npz)."""
import csv, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
base = "/home/claude/paper_v3"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#555555", "axes.labelcolor": "#222222", "xtick.color": "#444444",
                     "ytick.color": "#444444", "lines.linewidth": 1.6, "legend.frameon": False})
C = ["#2a78d6", "#eb6834", "#1baf7a", "#d55181", "#c98500"]
rows = list(csv.DictReader(open(f"{base}/materials/convergence_order.csv")))
cases = []
for r in rows:
    if r["case"] not in cases: cases.append(r["case"])
hl = {"p4": (C[0], "o", r"$t^4/4$"), "plog2_1": (C[1], "s", r"$t^2\log(1+t)$"),
      "twophase": (C[2], "^", r"two-phase, $a=\mathbf{1}_{|x|>1/2}$"), "expdead": (C[3], "D", r"$e^t-1$, $c=8$")}
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
keys = [("L2_rel", r"rel. $L^2$ error of $u$", 2), ("gradL2_rel", r"rel. $L^2$ error of $\nabla u$", 1),
        ("sup_nodes", r"max nodal error $|u-u_h|$", None)]
ymin = {k_: min(float(r[k_]) for r in rows) for k_, _, _ in keys}
for k, (key, ttl, slope) in enumerate(keys):
    a = ax[k]
    for c in cases:
        R = sorted([r for r in rows if r["case"] == c], key=lambda r: -float(r["h"]))
        h = np.array([float(r["h"]) for r in R]); y = np.array([float(r[key]) for r in R])
        if c in hl:
            col, mk, lab = hl[c]
            a.loglog(h, y, "-", color=col, marker=mk, ms=4.5, label=lab, zorder=3)
        else:
            a.loglog(h, y, "-", color="#b5b5b0", lw=1.0, marker="o", ms=2.5, zorder=2,
                     label="other smooth $G$ (6 cases)" if c == "p1.5" else None)
    if slope:
        hh = np.array([0.0125, 0.2]); y0 = ymin[key] * 0.4
        a.loglog(hh, y0 * (hh / hh[0]) ** slope, "k--", lw=0.8, label=f"slope {slope}")
    a.set_xlabel("$h$ (target mesh size)"); a.set_title(ttl, fontsize=9)
    a.set_xticks([0.0125, 0.025, 0.05, 0.1, 0.2]); a.set_xticklabels(["0.0125", "0.025", "0.05", "0.1", "0.2"])
    a.minorticks_off(); a.grid(True, which="major", color="#e6e6e3", lw=0.6)
ax[2].legend(fontsize=7.5, loc="lower right")
for j in (0, 1): ax[j].legend(handles=[l for l in ax[j].get_lines() if l.get_label().startswith("slope")], fontsize=7.5, loc="lower right")
plt.tight_layout(); plt.savefig(f"{base}/new/figures/convergence_order.pdf"); plt.close()

d = np.load(f"{base}/new/results/lyapunov_curves.npz")
lab = {"pq24": r"$t^2/2+t^4/4$, $r=1$", "p1.5": r"$t^{1.5}/1.5$, $r=1$", "p4": r"$t^4/4$, $r=1$",
       "exp": r"$e^t-1-t$, $r=1$", "expdead": r"$e^t-1$ ($c=8$), $r=30$"}
fig, a = plt.subplots(figsize=(5.2, 3.2))
for i, k in enumerate(["pq24", "exp", "p4", "p1.5", "expdead"]):
    y = d[k]; n = np.arange(1, len(y) + 1)
    m = (y > 1e-24) & (n <= 800)
    a.semilogy(n[m], y[m], color=C[i], label=lab[k])
a.set_xlabel("iteration $n$"); a.set_ylabel(r"$E_n/E_0$"); a.grid(True, color="#e6e6e3", lw=0.6)
a.legend(fontsize=7.5); plt.tight_layout(); plt.savefig(f"{base}/new/figures/lyapunov.pdf"); plt.close()
print("ok")
