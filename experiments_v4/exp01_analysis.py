"""Analysis and figures of EXPERIMENT 1 (reads results_v4/exp01_lavrentiev.json; no computation).

Figures (figures_v4/):
  F1_bracket_lambda3.pdf   NA (F1): P1 energy (upper) and best D(sigma_n) (lower) vs ndof, shaded bracket, lam = 3, both problems
  F2_width_vs_ndof.pdf     NA (F2): width Delta_h = J_h(P1) - best D vs ndof (log-log), all lam, slope ndof^-1 reference
  F5_nogap_controls.pdf    NA 2: no-gap controls (p(x) lam <= 0.75, two-phase lam <= 0.2), with P1 - E_inf and E_inf - CR
  F6_bos_thresholds.pdf    NA 4a-b: widths for lam near the BOS thresholds (p(x): 0.99 ... 1.08; two-phase 0.3 ... 0.6)
  F6b_width_vs_lambda.pdf  width on the finest mesh vs lam
  F7b_diagnostics_checkerboard.pdf   PDE section 0 / E4(i): N2, E2, Einf vs log(1/h), P1 and CR
Tables (results_v4/): exp01_widths.csv (ratios, slopes), exp01_extrapolation.csv (NA 3, HEURISTIC), printed as markdown.

    python experiments_v4/exp01_analysis.py
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_rows, write_rows, style, savefig, C, MK, PROBLEM_LABEL

rows = read_rows("exp01_lavrentiev")


def series(p, lam, s):
    return {x["n"]: x for x in rows if x["problem"] == p and abs(x["lam"] - lam) < 1e-12 and x["scheme"] == s}


def width(p, lam):
    P, Q = series(p, lam, "P1"), series(p, lam, "CR")
    ns = sorted(set(P) & set(Q))
    return ns, np.array([Q[n]["ndof"] for n in ns]), np.array([P[n]["J"] - Q[n]["best_D"] for n in ns]), P, Q


def lams(p):
    return sorted(set(x["lam"] for x in rows if x["problem"] == p))


def aitken(E):            # E ordered coarse -> fine, last three
    e1, e2, e3 = E[-3:]
    with np.errstate(all="ignore"):
        k = np.log2((e2 - e1) / (e3 - e2))
        return k, e3 + (e3 - e2) / (2 ** k - 1)     # NOTE: NA_REQUESTS.md writes "-" (sign typo)


def rich(E, k):
    return E[-1] + (E[-1] - E[-2]) / (2 ** k - 1)


def lsq(h, E):
    f = lambda hh, Einf, Cc, k: Einf + Cc * hh ** k
    try:
        d = E[-1] - E[-2]
        popt, _ = curve_fit(f, h[-4:], E[-4:], p0=(E[-1], d / (h[-1] - h[-2]), 1.0), maxfev=20000,
                            bounds=([-np.inf, -np.inf, 0.05], [np.inf, np.inf, 4.0]))
        return popt[2], popt[0]
    except Exception:
        return np.nan, np.nan


plt = style()
# ---------------- F1
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
for ax, p in zip(axs, ("px", "tp")):
    ns, nd, W, P, Q = width(p, 3.0)
    up = np.array([P[n]["J"] for n in ns]); lo = np.array([Q[n]["best_D"] for n in ns])
    ax.fill_between(nd, lo, up, color=C[5], alpha=0.25, lw=0, label="guaranteed bracket")
    ax.semilogx(nd, up, "-", color=C[1], marker="s", label=r"P1 energy $J_h(u_h^{P1})\geq m_H$")
    ax.semilogx(nd, lo, "-", color=C[0], marker="o", label=r"CR bound $D(\sigma_h)\leq m_W$")
    for i, n in enumerate(ns):
        for S, v in ((P, up[i]), (Q, lo[i])):
            if not S[n]["converged"]:
                ax.plot(nd[i], v, "x", color="k", ms=7)      # flag: stopping test not met
    ax.set_xlabel("CR free dofs"); ax.set_ylabel("energy"); ax.set_title(PROBLEM_LABEL[p] + r", $\lambda=3$")
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.08, 1, 1)); savefig(fig, "F1_bracket_lambda3.pdf")

# ---------------- F2
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
main = {"px": [0.5, 1.5, 2.0, 3.0, 4.0], "tp": [0.1, 0.2, 0.3, 0.5, 1.0, 3.0]}
for ax, p in zip(axs, ("px", "tp")):
    for j, lam in enumerate(main[p]):
        ns, nd, W, P, Q = width(p, lam)
        ax.loglog(nd, np.maximum(W, 1e-14), color=C[j], marker=MK[j], label=rf"$\lambda={lam:g}$")
    if p == "px":
        ns, nd, W, P, Q = width(p, 1.0)
        ax.text(nd[0], 3e-6, r"$\lambda=1$: $\Delta_h\leq%.0e$ (exact)" % max(abs(W).max(), 1e-16), fontsize=7)
    x = np.array([2e2, 1.6e6]); ax.loglog(x, 2e2 / x, "k:", lw=0.9, label=r"$\propto\mathrm{ndof}^{-1}$")
    ax.set_xlabel("CR free dofs"); ax.set_title(PROBLEM_LABEL[p])
axs[0].set_ylabel(r"width $\Delta_h=J_h(u_h^{P1})-D(\sigma_h)\ \geq m_H-m_W$")
for ax in axs:
    ax.legend(fontsize=6.5, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.2))
fig.tight_layout(); savefig(fig, "F2_width_vs_ndof.pdf")

# ---------------- width table (ratios, slopes)
wrows = []
for p in ("px", "tp"):
    for lam in lams(p):
        ns, nd, W, P, Q = width(p, lam)
        if len(ns) < 3:
            continue
        rat = W[:-1] / W[1:]
        k = min(4, len(ns))
        slope = np.polyfit(np.log(nd[-k:]), np.log(np.abs(W[-k:]) + 1e-300), 1)[0]
        wrows.append(dict(problem=p, lam=lam, n_max=ns[-1], width_finest=W[-1], P1_finest=P[ns[-1]]["J"],
                          CR_finest=Q[ns[-1]]["best_D"], ratios=" ".join(f"{x:.2f}" for x in rat), slope_ndof=slope,
                          all_converged=all(P[n]["converged"] and Q[n]["converged"] for n in ns),
                          max_jump=max(Q[n]["max_jump"] for n in ns),
                          min_gap=min(min(P[n]["gap"], Q[n]["gap"]) for n in ns)))
write_rows("exp01_widths", wrows, __file__)
print("\n| problem | lambda | n_max | P1 (upper) | D (lower) | width | ratios W_{n/2}/W_n | slope in ndof | converged |")
print("|---|---|---|---|---|---|---|---|---|")
for w in wrows:
    print(f"| {w['problem']} | {w['lam']:g} | {w['n_max']} | {w['P1_finest']:.8f} | {w['CR_finest']:.8f} | {w['width_finest']:.3e} | "
          f"{w['ratios']} | {w['slope_ndof']:.2f} | {w['all_converged']} |")

# ---------------- extrapolation (heuristic)
erows = []
for p, lam in (("px", 3.0), ("tp", 3.0), ("px", 0.5), ("tp", 0.2), ("px", 2.0)):
    ns, nd, W, P, Q = width(p, lam)
    h = 2.0 / np.array(ns, float)
    for name, E in (("P1", np.array([P[n]["J"] for n in ns])), ("CR", np.array([Q[n]["best_D"] for n in ns])), ("width", W)):
        k, Ea = aitken(E); kl, El = lsq(h, E)
        erows.append(dict(problem=p, lam=lam, quantity=name, n_max=ns[-1], finest=E[-1], aitken_kappa_h=k, aitken=Ea,
                          rich_ndof_m12=rich(E, 1.0), rich_ndof_m1=rich(E, 2.0), lsq_kappa_h=kl, lsq=El))
write_rows("exp01_extrapolation", erows, __file__, note="HEURISTIC extrapolations; only the bracket is a theorem")
print("\n| problem | lambda | quantity | finest | Aitken kappa(h) | Aitken | Rich. ndof^-1/2 | Rich. ndof^-1 | LSQ kappa | LSQ |")
print("|---|---|---|---|---|---|---|---|---|---|")
for e in erows:
    print(f"| {e['problem']} | {e['lam']:g} | {e['quantity']} | {e['finest']:.5f} | {e['aitken_kappa_h']:.2f} | {e['aitken']:.4f} | "
          f"{e['rich_ndof_m12']:.4f} | {e['rich_ndof_m1']:.4f} | {e['lsq_kappa_h']:.2f} | {e['lsq']:.4f} |")

# ---------------- F5 no-gap controls
ctrl = {"px": [0.25, 0.5, 0.75], "tp": [0.1, 0.2]}
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
j = 0
crows = []
for p in ("px", "tp"):
    for lam in ctrl[p]:
        ns, nd, W, P, Q = width(p, lam)
        Pv = np.array([P[n]["J"] for n in ns]); Cv = np.array([Q[n]["best_D"] for n in ns])
        Einf = rich(Cv, 2.0)                       # heuristic limit (Richardson in ndof^-1 on the CR sequence)
        lab = ("$p(x)$" if p == "px" else "two-phase") + rf", $\lambda={lam:g}$"
        axs[0].loglog(nd, W, color=C[j], marker=MK[j], label=lab)
        axs[1].loglog(nd, Pv - Einf, "-", color=C[j], marker=MK[j], label=lab)
        axs[1].loglog(nd[:-1], (Einf - Cv)[:-1], "--", color=C[j], marker=MK[j], mfc="white")
        crows.append(dict(problem=p, lam=lam, E_inf_heuristic=Einf, P1_minus_Einf=" ".join(f"{x:.2e}" for x in Pv - Einf),
                          Einf_minus_CR=" ".join(f"{x:.2e}" for x in Einf - Cv),
                          width_ratios=" ".join(f"{x:.2f}" for x in W[:-1] / W[1:])))
        j += 1
ns, nd, W, P, Q = width("px", 1.0)
axs[0].loglog(nd, np.maximum(np.abs(W), 1e-16), color="k", marker="x", lw=0.8, label=r"$p(x)$, $\lambda=1$ (exact; $|\Delta_h|$)")
for ax in axs:
    x = np.array([2e2, 1.6e6]); ax.loglog(x, 30 / x, "k:", lw=0.9)
    ax.set_xlabel("CR free dofs")
axs[0].set_ylabel(r"width $\Delta_h$ (rigorous bound for the gap)")
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=7, bbox_to_anchor=(0.5, -0.02))
axs[1].set_ylabel(r"$J_h^{P1}-E_\infty$ (solid), $E_\infty-D$ (dashed)")
axs[1].set_title(r"heuristic $E_\infty$ = Richardson (ndof$^{-1}$) of CR", fontsize=8)
fig.tight_layout(rect=(0, 0.13, 1, 1)); savefig(fig, "F5_nogap_controls.pdf")
write_rows("exp01_nogap_controls", crows, __file__, note="E_inf heuristic")
print("\n| problem | lambda | width ratios | P1 - E_inf | E_inf - CR |\n|---|---|---|---|---|")
for c in crows:
    print(f"| {c['problem']} | {c['lam']:g} | {c['width_ratios']} | {c['P1_minus_Einf']} | {c['Einf_minus_CR']} |")

# ---------------- F6 BOS thresholds
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
for ax, p, L in ((axs[0], "px", [0.99, 1.0, 1.02, 1.04, 1.06, 1.08]), (axs[1], "tp", [0.3, 0.4, 0.5, 0.6])):
    for j, lam in enumerate(L):
        ns, nd, W, P, Q = width(p, lam)
        ax.loglog(nd, np.maximum(np.abs(W), 1e-16), color=C[j], marker=MK[j], label=rf"$\lambda={lam:g}$")
    x = np.array([2e2, 1.6e6]); ax.loglog(x, 3 / x, "k:", lw=0.9, label=r"$\propto\mathrm{ndof}^{-1}$")
    ax.set_xlabel("CR free dofs"); ax.set_title(PROBLEM_LABEL[p])
    ax.legend(fontsize=6.5, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.2))
axs[0].set_ylabel(r"width $\Delta_h$")
fig.tight_layout(); savefig(fig, "F6_bos_thresholds.pdf")

fig, ax = plt.subplots(figsize=(3.8, 2.8))
for j, p in enumerate(("px", "tp")):
    L = lams(p); Wf = []; nn = []
    for lam in L:
        ns, nd, W, P, Q = width(p, lam); Wf.append(W[-1]); nn.append(ns[-1])
    ax.semilogy(L, np.maximum(np.abs(Wf), 1e-16), color=C[j], marker=MK[j], label=("$p(x)$" if p == "px" else "two-phase") + f", $n={max(nn)}$")
ax.set_xlabel(r"$\lambda$"); ax.set_ylabel(r"$\Delta_h$ on the finest mesh"); ax.legend()
fig.tight_layout(); savefig(fig, "F6b_width_vs_lambda.pdf")

# ---------------- F7b diagnostics (PDE section 0)
fig, axs = plt.subplots(2, 3, figsize=(7.4, 4.6), sharex=True)
for i, (p, L) in enumerate((("px", [0.5, 1.0, 1.5, 3.0, 4.0]), ("tp", [0.2, 0.5, 1.0, 3.0, 4.0]))):
    for j, lam in enumerate(L):
        for s, ls in (("P1", "-"), ("CR", "--")):
            S = series(p, lam, s); ns = sorted(S)
            for k, key in enumerate(("N2", "E2", "Einf")):
                (axs[i, k].semilogy if key == "Einf" else axs[i, k].plot)(np.log(np.array(ns) / 2.0), [S[n][key] for n in ns], ls, color=C[j], marker=MK[j],
                               mfc="white" if s == "CR" else C[j], label=rf"$\lambda={lam:g}$" if s == "P1" else None)
    for k, t in enumerate((r"$\|\nabla_h u_h\|_{L^2}$", r"$\|\eta_h\|_{L^2}$", r"$\max_T|\eta_h|$")):
        axs[i, k].set_title(("$p(x)$: " if p == "px" else "two-phase: ") + t)
    axs[i, 2].legend(fontsize=6, loc="upper left")
for k in range(3):
    axs[1, k].set_xlabel(r"$\log(1/h)$")
fig.suptitle("Checkerboards: size of the discrete saddle points (solid P1, dashed CR)", fontsize=9)
fig.tight_layout(); savefig(fig, "F7b_diagnostics_checkerboard.pdf")

# diagnostics table for the report (finest four meshes)
print("\n| problem | lambda | scheme | n | N2 | E2 | Einf | E0 | iters |\n|---|---|---|---|---|---|---|---|---|")
for p, L in (("px", [1.0, 3.0]), ("tp", [0.2, 3.0])):
    for lam in L:
        for s in ("P1", "CR"):
            S = series(p, lam, s)
            for n in sorted(S)[-4:]:
                x = S[n]
                print(f"| {p} | {lam:g} | {s} | {n} | {x['N2']:.4f} | {x['E2']:.4f} | {x['Einf']:.3f} | {x['E0']:.2f} | {x['iterations']} |")
