"""EXPERIMENT 2 (NA request 5; Theorem C): iterations vs the augmented-Lagrangian parameter r, for CR and P1.

Cases: p(x) checkerboard lam = 0.5 and 3; two-phase lam = 3; n = 32, 64 (and 128 for the best r of each case/scheme);
r in {0.1, 0.3, 1, 3, 10, 30, 100}; nu_0 = eta_1 = 0; stop: residual < 1e-10 and duality gap < 1e-9 rel; max_iter 50000.
Histories (residual, gap, D(sigma_n)) for p(x) lam = 3, n = 64, CR, r in {0.3, 3, 30} (figure F4).
Lyapunov functional E_n = ||eta_{n+1}-eta*||^2/(2r) + r||nu_n-nu*||^2/2 (Theorem C(iii)) for CR, p(x) lam = 3, n = 32,
r in {1, 3, 10}, with (nu*, eta*) from a run to residual 1e-13 (proxy for the discrete saddle point).

    python experiments_v4/exp02_r_study.py          -> results_v4/exp02_r_study.{csv,json}, exp02_histories.json
    python experiments_v4/exp02_r_study.py --plot   -> figures_v4/F3_iterations_vs_r.pdf, F4_history.pdf, F4b_lyapunov.pdf
"""
import argparse, os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import run_lav, write_rows, write_json, read_rows, style, savefig, C, MK, RESULTS, lav_disc, PROBLEMS
from glaplace.admm import ADMM
from glaplace.v4tools import lyapunov

RS = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
CASES = [("px", 0.5), ("px", 3.0), ("tp", 3.0)]


def compute():
    rows = []
    for p, lam in CASES:
        for n in (32, 64):
            for s in ("CR", "P1"):
                for r in RS:
                    row = run_lav(p, lam, n, s, r=r, check_every=1, jumps=False)
                    rows.append(row)
                    print(f"{p} lam={lam} n={n} {s} r={r:6g} it={row['iterations']:6d} conv={row['converged']} J={row['J']:.10f} "
                          f"({row['time']:.1f}s)", flush=True)
        for s in ("CR", "P1"):
            best = min((x for x in rows if x["problem"] == p and x["lam"] == lam and x["n"] == 64 and x["scheme"] == s),
                       key=lambda x: x["iterations"])
            row = run_lav(p, lam, 128, s, r=best["r"], check_every=5, jumps=False)
            rows.append(row)
            print(f"{p} lam={lam} n=128 {s} best r={best['r']} it={row['iterations']}", flush=True)
    write_rows("exp02_r_study", rows, __file__)
    # histories (F4)
    H = {}
    for r in (0.3, 3.0, 30.0):
        row = run_lav("px", 3.0, 64, "CR", r=r, check_every=1, jumps=True, record=True)
        h = row["_out"]["history"]
        H[str(r)] = {k: h[k].tolist() for k in ("it", "res", "dres", "J", "D", "gap", "jump")}
    # Lyapunov (F4b)
    L = {}
    for r in (1.0, 3.0, 10.0):
        m, d = lav_disc(32, 3.0, "CR")
        A = ADMM(d, PROBLEMS["px"](), r=r)
        ref = A.solve(tol_res=1e-14, tol_gap=1e-14, max_iter=20000, record=False)
        nus, etas = ref["nu"], ref["eta"]
        E = []
        A.solve(tol_res=0, tol_gap=0, max_iter=min(ref["iterations"], 3000), record=False,
                callback=lambda it, st: E.append(lyapunov(d.area, r, st["nu"], st["eta"], nus, etas)))
        L[str(r)] = dict(E=E, ref_iterations=ref["iterations"], ref_residual=ref["residual"])
        print("lyapunov r", r, "ref it", ref["iterations"], "res", ref["residual"], "E0", E[0], "monotone",
              bool(np.all(np.diff(E) <= 1e-14 * E[0])), flush=True)
    write_json("exp02_histories", dict(history=H, lyapunov=L), __file__)


def plot():
    plt = style()
    rows = read_rows("exp02_r_study")
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), sharey=True)
    for ax, (p, lam) in zip(axs, CASES):
        for j, n in enumerate((32, 64)):
            for s, ls in (("CR", "-"), ("P1", "--")):
                R = sorted([x for x in rows if x["problem"] == p and x["lam"] == lam and x["n"] == n and x["scheme"] == s],
                           key=lambda x: x["r"])
                ax.loglog([x["r"] for x in R], [x["iterations"] for x in R], ls, color=C[j], marker=MK[j],
                          mfc="white" if s == "P1" else C[j], label=f"{s}, $n={n}$")
                for x in R:
                    if not x["converged"]:
                        ax.plot(x["r"], x["iterations"], "x", color="k", ms=7)
        ax.set_title((r"$p(x)$" if p == "px" else "two-phase") + rf", $\lambda={lam:g}$")
        ax.set_xlabel("$r$")
    axs[0].set_ylabel("iterations to tolerance")
    axs[0].legend(loc="upper center")
    fig.tight_layout(); savefig(fig, "F3_iterations_vs_r.pdf")

    H = json.load(open(os.path.join(RESULTS, "exp02_histories.json")))["data"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for j, (r, h) in enumerate(H["history"].items()):
        it = np.array(h["it"])
        ls = ["-", "--", "-."][j]
        axs[0].semilogy(it, h["res"], color=C[j], ls=ls, label=f"$r={float(r):g}$")
        axs[1].semilogy(it, np.maximum(np.abs(h["gap"]), 1e-16), color=C[j], ls=ls, label=f"$r={float(r):g}$")
    axs[0].set_ylabel(r"scaled primal residual"); axs[1].set_ylabel(r"$|J_h(u_n)-D(\sigma_n)|$")
    for ax in axs:
        ax.set_xlabel("iteration $n$"); ax.legend()
    axs[0].axhline(1e-10, color="gray", lw=0.7, ls=":");
    fig.suptitle(r"CR, $p(x)$ checkerboard, $\lambda=3$, $n=64$", fontsize=9)
    fig.tight_layout(); savefig(fig, "F4_history.pdf")

    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    for j, (r, d) in enumerate(H["lyapunov"].items()):
        E = np.array(d["E"]); k = E / E[0] > 1e-24
        ax.semilogy(np.arange(1, len(E) + 1)[k], (E / E[0])[k], color=C[j], ls=["-", "--", "-."][j], label=f"$r={float(r):g}$")
    ax.set_xlabel("iteration $n$"); ax.set_ylabel(r"$E_n/E_1$"); ax.legend()
    ax.set_title(r"Lyapunov functional, CR, $p(x)$, $\lambda=3$, $n=32$")
    fig.tight_layout(); savefig(fig, "F4b_lyapunov.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--plot", action="store_true"); a = ap.parse_args()
    plot() if a.plot else compute()
