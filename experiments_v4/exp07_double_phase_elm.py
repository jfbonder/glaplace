"""EXPERIMENT 7 (PDE request E3): double phase with the Esposito-Leonetti-Mingione weight, in / out of the no-gap range.

Omega = (-1,1)^2, f = 0, data lam*u0, G = t^p/p + a(x) t^q/q, p = 1.5,
a(x) = max(x2^2 - x1^2, 0)/|x|^(2-alpha), alpha = 1 (Lipschitz, vanishes on the cones |x2| <= |x1|).
  q = 2.4 (<= p + alpha = 2.5: no gap, BGS Thm 2.3);  q = 2.8 (grey zone);  q = 4 (> d + alpha = 3: ELM Thm 3 gap).
lam in {0.5, 1, 2, 4, 8}; crisscross meshes n = 8 ... 128 (256 for lam in {1, 8}, P1-max and CR-min only); r = 3; stop as in experiment 1.
One-point quadratures of a on each triangle T:
  P1-bary  a(x_T)                 (PDE request; no bound proved)
  P1-max   max_T a  (majorant)    -> J_h(u_h) >= m_H >= m_W   (NA Prop. P1 with majorant sampling)
  CR-bary  a(x_T)                 (BOS Thm 1 NOT covered)
  CR-min   min_T a  (minorant)    -> D(sigma) <= m_W (Theorem D, Remark variants (a)); BOS assumption (B)
min_T a, max_T a: glaplace.v4tools.element_minmax (extrema on the edges, located to ~1e-8 of the edge length; this
introduces a relative error ~1e-15 in the min/max values -- immaterial but strictly not 'rigorous').
Rigorous bracket: [best D(CR-min), J_h(P1-max)].

    python experiments_v4/exp07_double_phase_elm.py [--jobs 1] [--plot]
"""
import argparse, os, sys, json, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import run_lav, write_rows, read_rows, style, savefig, C, MK, RESULTS, CODE_HASH
from glaplace.elementwise import PowerSum
from glaplace.v4tools import a_elm, element_minmax

QS = (2.4, 2.8, 4.0)
LAMS = (0.5, 1.0, 2.0, 4.0, 8.0)
VARIANTS = (("P1", "bary"), ("P1", "max"), ("CR", "bary"), ("CR", "min"))
CACHE = os.path.join(RESULTS, "exp07_cache.jsonl")


_MINMAX = {}


def _minmax(mesh):
    # element-wise min/max of a depend only on the mesh: computed once per mesh size and reused for all q, lam
    key = mesh.n_elems
    if key not in _MINMAX:
        _MINMAX[key] = element_minmax(a_elm(1.0), mesh)
    return _MINMAX[key]


def Gq(q, kind):
    def mk(d):
        a = a_elm(1.0)
        if kind == "bary":
            aT = a(d.xb, d.yb)
        else:
            mn, mx = _minmax(d.mesh); aT = mn if kind == "min" else mx
        return PowerSum(np.vstack([np.ones(d.Nt), aT]), np.vstack([np.full(d.Nt, 1.5), np.full(d.Nt, q)]))
    return mk


def job(args):
    q, lam, n, s, kind = args
    row = run_lav("px", lam, n, s, r=3.0, G=Gq(q, kind), check_every=1 if n <= 32 else 5)
    row.update(problem="elm", q=q, quad=kind, variant=f"{s}-{kind}", code_hash=CODE_HASH)
    return row


def compute(jobs_n):
    done = {}
    if os.path.exists(CACHE):
        for line in open(CACHE):
            r = json.loads(line)
            if r["code_hash"] == CODE_HASH:
                done[(r["q"], r["lam"], r["n"], r["variant"])] = r
    # n = 256 only for lam in {1, 8} and the two rigorous variants (P1-max, CR-min): compute budget
    jobs = [(q, lam, n, s, k) for q in QS for lam in LAMS for n in (8, 16, 32, 64, 128, 256) for s, k in VARIANTS
            if (n < 256 or (lam in (1.0, 8.0) and k in ("max", "min"))) and (q, lam, n, f"{s}-{k}") not in done]
    jobs.sort(key=lambda j: (-j[2], j[0], j[1]))
    t0 = time.time()
    with ProcessPoolExecutor(jobs_n) as ex, open(CACHE, "a") as fh:
        for row in ex.map(job, jobs):
            fh.write(json.dumps(row) + "\n"); fh.flush()
            done[(row["q"], row["lam"], row["n"], row["variant"])] = row
            print(f"[{time.time()-t0:6.0f}s] q={row['q']} lam={row['lam']} n={row['n']} {row['variant']:8s} it={row['iterations']} "
                  f"conv={row['converged']} J={row['J']:.8f} D={row['best_D']:.8f} ({row['time']:.0f}s)", flush=True)
    rows = sorted(done.values(), key=lambda r: (r["q"], r["lam"], r["n"], r["variant"]))
    write_rows("exp07_double_phase_elm", rows, __file__, p=1.5, alpha=1.0, r=3.0)


def plot():
    plt = style()
    rows = read_rows("exp07_double_phase_elm")
    get = lambda q, lam, v: {x["n"]: x for x in rows if x["q"] == q and x["lam"] == lam and x["variant"] == v}
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7), sharey=True)
    for ax, q in zip(axs, QS):
        for j, lam in enumerate(LAMS):
            P, Cm = get(q, lam, "P1-max"), get(q, lam, "CR-min")
            ns = sorted(set(P) & set(Cm))
            ax.loglog([Cm[n]["ndof"] for n in ns], [(P[n]["J"] - Cm[n]["best_D"]) / max(1.0, abs(P[n]["J"])) for n in ns],
                      color=C[j], marker=MK[j], label=rf"$\lambda={lam:g}$")
        x = np.array([1e2, 4e5]); ax.loglog(x, 3 / x, "k:", lw=0.8)
        ax.set_title(rf"$q={q:g}$" + {2.4: " (no gap: $q\\leq p+\\alpha$)", 2.8: " (grey zone)", 4.0: " ($q>d+\\alpha$)"}[q])
        ax.set_xlabel("CR free dofs")
    axs[0].set_ylabel(r"relative width $(J_h^{P1,\max}-D^{CR,\min})/J_h$")
    axs[0].legend(fontsize=6.5)
    fig.suptitle(r"Double phase $\frac{t^{3/2}}{3/2}+a(x)\frac{t^q}{q}$, ELM weight, $\alpha=1$: guaranteed bracket (dotted: ndof$^{-1}$)", fontsize=9)
    fig.tight_layout(); savefig(fig, "F13_elm_width.pdf")
    # diagnostics Einf(P1) and N2(CR)
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6), sharex=True)
    for ax, q in zip(axs, QS):
        for j, lam in enumerate(LAMS):
            for v, ls in (("P1-bary", "-"), ("CR-min", "--")):
                R = sorted(get(q, lam, v).values(), key=lambda x: x["n"])
                ax.semilogy([np.log(x["n"] / 2) for x in R], [x["Einf"] for x in R], ls, color=C[j], marker=MK[j],
                            mfc="white" if v.startswith("CR") else C[j], label=rf"$\lambda={lam:g}$" if v == "P1-bary" else None)
        ax.set_title(rf"$q={q:g}$: $\max_T|\eta_h|$ (solid P1, dashed CR)"); ax.set_xlabel(r"$\log(1/h)$")
    axs[0].legend(fontsize=6.5); fig.tight_layout(); savefig(fig, "F14_elm_diagnostics.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--plot", action="store_true"); ap.add_argument("--jobs", type=int, default=1)
    a = ap.parse_args()
    plot() if a.plot else compute(a.jobs)
