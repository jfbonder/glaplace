"""EXPERIMENT 1 (NA requests 1, 2, 4a-b; PDE request E4(i) + section 0): P1 vs CR with the guaranteed bracket on the
Zhikov checkerboards, lambda sweeps, criss-cross meshes n = 8 ... 256 (512 for selected lambda).

Setting: Omega = (-1,1)^2, f = 0, data lam*u0 (trace lam*x2), crisscross_mesh(n) (resolves |x1| = |x2|: G constant on every
triangle), r = 3, nu_0 = eta_1 = 0, stop when scaled primal residual < 1e-10 AND discrete duality gap
J_h(u_n) - D_h(sigma_n) < 1e-9 max(1,|J_h|), max_iter 50000.
  * P1 energy J_h(u_h^P1)  >= m_H >= m_W                   (NA Prop. P1; any P1 iterate)
  * CR lower bound D(sigma_n) <= m_W at every iteration    (Theorem D(iii)); best_D = max_n D(sigma_n)
  * width Delta_h = J_h(P1) - best_D(CR)  >= m_H - m_W (gap), rigorous
Recorded per run: scheme, n, ndof (free dofs), iterations, residual, gap, max interior normal jump of sigma_n (CR, every
check), energies, best bounds, PDE diagnostics N2 = ||grad_h u_h||, E2 = ||eta_h||, Einf = max|eta_h|, E0, run time.

Runs are cached one per line in results_v4/exp01_cache.jsonl (re-running skips finished runs); the final table is
results_v4/exp01_lavrentiev.csv/.json.

    python experiments_v4/exp01_lavrentiev_sweep.py [--jobs 2] [--maxn 256] [--big]      (--big adds the n = 512 runs)
"""
import argparse, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import run_lav, write_rows, RESULTS, CODE_HASH

LAMS = {"px": [0.25, 0.5, 0.75, 0.99, 1.0, 1.02, 1.04, 1.06, 1.08, 1.5, 2.0, 3.0, 4.0],
        "tp": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0]}
BIG = {"px": [0.5, 3.0], "tp": [0.2, 3.0]}            # n = 512
CACHE = os.path.join(RESULTS, "exp01_cache.jsonl")


def key(p, lam, n, s):
    return f"{p}|{lam:g}|{n}|{s}"


def job(args):
    p, lam, n, s = args
    row = run_lav(p, lam, n, s, r=3.0, check_every=1 if n <= 32 else 5)
    row["code_hash"] = CODE_HASH
    return row


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--maxn", type=int, default=256)
    ap.add_argument("--big", action="store_true")
    ap.add_argument("--only", default=None, help="px or tp")
    a = ap.parse_args()
    done = {}
    if os.path.exists(CACHE):
        for line in open(CACHE):
            r = json.loads(line)
            if r.get("code_hash") == CODE_HASH:
                done[key(r["problem"], r["lam"], r["n"], r["scheme"])] = r
    ns = [n for n in (8, 16, 32, 64, 128, 256) if n <= a.maxn]
    jobs = []
    for p, lams in LAMS.items():
        if a.only and p != a.only:
            continue
        for lam in lams:
            for n in ns + ([512] if a.big and lam in BIG[p] else []):
                for s in ("P1", "CR"):
                    if key(p, lam, n, s) not in done:
                        jobs.append((p, lam, n, s))
    jobs.sort(key=lambda j: -(j[2] ** 2) * (3 if j[3] == "CR" else 1))      # largest first
    print(f"{len(done)} cached, {len(jobs)} to run", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.jobs) as ex, open(CACHE, "a") as fh:
        for row in ex.map(job, jobs):
            fh.write(json.dumps(row) + "\n"); fh.flush()
            done[key(row["problem"], row["lam"], row["n"], row["scheme"])] = row
            print(f"[{time.time() - t0:7.0f}s] {row['problem']} lam={row['lam']:<5g} n={row['n']:4d} {row['scheme']} "
                  f"it={row['iterations']:6d} conv={row['converged']} J={row['J']:.10f} D={row['D']:.10f} "
                  f"gap={row['gap']:.1e} jump={row['max_jump']:.1e} ({row['time']:.0f}s)", flush=True)
    rows = sorted(done.values(), key=lambda r: (r["problem"], r["lam"], r["n"], r["scheme"]))
    write_rows("exp01_lavrentiev", rows, __file__, setting="crisscross, f=0, data lam*u0, r=3, tol_res=1e-10, tol_gap=1e-9 rel")
