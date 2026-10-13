"""Print the markdown tables used in reports/PROG_REPORT.md from the result files (no computation).

    python experiments_v4/report_tables.py > /tmp/tables.md
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_rows, RESULTS


def meta(name):
    return json.load(open(os.path.join(RESULTS, name + ".json")))["meta"]


def hdr(name):
    m = meta(name)
    print(f"\n### {name}  (code_hash {m['code_hash']}, {m['date']})\n")


def t_exp02():
    hdr("exp02_r_study")
    R = read_rows("exp02_r_study")
    rs = sorted(set(x["r"] for x in R if x["n"] != 128))
    print("| case | scheme | n | " + " | ".join(f"r={r:g}" for r in rs) + " |")
    print("|---|---|---|" + "---|" * len(rs))
    for p, lam in (("px", 0.5), ("px", 3.0), ("tp", 3.0)):
        for s in ("CR", "P1"):
            for n in (32, 64):
                d = {x["r"]: x for x in R if x["problem"] == p and x["lam"] == lam and x["scheme"] == s and x["n"] == n}
                print(f"| {p} λ={lam:g} | {s} | {n} | " + " | ".join(
                    f"{d[r]['iterations']}{'' if d[r]['converged'] else '*'} ({d[r]['time']:.1f}s)" for r in rs) + " |")
    print("\nn = 128 at the best r of n = 64:")
    for x in R:
        if x["n"] == 128:
            print(f"- {x['problem']} λ={x['lam']:g} {x['scheme']}: r={x['r']:g}, {x['iterations']} it, {x['time']:.1f}s, conv={x['converged']}")
    H = json.load(open(os.path.join(RESULTS, "exp02_histories.json")))["data"]
    for r, d in H["lyapunov"].items():
        E = np.array(d["E"]); dE = np.diff(E)
        print(f"- Lyapunov r={r}: ref run {d['ref_iterations']} it (res {d['ref_residual']:.1e}); max increase "
              f"{max(dE.max(), 0):.1e} (relative to E_1: {max(dE.max(), 0) / E[0]:.1e}); E_n/E_1 at n=100: {E[min(99, len(E)-1)] / E[0]:.2e}")
    for r, h in H["history"].items():
        jm = np.nanmax(h["jump"])
        print(f"- history r={r}: {len(h['it'])} it, max normal jump over all iterations {jm:.1e}, "
              f"D(sigma_n) range [{min(h['D']):.6f}, {max(h['D']):.6f}], final gap {h['gap'][-1]:.1e}")


def t_exp03():
    hdr("exp03_bos_quadrature")
    R = read_rows("exp03_bos_quadrature")
    ns = sorted(set(x["n"] for x in R))
    print("| variant | scheme | " + " | ".join(f"n={n}" for n in ns) + " |")
    print("|---|---|" + "---|" * len(ns))
    for v in ("exact", "pmax", "pmin"):
        for s in ("P1", "CR"):
            d = {x["n"]: x for x in R if x["variant"] == v and x["scheme"] == s}
            print(f"| {v} | {s} | " + " | ".join(f"{d[n]['J']:.4f}" for n in ns) + " |")


def t_exp04():
    hdr("exp04_exponential")
    R = read_rows("exp04_exponential")
    print("| case | λ | n | J_h(P1) | D(σ) (CR) | width | P1 it | CR it | conv |")
    print("|---|---|---|---|---|---|---|---|---|")
    keys = sorted(set((x["problem"], x["lam"], x["n"]) for x in R))
    for p, lam, n in keys:
        P = [x for x in R if (x["problem"], x["lam"], x["n"], x["scheme"]) == (p, lam, n, "P1")][0]
        Q = [x for x in R if (x["problem"], x["lam"], x["n"], x["scheme"]) == (p, lam, n, "CR")][0]
        print(f"| {p} | {lam:g} | {n} | {P['J']:.8f} | {Q['best_D']:.8f} | {P['J'] - Q['best_D']:.3e} | {P['iterations']} | "
              f"{Q['iterations']} | {P['converged'] and Q['converged']} |")


def fit(rows, key, k=3):
    rows = sorted(rows, key=lambda x: -x["h"])[-k:]
    return float(np.polyfit(np.log([x["h"] for x in rows]), np.log([x[key] for x in rows]), 1)[0])


def t_exp05():
    hdr("exp05_disk_radial")
    R = read_rows("exp05_disk_radial")
    cases = []
    for x in R:
        if x["case"] not in cases:
            cases.append(x["case"])
    print("| case | scheme | order L2(u) | order L2(∇u) | order L^1.5(∇u) | finest: L2 rel | ∇ L2 rel | modular | ‖η_h−η‖ | it (coarse→fine, r=1) | conv |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for c in cases:
        for s in ("P1", "CR"):
            S = sorted([x for x in R if x["case"] == c and x["scheme"] == s], key=lambda x: -x["h"])
            f = S[-1]
            print(f"| {c} | {s} | {fit(S,'L2_rel'):.2f} | {fit(S,'grad_L2_rel'):.2f} | {fit(S,'grad_L1.5_rel'):.2f} | {f['L2_rel']:.2e} | "
                  f"{f['grad_L2_rel']:.2e} | {f['modular']:.2e} | {f['eta_L2_err']:.2e} | {'→'.join(str(x['iterations']) for x in S)} | "
                  f"{all(x['converged'] for x in S)} |")
    print("\nDiagnostics (section 0) on the finest 4 meshes:\n")
    print("| case | scheme | h | N2 | E2 | Einf | E0 | iters |\n|---|---|---|---|---|---|---|---|")
    for c in cases:
        for s in ("P1", "CR"):
            S = sorted([x for x in R if x["case"] == c and x["scheme"] == s], key=lambda x: -x["h"])[-4:]
            for x in S:
                print(f"| {c} | {s} | {x['h']:.4f} | {x['N2']:.5f} | {x['E2']:.5f} | {x['Einf']:.4f} | {x['E0']:.4f} | {x['iterations']} |")
    H = json.load(open(os.path.join(RESULTS, "exp05_histories.json")))["data"]
    for h in H:
        mod = np.array(h["mod"]); eta = np.array(h["eta"])
        k = [10, 100, 1000]
        print(f"- history {h['case']} (h={h['h']}, ref {h['ref_iterations']} it): modular at n=" +
              ", ".join(f"{i}: {mod[i-1]:.1e}" for i in k if i <= len(mod)) + "; ‖η_n−η_h‖ at n=" +
              ", ".join(f"{i}: {eta[i-1]:.1e}" for i in k if i <= len(eta)))


def t_exp06():
    hdr("exp06_rough_load")
    R = read_rows("exp06_rough_load")
    print("| β | scheme | load | h | ndof | it | J(u_h) − J(u) | N2² | ΔN2² | E2² (exact) | Einf | E0 | ‖∇(u−u_h)‖_{1.5} rel | conv |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for beta in (1.5, 1.4, 1.6):
        for s in ("P1", "CR"):
            for load in ("accurate", "bary"):
                S = sorted([x for x in R if x["beta"] == beta and x["scheme"] == s and x["load"] == load], key=lambda x: -x["h"])
                prev = None
                for x in S:
                    dn = "" if prev is None else f"{x['N2sq'] - prev:.1f}"
                    prev = x["N2sq"]
                    print(f"| {beta} | {s} | {load} | {x['h']:.4f} | {x['ndof']} | {x['iterations']} | {x['J_minus_exact']:.4f} | "
                          f"{x['N2sq']:.1f} | {dn} | {x['E2']**2:.3f} ({x['E2sq_exact']:.3f}) | {x['Einf']:.2f} | {x['E0']:.1f} | "
                          f"{x['grad_L15_rel']:.3e} | {x['converged']} |")


def t_exp07():
    hdr("exp07_double_phase_elm")
    R = read_rows("exp07_double_phase_elm")
    print("| q | λ | n | P1-bary | P1-max (upper) | CR-bary | CR-min D (lower) | width rel | ratio | Einf P1-bary | N2 CR-min | conv |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for q in (2.4, 2.8, 4.0):
        for lam in (0.5, 1.0, 2.0, 4.0, 8.0):
            d = {}
            for x in R:
                if x["q"] == q and x["lam"] == lam:
                    d[(x["n"], x["variant"])] = x
            ns = sorted(set(n for n, v in d))
            prevW = None
            for n in ns:
                try:
                    a, b, c, e = (d[(n, v)] for v in ("P1-bary", "P1-max", "CR-bary", "CR-min"))
                except KeyError:
                    continue
                W = b["J"] - e["best_D"]
                rat = "" if prevW is None else f"{prevW / W:.2f}"
                prevW = W
                if n >= 32:
                    print(f"| {q} | {lam:g} | {n} | {a['J']:.6f} | {b['J']:.6f} | {c['J']:.6f} | {e['best_D']:.6f} | "
                          f"{W / max(1, abs(b['J'])):.2e} | {rat} | {a['Einf']:.2f} | {e['N2']:.3f} | "
                          f"{all(x['converged'] for x in (a, b, c, e))} |")


def t_exp08():
    for name in ("exp08_orders", "exp08_p_to_1"):
        hdr(name)
        R = read_rows(name)
        print("| case | order L2(u) | order L2(∇u) | order max | finest L2 | finest ∇L2 | it coarsest/finest | conv (all) |")
        print("|---|---|---|---|---|---|---|---|")
        cases = []
        for x in R:
            if x["case"] not in cases:
                cases.append(x["case"])
        for c in cases:
            S = sorted([x for x in R if x["case"] == c], key=lambda x: -x["h"])
            print(f"| {c} | {fit(S,'L2_rel'):.2f} | {fit(S,'gradL2_rel'):.2f} | {fit(S,'sup_nodes'):.2f} | {S[-1]['L2_rel']:.2e} | "
                  f"{S[-1]['gradL2_rel']:.2e} | {S[0]['iterations']}/{S[-1]['iterations']} | "
                  f"{''.join('y' if x['converged'] else 'n' for x in S)} |")
    hdr("exp08_manufactured")
    for x in read_rows("exp08_manufactured"):
        print(f"- {x['part']} n={x['n']} r={x['r']:g}: {x['iterations']} it, conv={x['converged']}, L2 rel {x['L2_rel']:.4e}, ∇L2 rel {x['gradL2_rel']:.4e}")
    hdr("exp08_deadcore_r")
    for x in read_rows("exp08_deadcore_r"):
        print(f"- r={x['r']:g} stop={x['stop']}: {x['iterations']} it, conv={x['converged']}, res {x['residual']:.1e}, dres "
              f"{x['dual_residual']:.1e}, gap {x['gap']:.1e}, L2 rel {x['L2_rel']:.4e}, core max|∇u_h| {x['core_max_grad']:.1e}")
    hdr("exp08_twophase_holder")
    for x in read_rows("exp08_twophase_holder"):
        print(f"- {x['a']} h={x['h']}: {x['iterations']} it, J={x['J']:.4f}, max u={x['max_u']:.4f}, max|∇u|={x['max_grad']:.3f}, "
              f"‖∇u‖={x['gradL2']:.4f}, ‖η‖={x['E2']:.4f}, max|η|={x['Einf']:.3f}")
    L = json.load(open(os.path.join(RESULTS, "exp08_lyapunov.json")))["data"]
    for k, v in L.items():
        E = np.array(v["E"])
        print(f"- Lyapunov {k}: strictly decreasing while E_n/E_1 > 1e-20: {v['monotone_until_1e_20']}; E_N/E_1 = {E[-1]/E[0]:.1e}")


if __name__ == "__main__":
    for f in (t_exp02, t_exp03, t_exp04, t_exp05, t_exp06, t_exp07, t_exp08):
        try:
            f()
        except FileNotFoundError as e:
            print(f"\n(missing: {e.filename})")
