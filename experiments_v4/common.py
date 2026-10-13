"""Shared helpers for the v4 experiment scripts (results with code hash, figure style, Lavrentiev test runs).

Every result file (CSV + JSON) records `code_hash` = sha256 of the glaplace/*.py sources (glaplace.admm.source_hash),
the date, the platform and the script that produced it.
"""
from __future__ import annotations
import csv, json, os, sys, time, platform, datetime
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")

from glaplace.admm import Discretisation, ADMM, source_hash          # noqa: E402
from glaplace.cr import CRSpace                                       # noqa: E402
from glaplace.lavrentiev import crisscross_mesh, G_pxlap, G_twophase, psi, u0   # noqa: E402

RESULTS = os.environ.get("GL_RESULTS", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "v4"))
FIGURES = os.environ.get("GL_FIGURES", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "v4", "figures"))
os.makedirs(RESULTS, exist_ok=True); os.makedirs(FIGURES, exist_ok=True)
CODE_HASH = source_hash()

# Okabe-Ito (colour-blind safe; validated: all checks pass, secondary encoding by markers/linestyles always used)
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
MK = ["o", "s", "^", "D", "v", "P", "X"]


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5, "legend.fontsize": 7.5,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "lines.linewidth": 1.4, "lines.markersize": 4.5,
        "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.5, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "savefig.bbox": "tight", "pdf.fonttype": 42,
        "mathtext.fontset": "cm", "axes.prop_cycle": matplotlib.cycler(color=C)})
    return plt


def savefig(fig, name):
    path = os.path.join(FIGURES, name)
    fig.savefig(path)
    print("figure ->", path)
    return path


def meta(script, **extra):
    return dict(code_hash=CODE_HASH, script=os.path.basename(script), date=datetime.datetime.now().isoformat(timespec="seconds"),
                python=platform.python_version(), numpy=np.__version__, machine=platform.machine(), **extra)


def _clean(v):
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, float) and not np.isfinite(v):
        return None
    return v


def write_rows(name, rows, script, **extra):
    """rows: list of flat dicts -> RESULTS/name.csv (with a '# code_hash=' comment line) and name.json (meta + rows)."""
    rows = [{k: _clean(v) for k, v in r.items()} for r in rows]
    m = meta(script, **extra)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(RESULTS, name + ".csv"), "w", newline="") as fh:
        fh.write(f"# code_hash={m['code_hash']} script={m['script']} date={m['date']}\n")
        w = csv.DictWriter(fh, keys); w.writeheader(); w.writerows(rows)
    with open(os.path.join(RESULTS, name + ".json"), "w") as fh:
        json.dump(dict(meta=m, rows=rows), fh, indent=1, default=_clean)
    print("results ->", os.path.join(RESULTS, name + ".csv"))


def read_rows(name):
    with open(os.path.join(RESULTS, name + ".json")) as fh:
        return json.load(fh)["rows"]


def write_json(name, obj, script, **extra):
    with open(os.path.join(RESULTS, name + ".json"), "w") as fh:
        json.dump(dict(meta=meta(script, **extra), data=obj), fh, indent=1, default=_clean)
    print("results ->", os.path.join(RESULTS, name + ".json"))


# ---------------------------------------------------------------------------------------------------------------
PROBLEMS = {"px": G_pxlap, "tp": G_twophase}
PROBLEM_LABEL = {"px": r"$p(x)$ checkerboard $(3/2,3)$", "tp": r"two-phase $\frac{t^{3/2}}{3/2}+\mathbf{1}_{|x_1|<|x_2|}\frac{t^3}{3}$"}


def lav_disc(n, lam, scheme, mesh=None):
    """Criss-cross mesh of (-1,1)^2, data lam*u0: P1 nodal values / CR exact edge means (trace lam*x2 affine per edge)."""
    m = crisscross_mesh(n) if mesh is None else mesh
    if scheme == "P1":
        gD = psi(lam)(*m.vertices[m.dirichlet].T)
    else:
        fem = CRSpace(m); gD = fem.edge_values(psi(lam), fem.bdry_edges)
    return m, Discretisation(m, scheme, gD)


def run_lav(problem, lam, n, scheme, r=3.0, G=None, tol_res=1e-10, tol_gap=1e-9, max_iter=50000, check_every=None,
            jumps=True, record=False, **kw):
    """One Lavrentiev run; returns a flat row (plus the raw output under key '_out' if record)."""
    m, d = lav_disc(n, lam, scheme)
    GG = PROBLEMS[problem]() if G is None else G(d)
    A = ADMM(d, GG, r=r)
    if check_every is None:
        check_every = 1 if n <= 32 else 5
    t0 = time.perf_counter()
    o = A.solve(tol_res=tol_res, tol_gap=tol_gap, max_iter=max_iter, check_every=check_every, jumps=jumps and scheme == "CR",
                record=record or jumps, **kw)
    row = dict(problem=problem, lam=lam, n=n, h=2.0 / n, scheme=scheme, ndof=len(d.free), Nt=d.Nt, r=r,
               iterations=o["iterations"], converged=o["converged"], residual=o["residual"], dual_residual=o["dual_residual"],
               gap=o["gap"], rel_gap=o["gap"] / max(1.0, abs(o["J"])), max_jump=o["max_jump"], J=o["J"], D=o["D"],
               best_D=o["best_D"], best_J=o["best_J"], N2=o["N2"], E2=o["E2"], Einf=o["Einf"], Ginf=o["Ginf"], E0=o["E0"],
               time=time.perf_counter() - t0, time_setup=A.time_setup)
    if record:
        row["_out"] = o
    return row


def fmt(x, d=10):
    return f"{x:.{d}f}" if isinstance(x, (float, np.floating)) and np.isfinite(x) else str(x)
