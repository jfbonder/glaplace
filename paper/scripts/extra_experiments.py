"""Extra experiments for the article (all with glaplace, PYTHONPATH=/home/claude/glaplace).
  sc   : true L2 error vs nodal-difference L2 norm ||I_h u - u_h|| (supercloseness) for radial torsion cases
  pto1 : p-Laplacian torsion for p = 1.2, 1.1 (iterations, convergence flag, true errors), max_iter 6000
  rman : manufactured x-dependent G: orders (n = 10..80, r = 100) and iterations vs r (n = 40)
  rdead: dead-core e^t - 1, c = 8, h = 0.05: iterations vs r
Usage: python extra_experiments.py sc|pto1|rman|rdead  > out.txt
"""
import sys, time, numpy as np
import glaplace as gl

ring = lambda x, y: np.where(np.hypot(x, y) > 0.5, 1.0, 0.0)

def sc():
    cases = {"p2": (gl.PowerLaw(2.0), 1.0, 1.0), "p4": (gl.PowerLaw(4.0), 1.0, 1.0), "pq24": (gl.pq(2, 4), 1.0, 1.0),
             "exp": (gl.Exp(), 1.0, 1.0)}
    for name, (G, c, r) in cases.items():
        u, gu, _ = gl.exact.radial_torsion(G, c)
        for h in (0.1, 0.05, 0.025, 0.0125):
            m = gl.disk_mesh(h); F = gl.P1Space(m)
            res = gl.GLaplaceSolver(m, G, c, 0.0, r=r).solve(tol=1e-10, max_iter=5000, stop="residual")
            ue = u(*m.vertices.T)
            true = F.l2_error(res.u, u)
            nodal = F.norm_L2(res.u - ue) / F.norm_L2(ue)
            print(f"{name} h={h} nodes={m.n_nodes} it={res.iterations} conv={res.converged} true_L2_rel={true:.4e} nodal_L2_rel={nodal:.4e}", flush=True)

def pto1():
    for p in (1.2, 1.1):
        G = gl.PowerLaw(p); u, gu, _ = gl.exact.radial_torsion(G, 1.0)
        for h in (0.2, 0.1, 0.05, 0.025):
            m = gl.disk_mesh(h); F = gl.P1Space(m); t0 = time.perf_counter()
            res = gl.GLaplaceSolver(m, G, 1.0, 0.0, r=1.0).solve(tol=1e-9, max_iter=6000, stop="residual")
            print(f"p={p} h={h} nodes={m.n_nodes} it={res.iterations} conv={res.converged} "
                  f"residual={res.history['residual'][-1]:.2e} L2_rel={F.l2_error(res.u, u):.4e} "
                  f"gradL2_rel={F.grad_error_Lp(res.u, gu, 2.0):.4e} gradLp_rel={F.grad_error_Lp(res.u, gu, p):.4e} "
                  f"umax_exact={u(0.,0.):.4e} time={time.perf_counter()-t0:.0f}", flush=True)

def rman():
    f, u, gu = gl.exact.manufactured("t**2/2 + (1+x)*t**4/4", "exp(x+y)")
    G = gl.pq(2, 4, lambda x, y: 1 + x)
    for n in (10, 20, 40, 80):
        m = gl.rectangle_mesh(n, n); F = gl.P1Space(m)
        res = gl.GLaplaceSolver(m, G, f, u, r=100.0, f_degree=2).solve(tol=1e-10, max_iter=5000, stop="residual")
        print(f"order n={n} nodes={m.n_nodes} it={res.iterations} conv={res.converged} L2_rel={F.l2_error(res.u, u):.4e} "
              f"gradL2_rel={F.grad_error_Lp(res.u, gu, 2.0):.4e}", flush=True)
    m = gl.rectangle_mesh(40, 40); F = gl.P1Space(m)
    for r in (1, 10, 30, 100, 300):
        res = gl.GLaplaceSolver(m, G, f, u, r=float(r), f_degree=2).solve(tol=1e-10, max_iter=5000, stop="residual")
        print(f"r={r} it={res.iterations} conv={res.converged} L2_rel={F.l2_error(res.u, u):.4e}", flush=True)

def rdead():
    G = gl.ExpLiteral(); u, gu, _ = gl.exact.radial_torsion(G, 8.0)
    m = gl.disk_mesh(0.05); F = gl.P1Space(m)
    for r in (1.0, 3.0, 10.0, 30.0, 100.0):
        res = gl.GLaplaceSolver(m, G, 8.0, 0.0, r=r).solve(tol=1e-8, max_iter=6000, stop="residual")
        gx, gy = F.gradient(res.u); core = np.hypot(*F.bary.T) < 0.2
        print(f"r={r} it={res.iterations} conv={res.converged} residual={res.history['residual'][-1]:.2e} "
              f"L2_rel={F.l2_error(res.u, u):.4e} max_grad_core={np.hypot(gx, gy)[core].max():.2e}", flush=True)


def pto1_hist():
    """L2 error after 3000 and after 6000 iterations (p = 1.2, 1.1; h = 0.05, 0.025)."""
    for p in (1.2, 1.1):
        G = gl.PowerLaw(p); u, gu, _ = gl.exact.radial_torsion(G, 1.0)
        for h in (0.05, 0.025):
            m = gl.disk_mesh(h); F = gl.P1Space(m); snap = {}
            def cb(it, uu, hist):
                if it in (1500, 3000, 4500, 6000): snap[it] = F.l2_error(uu, u)
            gl.GLaplaceSolver(m, G, 1.0, 0.0, r=1.0).solve(tol=1e-30, max_iter=6000, stop="residual", callback=cb)
            print(f"p={p} h={h} " + " ".join(f"L2@{k}={v:.6e}" for k, v in sorted(snap.items())), flush=True)

def meshinfo():
    for h in (0.2, 0.1, 0.05, 0.025, 0.0125):
        m = gl.disk_mesh(h); V, T = m.vertices, m.triangles
        e = np.concatenate([np.linalg.norm(V[T[:, i]] - V[T[:, (i + 1) % 3]], axis=1) for i in range(3)])
        print(f"h_target={h} nodes={m.n_nodes} elems={m.n_elems} h_max={e.max():.4f} h_mean={e.mean():.4f}", flush=True)

def exactcheck():
    r = np.linspace(0, 1, 2001)
    for p in (1.5, 2.0, 4.0):
        u, gu, _ = gl.exact.radial_torsion(gl.PowerLaw(p), 1.0)
        q = p / (p - 1); c = (p - 1) / p * 2 ** (-1 / (p - 1)); ue = c * (1 - r ** q)
        d = np.abs(u(r, 0 * r) - ue).max()
        print(f"p={p} max|u_num-u_closed|={d:.3e} relative={d/ue.max():.3e}", flush=True)

if __name__ == "__main__":
    {"sc": sc, "pto1": pto1, "rman": rman, "rdead": rdead, "pto1_hist": pto1_hist,
     "meshinfo": meshinfo, "exactcheck": exactcheck}[sys.argv[1]]()
