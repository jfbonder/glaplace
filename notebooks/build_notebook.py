"""Generates notebooks/glaplace_colab.ipynb (self-contained: it embeds the package source)."""
import nbformat as nbf, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
nb = nbf.v4.new_notebook()
C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n")))

md("""
# g-Laplacian solver — decomposition-coordination method for Young functions $G$, on general 2D domains

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jfbonder/glaplace/blob/main/notebooks/glaplace_colab.ipynb)

We solve, with $g=G'$ and $G(x,\\cdot)$ a Young function,
$$-\\mathrm{div}\\Big(g(x,|\\nabla u|)\\frac{\\nabla u}{|\\nabla u|}\\Big)=f\\ \\text{in }\\Omega,\\qquad u=u_D\\ \\text{on }\\partial\\Omega,$$
i.e. we minimise $\\int_\\Omega G(x,|\\nabla u|)-\\int_\\Omega fu$, with the iteration (augmented Lagrangian / ADMM, parameter $r>0$):

1. $r(-\\Delta u_n)=f+\\mathrm{div}(\\eta_n-r\\nu_{n-1})$  (linear, same matrix at every step → factorised once),
2. $\\nu_n=t\\,w/|w|$, $w=\\eta_n/r+\\nabla u_n$, where $t\\ge 0$ solves $g(x,t)+rt=r|w|$  (pointwise scalar equation),
3. $\\eta_{n+1}=\\eta_n+r(\\nabla u_n-\\nu_n)$.

For $G=t^p/p$ this is the method of Aragón–Fernández Bonder–Rubio for the $p$-Laplacian. **Run all cells.**
""")
md("## 1. Package source (self-contained)")
code("!pip -q install triangle meshio sympy 2>&1 | tail -1\n!mkdir -p glaplace")
for name in ("__init__", "mesh", "fem", "young", "nonlinear", "exact", "solver"):
    src = (root / "glaplace" / f"{name}.py").read_text()
    code(f"%%writefile glaplace/{name}.py\n{src}")
code("import sys; sys.path.insert(0, '.')\nimport numpy as np, matplotlib.pyplot as plt\nimport glaplace as gl\nfrom glaplace import *\nprint('ok')")

md("""
## 2. Young functions
`PowerLaw(p)` $=t^p/p$, `PLog(p, a)` $=t^p\\log^a(1+t)$, `Exp()` $=e^t-1-t$, `ExpLiteral()` $=e^t-1$, `Custom(G, g)`, `from_expr("t**3/3")` (sympy),
and sums with weights: `PowerLaw(2) + a * PowerLaw(4)` with `a` a number or a function `a(x, y) >= 0` (two-phase).
`check_young` tests numerically that $G(0)=0$, $g\\ge0$ and $g$ is non-decreasing.
""")
code("""
for Gf in (pq(2, 4), PLog(2, 1), PLog(3, -1), Exp(), ExpLiteral(), from_expr("t**3/3 + t**2/2")):
    check_young(Gf)
""")
md("## 3. Quick start: the $p$-$q$-Laplacian $G=t^2/2+t^4/4$ on the unit disk, compared with the exact radial solution")
code("""
G = pq(2, 4)
mesh = disk_mesh(0.04)
res = GLaplaceSolver(mesh, G, f=1.0, u_D=0.0).solve(tol=1e-9, stop='residual')
u_ex, grad_ex, _ = exact.radial_torsion(G, c=1.0)      # g(|u'|) = c r/2 solved numerically
fem = P1Space(mesh)
print(f"{res.iterations} iterations, converged={res.converged}, {res.time_iter:.1f}s")
print("true relative L2 error =", fem.l2_error(res.u, u_ex), "  gradient L2 error =", fem.grad_error_Lp(res.u, grad_ex, 2.0))
mesh.plot(res.u); plt.title('p-q-Laplacian torsion'); plt.show()
""")
md("""
`res.history` stores, per iteration, `rel_u` (relative change of $u$), `residual` ($\\|\\nabla u-\\nu\\|_{L^2}/(1+\\|\\nabla u\\|_{L^2})$) and `energy`.
**Choice of $r$:** the iteration is fastest when $r$ is of the order of $g(t)/t$ at the typical gradient size $t$ (for $G=t^p/p$ and $|\\nabla u|\\approx1$, $r\\approx1$; for larger gradients take larger $r$).
""")
code("""
h = res.history
fig, ax = plt.subplots(1, 3, figsize=(14, 3.3))
ax[0].semilogy(h['rel_u']); ax[0].set_title('rel_u'); ax[1].semilogy(h['residual']); ax[1].set_title('residual'); ax[2].plot(h['energy']); ax[2].set_title('energy'); plt.show()
""")
md("## 4. Torsion on the disk for several $G$ (true errors vs. the exact solution; order in $h$)")
code("""
cases = {'t^1.5/1.5': PowerLaw(1.5), 't^2/2+t^4/4': pq(2, 4), 't^2 log(1+t)': PLog(2, 1), 'e^t-1-t': Exp()}
print(f"{'G':14s} {'h':>7} {'it':>4} {'L2 rel err':>11} {'order':>6}")
for name, Gf in cases.items():
    u_ex, _, _ = exact.radial_torsion(Gf, 1.0); prev = None
    for h0 in (0.1, 0.05, 0.025):
        m = disk_mesh(h0); r_ = GLaplaceSolver(m, Gf, 1.0, 0.0).solve(tol=1e-9, stop='residual')
        e = P1Space(m).l2_error(r_.u, u_ex)
        print(f"{name:14s} {h0:7.3f} {r_.iterations:4d} {e:11.3e} {'' if prev is None else round(np.log(prev/e)/np.log(2), 2)}"); prev = e
""")
md("""
## 5. Two-phase problem: $G(x,t)=t^p/p+a(x)\\,t^q/q$ with $a\\ge0$ vanishing on a sector
""")
code("""
ang = lambda x, y: np.mod(np.arctan2(y, x), 2*np.pi)
a = lambda x, y: np.where((ang(x, y) > 0) & (ang(x, y) < np.pi/2), 0.0, 1.0)
G2 = pq(2, 6, a)
m = disk_mesh(0.03)
r2 = GLaplaceSolver(m, G2, 4.0, 0.0, r=3.0).solve(tol=1e-9, stop='residual')
fem = P1Space(m); gx, gy = fem.gradient(r2.u)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
m.plot(r2.u, ax=ax[0]); ax[0].set_title(f'u ({r2.iterations} it)')
c = ax[1].tripcolor(m.vertices[:, 0], m.vertices[:, 1], m.triangles, facecolors=np.hypot(gx, gy)); plt.colorbar(c, ax=ax[1]); ax[1].set_aspect('equal'); ax[1].set_title('|grad u|')
plt.show()
""")
md("""
## 6. $G(t)=e^t-1$: no $\\Delta_2$ condition and $g(0^+)=1$ → dead core
Torsion with load $c$: $g(|u'|)=cr/2$ has a solution only where $cr/2>1$; $u$ is flat for $r\\le 2/c$.
""")
code("""
Gx = ExpLiteral(); c = 8.0
u_ex, _, ur = exact.radial_torsion(Gx, c)
m = disk_mesh(0.05)
rx = GLaplaceSolver(m, Gx, c, 0.0, r=30.0).solve(tol=1e-8, max_iter=6000, stop='residual')
print(rx.iterations, rx.converged, 'true rel L2 error', P1Space(m).l2_error(rx.u, u_ex))
rr = np.linspace(0, 1, 300); R = np.hypot(*m.vertices.T)
plt.plot(rr, u_ex(rr, 0*rr), 'k', label='exact'); plt.plot(R, rx.u, '.', ms=2, label='computed'); plt.axvline(2/c, ls=':'); plt.legend(); plt.show()
""")
md("""
## 7. Other domains and your own data
```python
mesh = rectangle_mesh(100, 100)                                         # (0,1)^2
mesh = polygon_mesh(outer_xy, holes=[hole_xy], max_area=0.002)           # polygon with holes (needs `triangle`)
mesh = sdf_mesh(lambda P: np.hypot(P[:,0], P[:,1]) - 1, (-1,-1,1,1), h0=0.03)   # signed distance function
mesh = Mesh(vertices, triangles, dirichlet_node_indices)                 # your own arrays (0-based)
mesh = from_txt('vertex_coordinates.txt', 'elem_vertices.txt', 'dirichlet.txt')  # MatLab files (1-based)
mesh = from_meshio('domain.msh')                                         # Gmsh etc.
```
`f`, `u_D` and the weights in `G` can be constants or vectorised functions of `(x, y)`; $G(x,\\cdot)$ is sampled at the triangle barycentres.
""")
nb["cells"] = C
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["colab"] = {"provenance": []}
out = root / "notebooks" / "glaplace_colab.ipynb"
nbf.write(nb, out); print("wrote", out)
