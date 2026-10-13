"""Generates notebooks/glaplace_cr_lavrentiev.ipynb (self-contained; embeds the package source, including glaplace/cr.py and
glaplace/lavrentiev.py).  Does not touch glaplace_colab.ipynb / build_notebook.py."""
import nbformat as nbf, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
nb = nbf.v4.new_notebook()
C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n")))

md(r"""
# Lavrentiev gap: conforming P1 vs. nonconforming Crouzeix–Raviart (CR), with the decomposition–coordination method

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jfbonder/glaplace/blob/main/notebooks/glaplace_cr_lavrentiev.ipynb)

For integrands with a **Lavrentiev gap** every conforming finite element space (functions in $W^{1,\infty}$) converges to the infimum over
*Lipschitz/smooth* functions, $m_{\mathcal H}$, which is larger than the true minimum $m_{\mathcal W}$ over the energy space
(Carstensen–Ortner 2010). The **Crouzeix–Raviart** element (midpoint continuity, nonconforming) converges to the $\mathcal W$-minimiser
under suitable assumptions (Balci–Ortner–Storn 2022, one-point quadrature $\varphi(x_T,\cdot)$).

**Implementation.** The broken gradient of a CR function is constant on each triangle, exactly like the P1 gradient. Hence $\nu,\eta$ (piecewise constants)
and the pointwise nonlinear step are *identical* to the P1 code; only the space of $u$ changes (`glaplace/cr.py`, the P1 code is untouched).
At each mesh the iteration is a finite-dimensional ADMM, so it converges for every $r>0$ (Theorem C of the paper).
What happens as $h\to0$ is **not** covered by our theory: it is the Lavrentiev question.

**Test problems** (Zhikov checkerboard; $\Omega=(-1,1)^2$, $f=0$, boundary data $\lambda u_0$, $u_0=1$ for $|x_1|<x_2$, $-1$ for $|x_1|<-x_2$, $x_2/|x_1|$ otherwise):

* **p(x)-Laplacian with discontinuous $p$:** $G=t^{p(x)}/p(x)$, $p=3$ if $|x_1|<|x_2|$, $p=3/2$ otherwise;
* **two-phase with discontinuous $a$:** $G=t^{3/2}/(3/2)+a(x)\,t^3/3$, $a=1$ if $|x_1|<|x_2|$, $a=0$ otherwise.

We compare the **energies** $F_h(u_h)=\sum_T|T|\,G(x_T,\nabla_h u_h|_T)$ (the same formula for both schemes). No exact solution is available except for $\lambda=1$, $p(x)$ problem
($u=x_2$, $F=2$). **Run all cells.**
""")
md("## 1. Package source (self-contained)")
code("!pip -q install triangle meshio sympy 2>&1 | tail -1\n!mkdir -p glaplace")
for name in ("__init__", "mesh", "fem", "young", "nonlinear", "exact", "solver", "cr", "lavrentiev"):
    src = (root / "glaplace" / f"{name}.py").read_text()
    code(f"%%writefile glaplace/{name}.py\n{src}")
code("""
import sys; sys.path.insert(0, '.')
import numpy as np, matplotlib.pyplot as plt, matplotlib.tri as mtri
import glaplace as gl
from glaplace.cr import GLaplaceCRSolver, CRSpace
from glaplace.solver import GLaplaceSolver
from glaplace.lavrentiev import (crisscross_mesh, psi, G_pxlap, G_twophase, compare, energy_exact_pxlap_lambda1, p_zhikov, a_checker)
print('ok')
""")

md(r"""
## 2. Sanity check: $\lambda=1$, $p(x)$-Laplacian — the exact minimiser is $u=x_2$, $F=2/p_-+2/p_+=2$
Both schemes must reproduce it (to rounding) on every mesh: the flux $(0,1)$ is divergence-free in both regions, and $x_2$ is in both spaces.
""")
code("""
print('exact F =', energy_exact_pxlap_lambda1())
rows = compare(G_pxlap(), 1.0, ns=(8, 16, 32), r=1.0)
""")

md(r"""
## 3. $p(x)$-Laplacian with discontinuous $p$ ($p_-=3/2$, $p_+=3$), $\lambda=3$ (large: gap expected)
""")
code("""
rows_px3 = compare(G_pxlap(), 3.0, ns=(8, 16, 32, 64), r=3.0, max_iter=8000)
""")
code("""
def plot_energies(rows, title):
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for s, c in (('P1', 'C3'), ('CR', 'C0')):
        R = [x for x in rows if x['scheme'] == s]
        h = np.array([1/np.sqrt(x['ndof']) for x in R]); E = np.array([x['energy'] for x in R])
        ax[0].plot(h, E, 'o-', c=c, label=s)
    ax[0].set_xlabel('ndof$^{-1/2}$'); ax[0].set_ylabel('$F_h(u_h)$'); ax[0].legend(); ax[0].set_title(title); ax[0].set_xlim(0)
    P = [x['energy'] for x in rows if x['scheme'] == 'P1']; Q = [x['energy'] for x in rows if x['scheme'] == 'CR']
    ns = [x['n'] for x in rows if x['scheme'] == 'P1']
    ax[1].plot(ns, np.array(P) - np.array(Q), 'ks-'); ax[1].set_xscale('log', base=2); ax[1].set_xlabel('n (squares per side)')
    ax[1].set_ylabel('$F_h(P1)-F_h(CR)$'); ax[1].set_title('gap between the two discrete minima'); ax[1].set_ylim(0)
    plt.tight_layout(); plt.show()
plot_energies(rows_px3, 'p(x), $\\\\lambda=3$')
""")
md(r"""
The P1 energies decrease and the CR energies increase as $h\to0$, from opposite sides. The difference $F_h(\mathrm{P1})-F_h(\mathrm{CR})$ decreases only slowly (about 2.9, 2.2, 1.8, 1.6 on the four meshes), so these meshes cannot show whether it tends to a positive limit. The theory of the Lavrentiev gap predicts that it does: the P1 minima converge to $m_{\mathcal H}$ and the CR minima to $m_{\mathcal W}<m_{\mathcal H}$ (Carstensen–Ortner; Balci–Ortner–Storn). The convergence of each sequence is slow (Balci–Ortner–Storn report $\mathrm{ndof}^{-1/2}$ for the energies), so refining much further, or Richardson extrapolation in $\mathrm{ndof}^{-1/2}$, would be needed to see the limits.
""")
code("""
# Solutions and gradient moduli on the finest mesh
mesh = crisscross_mesh(64); G = G_pxlap(); lam = 3.0
sP = GLaplaceSolver(mesh, G, 0.0, psi(lam), r=3.0); rP = sP.solve(tol=1e-8, stop='residual', max_iter=8000)
sC = GLaplaceCRSolver(mesh, G, 0.0, psi(lam), r=3.0); rC = sC.solve(tol=1e-8, stop='residual', max_iter=8000)
tri = mtri.Triangulation(mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.triangles)
def gradmod(S, res): gx, gy = S.fem.gradient(res.u); return np.hypot(gx, gy)
fig, ax = plt.subplots(2, 3, figsize=(14, 8.4))
uP = rP.u; uC = sC.fem.nodal_average(rC.u)
for j, (u, nm) in enumerate(((uP, 'P1'), (uC, 'CR (vertex average)'))):
    c = ax[j, 0].tripcolor(tri, u, shading='gouraud'); plt.colorbar(c, ax=ax[j, 0]); ax[j, 0].set_title(f'u, {nm}')
for j, (S, res, nm) in enumerate(((sP, rP, 'P1'), (sC, rC, 'CR'))):
    g = gradmod(S, res)
    c = ax[j, 1].tripcolor(tri, facecolors=np.log10(g + 1e-12), vmin=-2, vmax=1.5); plt.colorbar(c, ax=ax[j, 1]); ax[j, 1].set_title(f'log10 |grad u|, {nm}')
c = ax[0, 2].tripcolor(tri, facecolors=p_zhikov()(sP.xb, sP.yb)); plt.colorbar(c, ax=ax[0, 2]); ax[0, 2].set_title('p(x)')
ax[1, 2].axis('off')
for a in ax.ravel(): a.set_aspect('equal')
plt.tight_layout(); plt.show()
print('F_h: P1 =', sP.energy(rP.u), '  CR =', sC.energy(rC.u))
""")

md(r"""
## 4. Two-phase with discontinuous $a(x)$: $G=t^{3/2}/(3/2)+a(x)t^3/3$, $\lambda=1$ and $\lambda=3$
""")
code("""
rows_tp1 = compare(G_twophase(), 1.0, ns=(8, 16, 32, 64), r=3.0, max_iter=8000)
plot_energies(rows_tp1, 'two-phase, $\\\\lambda=1$')
""")
code("""
rows_tp3 = compare(G_twophase(), 3.0, ns=(8, 16, 32, 64), r=3.0, max_iter=8000)
plot_energies(rows_tp3, 'two-phase, $\\\\lambda=3$')
""")

md(r"""
## 5. Dependence on the scaling $\lambda$ ($p(x)$ problem, $n=32$)
For small $\lambda$ the two schemes should agree (no gap); the difference should grow for large $\lambda$. This is **evidence at fixed $h$, not a proof**:
at a fixed mesh the two discrete minima always differ, so the statement concerns the limit $h\to0$ (compare the previous sections).
""")
code("""
lams = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
diff = []
for lam in lams:
    r = compare(G_pxlap(), lam, ns=(32,), r=max(1.0, lam), max_iter=10000, verbose=False)
    eP, eC = [x['energy'] for x in r]; diff.append(eP - eC)
    print(f'lambda={lam:4.2f}: F_h(P1)={eP:.5f}  F_h(CR)={eC:.5f}  difference={eP-eC:.5f}  relative={(eP-eC)/eC:.2e}')
plt.figure(figsize=(5, 3.5)); plt.plot(lams, diff, 'o-'); plt.xlabel('$\\\\lambda$'); plt.ylabel('$F_h(P1)-F_h(CR)$, n=32'); plt.show()
""")

md(r"""
## 6. Remarks / caveats

* The comparison is about the **energies of the discrete minimisers**. Theorem C of the paper only gives convergence of the iteration at each fixed mesh; it says nothing about $h\to0$.
* That CR converges to the $\mathcal W$-minimiser is the theorem of Balci–Ortner–Storn (2022) under their assumptions (A) (growth conditions; $x_T$ chosen so that the integrand is evaluated at a point of $T$, with $p(x)$ constant on each triangle here because the mesh resolves the interface); for the *two-phase* problem with discontinuous $a$ and $e^t$-type integrands we have not checked that the assumptions hold.
* Observed convergence of $F_h$ is slow and the difference between the two schemes still decreases on the finest meshes: the limits of the two sequences cannot be read off from four meshes.
* For discontinuous Galerkin we know of no convergence theorem in the presence of a gap; CR is the method with theory.
""")
nb["cells"] = C
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
out = root / "notebooks" / "glaplace_cr_lavrentiev.ipynb"
nbf.write(nb, out)
print("wrote", out)
