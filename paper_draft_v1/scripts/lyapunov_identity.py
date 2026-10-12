"""Numerical check of the Lyapunov identity (Theorem A, eq. (lyap)) for the discrete iteration.

For each case we run L iterations of the algorithm (nu_0 = eta_1 = 0) on a disk mesh and build a discrete
saddle-point proxy from the last iterate:  u* = u_L,  nu* = grad u_L,  eta* = eta_{L+1} + r (nu_L - nu_{L-1}).
By construction nu* = grad_h u* and (eta*, grad w)_h = <f, w>_h hold EXACTLY (up to round-off), which is all the
identity uses; eta* in dF(nu*) holds only approximately (this affects the sign of m_n only at the proxy-error level).
We report, over n = 1..L-500:
  - max |identity residual| / E_0,
  - max (E_n - E_{n-1}) / E_0 (should be <= 0),
  - min cross term  r (s_n, nu_n - nu_{n-1}) / E_0 (>= 0 exactly by monotonicity of the prox),
  - min m_n / E_0.
Run: PYTHONPATH=/home/claude/glaplace python lyapunov_identity.py
"""
import numpy as np, json, sys
import glaplace as gl
from glaplace.nonlinear import update_nu

def run(G, f, h, r, L):
    m = gl.disk_mesh(h); S = gl.GLaplaceSolver(m, G, f, 0.0, r=r); fem = S.fem
    w = fem.area
    ip = lambda a, b: float(np.sum(w * (a[0] * b[0] + a[1] * b[1])))
    sub = lambda a, b: (a[0] - b[0], a[1] - b[1])
    Nt = fem.Nt
    nu = (np.zeros(Nt), np.zeros(Nt)); eta = (np.zeros(Nt), np.zeros(Nt))
    H = []   # (nu_{n-1}, nu_n, eta_n, eta_{n+1}, grad u_n)
    for n in range(1, L + 1):
        rhs = S.F + fem.weak_div(r * nu[0] - eta[0], r * nu[1] - eta[1]); u = S._solve_linear(rhs)
        gu = fem.gradient(u)
        nun = update_nu(S._gfun, eta[0] / r + gu[0], eta[1] / r + gu[1], r)
        etan = (eta[0] + r * (gu[0] - nun[0]), eta[1] + r * (gu[1] - nun[1]))
        H.append((nu, nun, eta, etan, gu))
        nu, eta = nun, etan
    nuL_prev, nuL, _, etaL1, guL = H[-1]
    nus = guL
    etas = (etaL1[0] + r * (nuL[0] - nuL_prev[0]), etaL1[1] + r * (nuL[1] - nuL_prev[1]))
    # check of the discrete equation for eta*: (eta*, grad phi_i) = F_i on free nodes
    eq = fem.weak_div(*etas) - S.F
    eqres = np.abs(eq[S.free]).max() / np.abs(S.F[S.free]).max()
    E = lambda nu_, eta1: ip(sub(eta1, etas), sub(eta1, etas)) / (2 * r) + r * ip(sub(nu_, nus), sub(nu_, nus)) / 2
    E0 = E(H[0][0], H[0][2])          # E_0 = E(nu_0, eta_1)
    En, res, cross, mn = [], [], [], []
    Eprev = E0
    for (nup, nun, etan, etan1, gu) in H:
        s = sub(gu, nun)
        Ecur = E(nun, etan1)
        m_ = ip(sub(etan1, etas), sub(nun, nus))
        cr = r * ip(s, sub(nun, nup))
        lhs = Ecur + r / 2 * ip(s, s) + r / 2 * ip(sub(nun, nup), sub(nun, nup)) + m_ + cr
        res.append(lhs - Eprev); En.append(Ecur); cross.append(cr); mn.append(m_)
        Eprev = Ecur
    En, res, cross, mn = map(np.array, (En, res, cross, mn))
    k = L - 500
    dE = np.diff(np.concatenate([[E0], En]))[:k]
    En_prev = np.concatenate([[E0], En[:-1]])
    act = np.flatnonzero(En_prev / E0 > 1e-20)          # iterations n with E_{n-1}/E_0 > 1e-20
    na = int(act.max()) + 1
    restricted = dict(n_active=na,
                      max_identity_res_active=float(np.abs(res[:na]).max() / E0),
                      max_increase_active=float(dE[:na].max() / E0),
                      min_cross_active=float(cross[:na].min() / E0),
                      min_m_active=float(mn[:na].min() / E0))
    return dict(**restricted, nodes=m.n_nodes, E0=E0, E_k=float(En[k - 1]), k=k,
                eq_residual=float(eqres),
                max_identity_res=float(np.abs(res[:k]).max() / E0),
                max_increase=float(dE.max() / E0),
                min_cross=float(cross[:k].min() / E0),
                min_m=float(mn[:k].min() / E0),
                final_residual=float(np.sqrt(ip(sub(H[-1][4], H[-1][1]), sub(H[-1][4], H[-1][1]))))), En / E0

if __name__ == "__main__":
    cases = [("pq24", "t^2/2+t^4/4", gl.pq(2, 4), 1.0, 1.0, 1500),
             ("p1.5", "t^{1.5}/1.5", gl.PowerLaw(1.5), 1.0, 1.0, 2500),
             ("p4", "t^4/4", gl.PowerLaw(4.0), 1.0, 1.0, 1500),
             ("exp", "e^t-1-t", gl.Exp(), 1.0, 1.0, 1500),
             ("expdead", "e^t-1, c=8", gl.ExpLiteral(), 8.0, 30.0, 2500)]
    out = {}; curves = {}
    for name, lab, G, f, r, L in cases:
        d, curve = run(G, f, 0.1, r, L)
        d["label"] = lab; d["r"] = r; d["L"] = L
        out[name] = d; curves[name] = curve
        print(name, json.dumps(d), flush=True)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "lyapunov.json", "w"), indent=1)
    np.savez(sys.argv[2] if len(sys.argv) > 2 else "lyapunov_curves.npz", **curves)
