"""Tests of the v4 modules (glaplace.elementwise, glaplace.admm): fast pointwise step, conjugates, discrete duality gap,
guaranteed lower bound of Theorem D (NA_REQUESTS.md section 0, tests (a)-(c))."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
import pytest
import glaplace as gl
from glaplace.cr import CRSpace, GLaplaceCRSolver
from glaplace.lavrentiev import crisscross_mesh, G_pxlap, G_twophase, psi
from glaplace.elementwise import PowerSum, ExpMinus1MinusT, ExpMinus1, Generic, from_young
from glaplace.nonlinear import solve_scalar
from glaplace.admm import Discretisation, ADMM, prolong


def _disc(n, lam, scheme):
    m = crisscross_mesh(n)
    if scheme == "P1":
        gD = psi(lam)(*m.vertices[m.dirichlet].T)
    else:
        fem = CRSpace(m); gD = fem.edge_values(psi(lam), fem.bdry_edges)
    return Discretisation(m, scheme, gD)


def test_prox_matches_bisection():
    rng = np.random.default_rng(1); N = 5000
    s = 10 ** rng.uniform(-10, 5, N)
    fams = [PowerSum(np.ones((1, N)), rng.uniform(1.05, 6, (1, N))),
            PowerSum(np.vstack([np.ones(N), rng.integers(0, 2, N)]), np.vstack([np.full(N, 1.5), np.full(N, 3.0)])),
            ExpMinus1MinusT(N), ExpMinus1(N)]
    for F in fams:
        for r in (0.1, 3.0, 100.0):
            t = F.prox(s, r)
            tb = solve_scalar(lambda tt, ii: F.subset(ii).g(tt), s, r)
            assert np.all((t == 0) == (tb == 0))
            pos = t > 0
            assert np.max(np.abs(t[pos] - tb[pos]) / tb[pos]) < 1e-11
            # warm start gives the same root
            tw = F.prox(s, r, t0=t * rng.uniform(0.01, 100, N))
            assert np.max(np.abs(tw[pos] - t[pos]) / t[pos]) < 1e-12


def test_conjugates_are_upper_bounds_and_sharp():
    rng = np.random.default_rng(2); N = 400
    s = 10 ** rng.uniform(-3, 2, N)
    tgrid = np.concatenate([[0.0], np.logspace(-8, 5, 40001)])
    fams = [PowerSum(np.ones((1, N)), rng.uniform(1.5, 5, (1, N))),
            PowerSum(np.vstack([np.ones(N), rng.integers(0, 2, N) * 1.0]), np.vstack([np.full(N, 1.5), np.full(N, 3.0)])),
            ExpMinus1MinusT(N), ExpMinus1MinusT(N, w=rng.uniform(0, 2, N), beta=1.0), ExpMinus1(N),
            Generic(gl.pq(1.5, 3.0), np.zeros(N), np.zeros(N))]
    for F in fams:
        c = F.conj(s)
        brute = np.array([np.max(s[i] * tgrid - F.subset(np.full(tgrid.size, i)).G(tgrid)) for i in range(N)]) \
            if not isinstance(F, Generic) else np.array([np.max(s[i] * tgrid - F.Gf.G(tgrid)) for i in range(N)])
        assert np.all(c >= brute - 1e-12 * np.maximum(1, np.abs(brute)))          # upper bound (rigour)
        assert np.max(np.abs(c - brute) / np.maximum(1, np.abs(brute))) < 1e-4    # and sharp (grid resolution)


def test_admm_matches_existing_solvers():
    m = crisscross_mesh(8)
    for Solver, scheme in ((gl.GLaplaceSolver, "P1"), (GLaplaceCRSolver, "CR")):
        S = Solver(m, G_twophase(), 0.0, psi(2.0), r=3.0)
        ref = S.solve(tol=1e-12, stop="residual", max_iter=20000)
        o = ADMM(_disc(8, 2.0, scheme), G_twophase(), r=3.0).solve(tol_res=1e-12, tol_gap=1e-12)
        assert np.max(np.abs(o["u"] - ref.u)) < 1e-9
        assert abs(o["J"] - S.energy(ref.u)) < 1e-10


def test_a_exact_value_two():
    """(a) p(x), lambda = 1, n = 8: D(sigma_final) = 2 to 1e-12 and F_h(P1) = 2."""
    for scheme in ("P1", "CR"):
        o = ADMM(_disc(8, 1.0, scheme), G_pxlap(), r=3.0).solve(jumps=True)
        assert o["converged"]
        assert abs(o["J"] - 2.0) < 1e-12
        assert abs(o["D"] - 2.0) < 1e-12


@pytest.mark.parametrize("G,lam", [(G_pxlap(), 3.0), (G_twophase(), 3.0), (G_twophase(), 0.5)])
def test_b_c_jumps_and_bounds_every_iteration(G, lam):
    """(b) max normal jump < 1e-11 at every iteration; (c) D(sigma_n) <= F_h(u_n^CR) + 1e-9 and D(sigma_n) <= J_h(u^P1)."""
    dc = _disc(8, lam, "CR")
    A = ADMM(dc, G, r=3.0)
    o = A.solve(jumps=True)
    H = o["history"]
    assert o["converged"]
    assert np.all(H["jump"] < 1e-11)
    assert np.all(H["D"] <= H["J"] + 1e-9)
    # Theorem D boundary-flux formula == algebraic discrete dual
    sx, sy = o["sigma"]
    Dna = dc.boundary_flux_term(sx, sy) - np.sum(dc.area * A.G.conj(np.hypot(sx, sy)))
    assert abs(Dna - o["D"]) < 1e-11
    p1 = ADMM(_disc(8, lam, "P1"), G, r=3.0).solve()
    assert np.all(H["D"] <= p1["J"] + 1e-12)
    assert np.all(p1["history"]["J"] >= p1["J"] - 1e-9)   # all P1 iterates are upper bounds; final one is the smallest


def test_table0_reference_values():
    """Reproduce two entries of NA Table 0 (lb_check_reference.py)."""
    o = ADMM(_disc(16, 3.0, "CR"), G_pxlap(), r=3.0).solve()
    assert abs(o["D"] - 17.3576790) < 1e-6
    o = ADMM(_disc(16, 3.0, "P1"), G_twophase(), r=3.0).solve()
    assert abs(o["J"] - 24.0669252) < 1e-6


def test_prolongation_warm_start():
    mc, mf = crisscross_mesh(8), crisscross_mesh(16)
    oc = ADMM(_disc(8, 3.0, "CR"), G_pxlap(), r=3.0).solve()
    nu = prolong(mc, mf, list(oc["nu"])); eta = prolong(mc, mf, list(oc["eta"]))
    A = ADMM(_disc(16, 3.0, "CR"), G_pxlap(), r=3.0)
    cold, warm = A.solve(), A.solve(nu0=nu, eta1=eta)
    assert warm["converged"] and abs(warm["D"] - cold["D"]) < 1e-9
