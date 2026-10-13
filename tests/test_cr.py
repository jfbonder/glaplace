"""Tests of the Crouzeix-Raviart variant (glaplace.cr, glaplace.lavrentiev)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
import glaplace as gl
from glaplace.cr import CRSpace, GLaplaceCRSolver, VarPower
from glaplace.lavrentiev import crisscross_mesh, G_pxlap, G_twophase, psi, energy_exact_pxlap_lambda1, u0


def test_cr_reproduces_linear_functions():
    sp = CRSpace(gl.rectangle_mesh(6))
    u = 2 * sp.mid[:, 0] - 3 * sp.mid[:, 1] + 1
    gx, gy = sp.gradient(u)
    assert np.allclose(gx, 2) and np.allclose(gy, -3)


def test_cr_poisson_order_two():
    ue = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)
    err = []
    for n in (8, 16, 32):
        S = GLaplaceCRSolver(gl.rectangle_mesh(n), 2.0, f=lambda x, y: 2 * np.pi**2 * ue(x, y), u_D=0.0, r=1.0)
        res = S.solve(tol=1e-11, stop="residual", max_iter=3000)
        assert res.converged
        err.append(S.fem.l2_error(res.u, ue))
    orders = np.log2(np.array(err[:-1]) / np.array(err[1:]))
    assert np.all(orders > 1.9)


def test_cr_energy_is_minimal():
    mesh = crisscross_mesh(8)
    S = GLaplaceCRSolver(mesh, G_twophase(), 0.0, psi(1.5), r=3.0)
    res = S.solve(tol=1e-10, stop="residual", max_iter=5000)
    assert res.converged
    J0 = S.energy(res.u)
    rng = np.random.default_rng(0)
    for _ in range(5):
        w = np.zeros_like(res.u); w[S.free] = rng.standard_normal(len(S.free))
        assert S.energy(res.u + 1e-3 * w) >= J0 - 1e-12


def test_lambda_one_exact_energy_both_schemes():
    G = G_pxlap()
    mesh = crisscross_mesh(8)
    for Solver in (gl.GLaplaceSolver, GLaplaceCRSolver):
        S = Solver(mesh, G, 0.0, psi(1.0), r=1.0)
        res = S.solve(tol=1e-10, stop="residual", max_iter=5000)
        assert res.converged
        assert abs(S.energy(res.u) - energy_exact_pxlap_lambda1()) < 1e-8


def test_gap_P1_energy_above_CR_energy():
    """lambda = 3: the conforming (P1) energy stays strictly above the nonconforming (CR) one (Lavrentiev gap)."""
    G = G_pxlap()
    mesh = crisscross_mesh(16)
    e = {}
    for name, Solver in (("P1", gl.GLaplaceSolver), ("CR", GLaplaceCRSolver)):
        S = Solver(mesh, G, 0.0, psi(3.0), r=3.0)
        res = S.solve(tol=1e-8, stop="residual", max_iter=8000)
        assert res.converged
        e[name] = S.energy(res.u)
    assert e["P1"] > e["CR"] + 1.0
