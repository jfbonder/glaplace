import numpy as np, pytest
import glaplace as gl
from glaplace import (rectangle_mesh, disk_mesh, polygon_mesh, P1Space, GLaplaceSolver, PowerLaw, PLog, Exp, ExpLiteral,
                      pq, check_young, from_expr, exact)
from glaplace.nonlinear import solve_scalar, update_nu

ring = lambda x, y: np.where(np.hypot(x, y) > 0.5, 1.0, 0.0)


def test_young_algebra_and_checks():
    G = pq(2, 4)
    t = np.array([0.0, 0.5, 2.0])
    assert np.allclose(G.G(t), t**2 / 2 + t**4 / 4) and np.allclose(G.g(t), t + t**3)
    a = lambda x, y: x                                           # x-dependent weight
    H = PowerLaw(2) + a * PowerLaw(3)
    x = np.array([0.0, 2.0, 3.0]); T = np.array([1.0, 1.0, 2.0])
    assert np.allclose(H.g(T, x, 0 * x), T + x * T**2)
    assert np.allclose((2.0 * PowerLaw(2) + np.float64(3.0) * PowerLaw(3)).g(T), 2 * T + 3 * T**2)
    for Gf in (PowerLaw(1.3), PLog(2, 1), PLog(3, -1), PLog(1.2, -0.1), Exp(), ExpLiteral(), pq(1.5, 3, 0.3),
               from_expr("t**3/3 + t**2/2")):
        assert check_young(Gf, verbose=False), Gf
    assert not check_young(from_expr("-t**2/2 + t"), verbose=False)       # g not >= 0 / not monotone


@pytest.mark.parametrize("Gf", [PowerLaw(1.05), PowerLaw(50), pq(2, 4), PLog(2, 1), PLog(3, -1), Exp(), ExpLiteral()])
def test_scalar_equation(Gf):
    s = np.logspace(-12, 6, 60)
    gfun = lambda t, idx=None: Gf.g(t)
    t = solve_scalar(gfun, s, 1.0)
    pos = t > 0
    assert np.allclose((Gf.g(t) + t)[pos], s[pos], rtol=1e-9, atol=0)
    if isinstance(Gf, ExpLiteral):                                        # g(0+) = 1: root is 0 for s <= 1
        assert np.all(t[s <= 1.0] == 0) and np.all(t[s > 1.0 + 1e-6] > 0)
    else:
        assert pos.all()


def test_p2_is_poisson():
    m = rectangle_mesh(40)
    ex = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)
    res = GLaplaceSolver(m, 2.0, lambda x, y: 2 * np.pi**2 * ex(x, y), 0.0).solve(tol=1e-10)
    assert res.converged and np.abs(res.u - ex(*m.vertices.T)).max() < 2e-3


CASES = [("p4", PowerLaw(4), 1.0, 1.0), ("pq", pq(2, 4), 1.0, 1.0), ("plog", PLog(2, 1), 1.0, 1.0),
         ("exp", Exp(), 1.0, 1.0), ("two-phase", pq(2, 4, ring), 1.0, 1.0)]


@pytest.mark.parametrize("name,Gf,c,r", CASES)
def test_radial_torsion(name, Gf, c, r):
    u, gu, _ = exact.radial_torsion(Gf, c)
    errs = []
    for h in (0.1, 0.05):
        m = disk_mesh(h)
        res = GLaplaceSolver(m, Gf, c, 0.0, r=r).solve(tol=1e-9, max_iter=2000, stop="residual")
        assert res.converged
        errs.append(P1Space(m).l2_error(res.u, u))
    assert errs[1] < 2e-3 and np.log(errs[0] / errs[1]) / np.log(2) > 1.5          # ~ order 2


def test_dead_core_exp_literal():
    Gf = ExpLiteral(); u, gu, ur = exact.radial_torsion(Gf, 8.0)
    assert ur(np.array([0.2]))[0] == 0 and ur(np.array([0.4]))[0] > 0              # dead core r < 1/4
    m = disk_mesh(0.05)
    res = GLaplaceSolver(m, Gf, 8.0, 0.0, r=3.0).solve(tol=1e-7, max_iter=3000, stop="residual")
    assert res.converged and P1Space(m).l2_error(res.u, u) < 3e-3
    gx, gy = P1Space(m).gradient(res.u)
    core = np.hypot(*P1Space(m).bary.T) < 0.15
    assert np.hypot(gx, gy)[core].max() < 1e-4                                      # flat in the core


def test_manufactured_x_dependent_G():
    f, u, gu = exact.manufactured("t**2/2 + (1+x)*t**4/4", "exp(x+y)")
    Gf = pq(2, 4, lambda x, y: 1 + x)
    m = rectangle_mesh(20, 20)
    res = GLaplaceSolver(m, Gf, f, u, r=30.0, f_degree=2).solve(tol=1e-9, max_iter=2000, stop="residual")
    assert res.converged and P1Space(m).l2_error(res.u, u) < 5e-4


def test_minimiser_of_discrete_energy():
    m = disk_mesh(0.1); Gf = PLog(2, 1)
    S = GLaplaceSolver(m, Gf, 1.0, 0.0)
    res = S.solve(tol=1e-11, max_iter=3000, stop="residual")
    E0 = S.energy(res.u)
    rng = np.random.default_rng(0)
    for _ in range(5):
        v = np.zeros(m.n_nodes); v[S.free] = rng.standard_normal(len(S.free))
        assert S.energy(res.u + 1e-3 * v) >= E0 - 1e-12 and S.energy(res.u - 1e-3 * v) >= E0 - 1e-12


def test_two_phase_sector_maximum_principle_and_init():
    th = np.linspace(0, 2 * np.pi, 60, endpoint=False)
    m = polygon_mesh(np.c_[np.cos(th), np.sin(th)], max_area=0.01)
    a = lambda x, y: np.where(np.arctan2(y, x) > 0.5, 1.0, 0.0)      # vanishes on a sector
    Gf = pq(2, 4, a)
    S = GLaplaceSolver(m, Gf, 1.0, 0.0)
    r1 = S.solve(tol=1e-9, max_iter=3000, stop="residual")
    r2 = S.solve(tol=1e-9, max_iter=3000, stop="residual", init="random")
    assert r1.converged and r2.converged and np.abs(r1.u - r2.u).max() < 1e-6 and r1.u.min() > -1e-8


def test_true_errors_quadrature():
    m = rectangle_mesh(10)
    F = P1Space(m)
    ex = lambda x, y: x**2
    assert F.l2_error(ex(*m.vertices.T), ex, relative=False) > 0         # P1 interpolation of x^2 is not exact
