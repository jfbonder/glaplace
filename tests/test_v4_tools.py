"""Tests of glaplace.v4tools (meshes resolving circles, singular quadrature, element min/max)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
from glaplace.v4tools import disk_mesh_resolved, mesh_h, SingularQuadrature, element_minmax, a_elm
from glaplace.lavrentiev import crisscross_mesh


def test_disk_mesh_resolves_circle_and_refines():
    hs, nn = (0.1, 0.05, 0.0125), []
    for h in hs:
        m = disk_mesh_resolved(h)
        r = np.hypot(*m.vertices.T)
        assert np.any(r < 1e-14)                                        # vertex at the origin
        on = np.abs(r - 0.5) < 1e-12
        assert on.sum() >= int(np.ceil(np.pi / h))                       # vertices on r = 1/2
        # every triangle lies on one side of the polygon r = 1/2
        rv = r[m.triangles]
        assert np.all((rv <= 0.5 + 1e-12).all(1) | (rv >= 0.5 - 1e-12).all(1))
        assert mesh_h(m) < 1.6 * h
        nn.append(m.n_nodes)
    assert nn[2] > 10 * nn[1]                                            # the area switch is honoured on fine meshes


def test_singular_quadrature():
    m = disk_mesh_resolved(0.1, circles=())
    k = int(np.ceil(2 * np.pi / 0.1))
    Q0 = SingularQuadrature(m, gamma=0.0)
    assert abs(Q0.integrate(np.ones_like(Q0.W)).sum() - 0.5 * k * np.sin(2 * np.pi / k)) < 1e-12
    Q = SingularQuadrature(m, gamma=1.5)
    val = Q.integrate(np.hypot(Q.X[..., 0], Q.X[..., 1]) ** -1.5).sum()
    th = np.linspace(-np.pi / k, np.pi / k, 200001); R = np.cos(np.pi / k) / np.cos(th)
    assert abs(val - k * np.trapezoid(2 * np.sqrt(R), th)) < 1e-8
    u = 1 + m.vertices[:, 0] - 2 * m.vertices[:, 1]
    assert np.allclose(Q.eval_P1(m, u), 1 + Q.X[..., 0] - 2 * Q.X[..., 1])


def test_element_minmax_brackets_samples():
    cm = crisscross_mesh(16); a = a_elm(1.0)
    mn, mx = element_minmax(a, cm)
    P = cm.vertices[cm.triangles]
    L = np.random.default_rng(0).dirichlet([1, 1, 1], 3000)
    A = a(*np.einsum("qi,tid->dtq", L, P))
    assert np.all(A.min(1) >= mn - 1e-13) and np.all(A.max(1) <= mx + 1e-13)
