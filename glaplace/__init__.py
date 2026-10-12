"""glaplace: g-Laplacian solver (decomposition-coordination / ADMM + P1 FEM) on general 2D domains.

-div( g(x,|grad u|) grad u/|grad u| ) = f,  g = dG/dt,  G a Young function (p-q, p log^a, exp, two-phase, ...).
Built from the p(x)-Laplacian code `pxlaplace` (Aragon, Fernandez Bonder, Rubio).
"""
from .mesh import (Mesh, rectangle_mesh, polygon_mesh, sdf_mesh, disk_mesh,
                   from_txt, from_meshio, boundary_nodes_from_triangles)
from .fem import P1Space
from .young import (Young, PowerLaw, PLog, Exp, ExpLiteral, Custom, Sum, from_expr, pq, check_young)
from .solver import GLaplaceSolver, solve, SolveResult
from .nonlinear import solve_scalar, update_nu
from . import exact

__all__ = ["Mesh", "rectangle_mesh", "polygon_mesh", "sdf_mesh", "disk_mesh", "from_txt", "from_meshio",
           "boundary_nodes_from_triangles", "P1Space", "Young", "PowerLaw", "PLog", "Exp", "ExpLiteral", "Custom",
           "Sum", "from_expr", "pq", "check_young", "GLaplaceSolver", "solve", "SolveResult", "solve_scalar",
           "update_nu", "exact"]
