import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import numpy as np

# Begin by defining the geometry of the problem

L = 5
a = 0.1 * L
b = 0.1 * L
c = 0.3 * L
d = 0.05 * L
h = 0.15 * L
t = 0.05 * L

g = cfg.Geometry()

# Begin with the points for the nylon part
g.point([0, 0])  # 0
g.point([c + d, 0])  # 1
g.point([c + d, 0.3 * L])  # 2
g.point([a + t, 0.3 * L])  # 3
g.point([a + t, 0.3 * L - h])  # 4
g.point([a, 0.3 * L - h])  # 5
g.point([a, 0.3 * L])  # 6
g.point([0, 0.3 * L])  # 7

# Now the points for the copper part
g.point([a + c + d, 0])  # 8
g.point([L - 2 * d, 0.3 * L - d])  # 9
g.point([L, 0.3 * L - d])  # 10
g.point([L, 0.3 * L])  # 11
g.point([L - 2 * d, 0.3 * L])  # 12
g.point([a + c + d, d])  # 13
g.point([a + c + d, 0.5 * L - a - d])  # 13
g.point([a + c, 0.5 * L - a])  # 14
g.point([a, 0.5 * L - a])  # 15
g.point([a, 0.5 * L])  # 16
g.point([0, 0.5 * L])  # 17
g.point([0, 0.5 * L - a])  # 18

# Now define the markers and create the splines
id_nylon_copper = 10  # Nylon copper boundary
id_copper_outer = 20  # Copper outer boundary
id_sym = 30  # Symmetry boundary
id_isolated = 40  # No heat flow boundary
id_thermal_load = 50  # Thermal load boundary

# Create splines for the nylon area
g.spline([1, 0], marker=id_isolated)
g.spline([2, 1], marker=id_nylon_copper)
g.spline([3, 2], marker=id_nylon_copper)
g.spline([4, 3], marker=id_nylon_copper)
g.spline([5, 4], marker=id_nylon_copper)
g.spline([6, 5], marker=id_nylon_copper)
g.spline([7, 6], marker=id_nylon_copper)
g.spline([0, 7], marker=id_isolated)

# Create splines for the copper area
g.spline([1, 8], marker=id_isolated)
g.spline([8, 9], marker=id_copper_outer)
g.spline([9, 10], marker=id_copper_outer)
g.spline([10, 11], marker=id_sym)
g.spline([11, 12], marker=id_copper_outer)
g.spline([12, 13], marker=id_copper_outer)
g.spline([13, 14], marker=id_copper_outer)
g.spline([14, 15], marker=id_copper_outer)
g.spline([15, 16], marker=id_copper_outer)
g.spline([16, 17], marker=id_copper_outer)
g.spline([17, 18], marker=id_sym)
g.spline([18, 19], marker=id_thermal_load)
g.spline([19, 7], marker=id_isolated)

g.surface([0, 1, 2, 3, 4, 5, 6, 7])
g.surface([6, 5, 4, 3, 2, 1, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

cfv.draw_geometry(g)

mesh = cfm.GmshMesh(g)
mesh.el_type = 2
mesh.dofs_per_node = 1  # Degrees of freedom per node.
mesh.el_size_factor = 0.05  # Factor that changes element sizes.

coords, edof, dofs, bdofs, element_markers = mesh.create()

ndofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros((ndofs, ndofs))
cfv.draw_mesh(
    coords=coords,
    edof=edof,
    dofs_per_node=mesh.dofs_per_node,
    el_type=mesh.el_type,
    filled=True,
    title="Example 01",
)

cfv.show()
