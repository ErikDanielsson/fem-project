import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import mesh
import numpy as np
import plantml

# Define the material properties
E_copper = 128
E_nylon = 3.00
ny_copper = 0.36
ny_nylon = 0.39
alpha_copper = 17.6e-6
alpha_nylon = 80e-6
rho_copper = 8930
rho_nylon = 1100
cp_copper = 386
cp_nylon = 1500
k_copper = 385
k_nylon = 0.26

D_copper = k_copper * np.eye(2)
D_nylon = k_nylon * np.eye(2)

# Now define the properties of the boundaries
heat_flux = 1e5
T_inf = 18

ep = [1, 1]

# Define the geometry and the mesh
g = mesh.define_geometry()
m = mesh.define_mesh(g, 0.1)
coords, edof, dofs, bdofs, element_markers, boundary_elements = m.create()


ndofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros((ndofs, ndofs))

# Assemble the stiffness matrix
for eldof, elx, ely, material_marker in zip(edof, ex, ey, element_markers):
    D = D_copper if material_marker == mesh.copper_marker else D_nylon
    Ke = cfc.flw2te(elx, ely, ep, D)
    cfc.assem(eldof, K, Ke)

# Apply convection
# Extract the boundary elements
convection_elements = boundary_elements[mesh.id_copper_outer]
for cel in convection_elements:
    nodes = cel["node-number-list"]


# Create the force vector
f = np.zeros((ndofs, 1))
cfu.apply_force(bdofs, f, mesh.id_thermal_load, heat_flux)
cfu.apply_force(bdofs, f, mesh.id_isolated, 0)
cfu.apply_force(bdofs, f, mesh.id_sym, 0)

# Create the boundary vectors assuming T_inf at boundaries
bc = np.array([], "i")
bcVal = np.array([], "f")
bc, bcVal = cfu.apply_bc(bdofs, bc, bcVal, mesh.id_copper_outer, 0)

a, q = cfc.solveq(K, f, bc, bcVal)

if __name__ == "__main__":
    cfv.figure(fig_size=(40, 30))
    cfv.subplot(221)
    cfv.draw_geometry(g)
    cfv.subplot(222)
    mesh.draw_mesh(m, coords, edof)
    cfv.subplot(223)
    cfv.draw_nodal_values_shaded(a, coords, edof, title="Temperature")
    cfv.colorbar()
    cfv.subplot(224)
    cfv.draw_nodal_values_shaded(q, coords, edof, title="Flux")
    cfv.colorbar()
    cfv.show()
