
import meshio
msh = meshio.read("CoaxWaveguideTransition.msh")
for cell in msh.cells:
    if  cell.type == "tetra":
        tetra_cells = cell.data

for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]

tetra_mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells},
                           cell_data={"VolumeRegions":[tetra_data]})

meshio.write("mesh.xdmf", tetra_mesh)



from dolfin import *
from vtkplotter.dolfin import datadir, plot

mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 3)
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mvc, "VolumeRegions")
cf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

info(mesh)

# Volume domains
File("VolSubDomains.pvd").write(cf)
File("Mesh.pvd").write(mesh)

plot(mesh)
