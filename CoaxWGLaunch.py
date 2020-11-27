from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os, sys, traceback

a = 2.286
b = 1.016
w = 0.554
pi = 3.1415926536
tol = 1.0e-12
eta = 377.0
#k0 = 2.178  # 10.4GHz
k0 = 1.92
beta = sqrt(k0 * k0 - (pi / a) **2.0)
Zm = eta * k0 / beta
eps0 = 1.0  # WG air filled
eps_c = 2.2 # Teflon coax dielectric

class PEC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class InputBC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], -1.0, tol)

class OutputBC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 4.0, tol)


mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 3)

info(mesh)
#plot(mesh)
#plt.show()

# Mark boundaries
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
sub_domains.set_all(3)
pec = PEC()
pec.mark(sub_domains, 0)
in_port = InputBC()
in_port.mark(sub_domains, 1)
out_port = OutputBC()
out_port.mark(sub_domains, 2)
File("BoxSubDomains.pvd").write(sub_domains)

Dk = Expression('x[1] <= 0.0 + tol ? eps_c : eps0', degree = 0, tol = 1.0e-12, eps_c = eps_c, eps0 = eps0) # Dielectric subdomains

# Set up function spaces
cell = tetrahedron
ele_type = FiniteElement('N1curl', cell, 3) # H(curl) element for EM
V2 = FunctionSpace(mesh, ele_type * ele_type)
V = FunctionSpace(mesh, ele_type)
(u_r, u_i) = TrialFunctions(V2)
(v_r, v_i) = TestFunctions(V2)

#surface integral definitions from boundaries
ds = Measure('ds', domain = mesh, subdomain_data = sub_domains)
# with source and sink terms
u0 = Constant((0.0, 0.0, 0.0)) #PEC definition

h_src = Expression(('-(x[2] - w) / (2.0 * pi * (pow(x[0] - a/2.0, 2.0) + pow(x[2] - w,2.0)))', 0.0, '(x[0] - a/2.0) / (2.0 * pi *(pow(x[0] - a/2.0,2.0) + pow(x[2] - w,2.0)))'), degree = 2, a = a, w = w)
e_src = Expression(('(x[1]<0.0)?(x[0] - w) / (2.0 * pi * (pow(x[0] - a/2.0, 2.0) + pow(x[2] - w,2.0))):0.0', 0.0, '(x[1]<0.0)?(x[2] - a/2.0) / (2.0 * pi *(pow(x[0] - a/2.0,2.0) + pow(x[2] - w,2.0))):0.0', 0.0, 0.0, 0.0), degree = 2, a = a, w = w)
#Boundary condition dictionary
boundary_conditions = {0: {'PEC' : u0},
                       1: {'InputBC': (h_src, eps_c)},
                       2: {'OutputBC': 1.0}}

n = FacetNormal(mesh)

#Build PEC boundary conditions for real and imaginary parts
bcs = []
for i in boundary_conditions:
    if 'PEC' in boundary_conditions[i]:
        bc = DirichletBC(V2.sub(0), boundary_conditions[i]['PEC'], sub_domains, i)
        bcs.append(bc)
        bc = DirichletBC(V2.sub(1), boundary_conditions[i]['PEC'], sub_domains, i)
        bcs.append(bc)

# Build input BC source term and loading term
integral_source = []
integrals_load =[]
for i in boundary_conditions:
    if 'InputBC' in boundary_conditions[i]:
        r, s = boundary_conditions[i]['InputBC']
        bb1 = 2.0 * (k0 * eta) * inner(v_i, cross(n, r)) * ds(i) #Factor of two from field equivalence principle
        integral_source.append(bb1)
        bb2 = inner(cross(n, v_i), cross(n, u_r)) * k0 * sqrt(Dk) * ds(i)
        integrals_load.append(bb2)
        bb2 = inner(-cross(n, v_r), cross(n, u_i)) * k0 * sqrt(Dk) * ds(i)
        integrals_load.append(bb2)

for i in boundary_conditions:
    if 'OutputBC' in boundary_conditions[i]:
        r = boundary_conditions[i]['OutputBC']
        bb2 = inner(cross(n, v_i), cross(n, u_r)) * beta * ds(i)
        integrals_load.append(bb2)
        bb2 = inner(-cross(n, v_r), cross(n, u_i)) * (beta) * ds(i)
        integrals_load.append(bb2)

a = (inner(curl(v_r), curl(u_r)) + inner(curl(v_i), curl(u_i)) - Dk * k0 * k0 * (inner(v_r, u_r) + inner(v_i, u_i))) * dx + sum(integrals_load)
L = sum(integral_source)

u1 = Function(V2)
#u1.vector()[:] = interpolate(e_src, V2).vector()
vdim = u1.vector().size()
print("Solution vector size =", vdim)
#u1.vector()[:] = np.random.uniform(-1, 1, vdim)
solve(a == L, u1, bcs, solver_parameters = {'linear_solver' : 'mumps'}) 

#A, b = assemble_system(a, L, bcs)
#solver = PETScKrylovSolver("gmres","none")

#solver.parameters['absolute_tolerance'] = 1e-10
#solver.parameters['relative_tolerance'] = 1e-6
#solver.parameters['maximum_iterations'] = 10000
#solver.parameters['monitor_convergence'] = True
#solver.parameters['nonzero_initial_guess'] = True 
#solver.ksp().setGMRESRestart(1000)
#solver.solve(A,u1.vector(),b)

u1_r, u1_i = u1.split(True)
fp = File("EField_r.pvd")
fp << u1_r
fp = File("EField_i.pvd")
fp << u1_i

fp = File('WaveFile.pvd')

ut = u1_r.copy(deepcopy=True)
for i in range(50):
    ut.vector().zero()
    ut.vector().axpy(cos(pi * i / 25.0 + pi / 2.0), u1_i.vector())
    ut.vector().axpy(cos(pi * i / 25.0), u1_r.vector()) 
    fp << (ut, i)
H = interpolate(h_src, V) # Get input field
P =  assemble((-dot(u1_r,cross(curl(u1_i),n))+dot(u1_i,cross(curl(u1_r),n))) * ds(2))
#P_refl = assemble((dot((u1_i-sqrt(1.0 / Dk)*eta*cross(n,H)),cross(curl(u1_r), n)) + dot(u1_r, cross(curl(u1_i - sqrt(1.0 / Dk)*eta*cross(n, H)), n))) * ds(1))
P_refl = assemble((-dot(u1_i,cross(curl(u1_r), n)) + dot(u1_r, cross(curl(u1_i), n))) * ds(1))
P_inc = assemble((dot(H, H) * eta / (2.0 * sqrt(eps_c))) * ds(1))
print("Integrated power on port 2:", P/(2.0 * k0 * eta))
print("Incident power at port 1:", P_inc)
print("Integrated reflected power on port 1:", P_inc - P_refl / (2.0 * k0 * eta))
sys.exit(0)

