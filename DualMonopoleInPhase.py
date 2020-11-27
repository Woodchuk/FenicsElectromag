import meshio
from dolfin import *
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import os, sys, traceback
a = 2.25
b = 0.5;
d = 0.3;
s = 0.5;
l = 10.0;
w = 3.0;
pi = 3.1415926536
tol = 1.0e-12
eta = 377.0
Dk = 1.0 
qoffset = 0.0
poffset = 1.75

class PEC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class InputBC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], -1, tol)

class OutputBC(SubDomain):
    def inside(self, x, on_boundary):
        rb = sqrt(x[0] * x[0] + (x[1] + 1.25)* (x[1] + 1.25) + x[2] * x[2])
        return on_boundary and near(rb, 10.0, 5.0e-2)

class PMC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], 0.0, tol) or near(x[1], -1.25, tol))

mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 3)

info(mesh)
#plot(mesh)
#plt.show()

# Mark boundaries
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
sub_domains.set_all(4)
pec = PEC()
pec.mark(sub_domains, 0)
in_port = InputBC()
in_port.mark(sub_domains, 1)
out_port = OutputBC()
out_port.mark(sub_domains, 2)
pmc = PMC()
pmc.mark(sub_domains, 3)
File("BoxSubDomains.pvd").write(sub_domains)

dk = 0.025
#for m in range(40):
for m in range(1):
#   K = 0.500 + m * dk
   K = 0.5
   k0 = Constant(K)  
   beta = K
   b0 = Constant(beta)
   Kr = Dk

# Set up function spaces
# For low order problem
   cell = tetrahedron
   ele_type = FiniteElement('N1curl', cell, 2) # H(curl) element for EM
   V2 = FunctionSpace(mesh, MixedElement([ele_type, ele_type]))
   V = FunctionSpace(mesh, ele_type)
   u_r, u_i = TrialFunctions(V2)
   v_r, v_i = TestFunctions(V2)

#surface integral definitions from boundaries
   ds = Measure('ds', domain = mesh, subdomain_data = sub_domains)
# with source and sink terms
   u0 = Constant((0.0, 0.0, 0.0)) #PEC definition
   h_src = Expression(('-(x[1] - poffset) / (2.0 * pi * (pow((x[0]-qoffset), 2.0) + pow(x[1] - poffset,2.0)))', '(x[0]-qoffset) / (2.0 * pi *(pow((x[0]-qoffset),2.0) + pow(x[1] - poffset,2.0)))', 0.0), degree = 2, poffset = poffset, qoffset = qoffset)
   e_src = Expression(('(x[0] - qoffset) / (2.0 * pi * (pow((x[0]-qoffset), 2.0) + pow(x[1] - poffset,2.0)))', '(x[1]-poffset) / (2.0 * pi *(pow((x[0]-qoffset),2.0) + pow(x[1] - poffset,2.0)))', 0.0), degree = 2, poffset = poffset, qoffset = qoffset)
   Rrad = Expression(('sqrt(x[0] * x[0] + (x[1] + 1.25) * (x[1] + 1.25) + x[2] * x[2])'), degree = 2)
#Boundary condition dictionary
   boundary_conditions = {0: {'PEC' : u0},
                       1: {'InputBC': (h_src, Dk)},
                       2: {'OutputBC': Rrad},
                       3: {'PMC': 0.0}}

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
          bb1 = 2.0 * k0 * eta * inner(v_i, cross(n, r)) * ds(i) #Factor of two from field equivalence principle
          integral_source.append(bb1)
          bb2 = inner(cross(n, v_i), cross(n, u_r)) * b0 * ds(i)
          integrals_load.append(bb2)
          bb2 = inner(-cross(n, v_r), cross(n, u_i)) * b0 * ds(i)
          integrals_load.append(bb2)

   for i in boundary_conditions:
      if 'OutputBC' in boundary_conditions[i]:
        r = boundary_conditions[i]['OutputBC']
        bb2 = (inner(cross(n, v_i), cross(n, u_r)) * k0 + 0.5 * inner(cross(n, v_i), cross(n, u_i)) / r)* ds(i)
        integrals_load.append(bb2)
        bb2 = (inner(-cross(n, v_r), cross(n, u_i)) * k0 + 0.5 * inner(cross(n, v_r), cross(n, u_r)) / r)* ds(i)
        integrals_load.append(bb2)
# for PMC, do nothing. Natural BC.

   af = (inner(curl(v_r), curl(u_r)) + inner(curl(v_i), curl(u_i)) - Kr * k0 * k0 * (inner(v_r, u_r) + inner(v_i, u_i))) * dx + sum(integrals_load)
   L = sum(integral_source)

   u1 = Function(V2)
   vdim = u1.vector().size()
   print("Vdim = ", vdim)

   solve(af == L, u1, bcs, solver_parameters = {'linear_solver' : 'mumps'}) 

   u1_r, u1_i = u1.split(True)

   fp = File("EField_r.pvd")
   fp << u1_r
   fp = File("EField_i.pvd")
   fp << u1_i
   fp = File('SIWWaveFile.pvd')

   ut = u1_r.copy(deepcopy=True)
   for i in range(50):
      ut.vector().zero()
      ut.vector().axpy(-sin(pi * i / 25.0), u1_i.vector())
      ut.vector().axpy(cos(pi * i / 25.0), u1_r.vector()) 
      fp << (ut, i)

   H = interpolate(h_src, V) # Get input field
   P =  assemble((-dot(u1_r,cross(curl(u1_i),n))+dot(u1_i,cross(curl(u1_r),n))) * ds(2))
   P_refl = assemble((-dot(u1_i,cross(curl(u1_r), n)) + dot(u1_r, cross(curl(u1_i), n))) * ds(1))
   P_inc = assemble((dot(H, H) * 0.5 * eta * b0 / k0) * ds(1))
   print("k0 = ", K)
   print("Beta = ", beta)
   print("Integrated power on rad boundary:", P/(2.0 * K * eta))
   print("Incident power at port 1:", P_inc)
   print("Integrated reflected power on port 1:", P_inc - P_refl / (2.0 * K * eta))
   eps_c = 1.0
   lc = 1.0
   E = interpolate(e_src, V) # Incident E field
   ccr = assemble(-dot(u1_r - E * (eta / sqrt(eps_c)), E * (eta / sqrt(eps_c))) * ds(1))
   cci = assemble(dot(u1_i, E) * ds(1)) * eta / sqrt(eps_c)
   cc = assemble(dot(E, E) * ds(1)) * eta * eta / eps_c
   Zo = 50.0
   rho = complex(ccr / cc, cci / cc)
   print("Input port reflection coefficient: {0:<f}+j{1:<f}".format(rho.real, rho.imag))
   Zin = Zo * (1.0 + rho) / (1.0 - rho)
   print("Input port impedance: {0:<f} + j{1:<f}".format(Zin.real, Zin.imag))
   Zl = Zo * (Zin - (1j) * Zo * tan(K * sqrt(eps_c) * lc)) / (Zo - (1j) * Zin * tan(K * sqrt(eps_c) * lc))
   print("Antenna feedpoint impedance: {0:<f} + j{1:<f}".format(Zl.real, Zl.imag))

# Generate radiation pattern!
print("Generate radiation pattern.")
metadata = {"quadrature_degree": 6, "quadrature_scheme": "default"}
dsm = ds(metadata=metadata)
NumTheta = 25
NumPhi = 100
TwoPi = 6.2831853071
PiOverTwo = 1.5707963268
fp = open("MonoPattern1.txt", "w")
print("#Elevation    Azimuth     Pvert      Phoriz", file = fp)
# Reflection transformations
# PMC
JCx = Constant(((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))) # PMC Reflection thru x = 0 plane
#JCy = Constant(((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0))) # PEC Reflection thru y = -5 plane
JCy = Constant(((1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0))) # PMC Reflection thru y = -5 plane
MCx = Constant(((1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)))
#MCy = Constant(((1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)))
MCy = Constant(((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)))
# PEC ground plane
JCz = Constant(((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0))) #PEC reflection thru z = 0 plane
MCz = Constant(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)))

# Surface currents on external boundary
M_r = -cross(n, u1_r)
M_i = -cross(n, u1_i)
J_r = -cross(n, curl(u1_i)) / (k0 * eta)
J_i = cross(n, curl(u1_r)) / (k0 * eta)

for m in range(NumTheta+1):
    theta = m * PiOverTwo / NumTheta
    print(" ", file = fp) # for Gnuplot
    for nn in range(NumPhi+1):
        L_r = [] #List objects for integrals
        L_i = []
        N_r = []
        N_i = []
# Do NFF transformation
        phi = nn * TwoPi / NumPhi
        rr = Expression(('sin(theta)*cos(phi)', 'sin(theta)*sin(phi)', 'cos(theta)'), degree = 3, phi = phi, theta = theta)
        rtheta = Expression(('cos(theta)*cos(phi)', 'cos(theta)*sin(phi)', '-sin(theta)'), degree = 3, phi = phi, theta = theta)
        rphi = Expression(('-sin(phi)', 'cos(phi)', '0.0'), degree = 3, phi = phi)
# Sum up all the image sources taking into account the proper symmetries
# First octant
        rp1 = Expression(('x[0]', 'x[1]+1.25', 'x[2]'), degree = 1)
        sr = J_r * cos(k0 * dot(rr, rp1)) - J_i * sin(k0 * dot(rr, rp1))
        si = J_i * cos(k0 * dot(rr, rp1)) + J_r * sin(k0 * dot(rr, rp1))
        N_r.append(sr)
        N_i.append(si)
        qr = M_r * cos(k0 * dot(rr, rp1)) - M_i * sin(k0 * dot(rr, rp1))
        qi = M_i * cos(k0 * dot(rr, rp1)) + M_r * sin(k0 * dot(rr, rp1))
        L_r.append(qr)
        L_i.append(qi)
# Second octant x < 0, y > 0, z > 0
        rp2 = Expression(('-x[0]', 'x[1]+1.25', 'x[2]'), degree = 1)
        sr = JCx * (J_r * cos(k0 * dot(rr, rp2)) - J_i * sin(k0 * dot(rr, rp2)))
        si = JCx * (J_i * cos(k0 * dot(rr, rp2)) + J_r * sin(k0 * dot(rr, rp2)))
        N_r.append(sr)
        N_i.append(si)
        qr = MCx * (M_r * cos(k0 * dot(rr, rp2)) - M_i * sin(k0 * dot(rr, rp2)))
        qi = MCx * (M_i * cos(k0 * dot(rr, rp2)) + M_r * sin(k0 * dot(rr, rp2)))
        L_r.append(qr)
        L_i.append(qi)
# third octant x < 0, y < 0, z > 0
        rp3 = Expression(('-x[0]', '-x[1]-1.25', 'x[2]'), degree = 1)
        sr = JCy * JCx * (J_r * cos(k0 * dot(rr, rp3)) - J_i * sin(k0 * dot(rr, rp3)))
        si = JCy * JCx * (J_i * cos(k0 * dot(rr, rp3)) + J_r * sin(k0 * dot(rr, rp3)))
        N_r.append(sr)
        N_i.append(si)
        qr = MCy * MCx * (M_r * cos(k0 * dot(rr, rp3)) - M_i * sin(k0 * dot(rr, rp3)))
        qi = MCy * MCx * (M_i * cos(k0 * dot(rr, rp3)) + M_r * sin(k0 * dot(rr, rp3)))
        L_r.append(qr)
        L_i.append(qi)
# fourth octant x > 0, y < 0, z > 0
        rp4 = Expression(('x[0]', '-x[1]-1.25', 'x[2]'), degree = 1)
        sr = JCy * (J_r * cos(k0 * dot(rr, rp4)) - J_i * sin(k0 * dot(rr, rp4)))
        si = JCy * (J_i * cos(k0 * dot(rr, rp4)) + J_r * sin(k0 * dot(rr, rp4)))
        N_r.append(sr)
        N_i.append(si)
        qr = MCy * (M_r * cos(k0 * dot(rr, rp4)) - M_i * sin(k0 * dot(rr, rp4)))
        qi = MCy * (M_i * cos(k0 * dot(rr, rp4)) + M_r * sin(k0 * dot(rr, rp4)))
        L_r.append(qr)
        L_i.append(qi)
# Fifth octant x > 0, y > 0, z < 0
        rp5 = Expression(('x[0]', 'x[1]+1.25', '-x[2]'), degree = 1)
        sr = JCz * (J_r * cos(k0 * dot(rr, rp5)) - J_i * sin(k0 * dot(rr, rp5)))
        si = JCz * (J_i * cos(k0 * dot(rr, rp5)) + J_r * sin(k0 * dot(rr, rp5)))
        N_r.append(sr)
        N_i.append(si)
        qr = MCz * (M_r * cos(k0 * dot(rr, rp5)) - M_i * sin(k0 * dot(rr, rp5)))
        qi = MCz * (M_i * cos(k0 * dot(rr, rp5)) + M_r * sin(k0 * dot(rr, rp5)))
        L_r.append(qr)
        L_i.append(qi)
# Sixth octant x < 0, y > 0, z < 0
        rp6 = Expression(('-x[0]', 'x[1]+1.25', '-x[2]'), degree = 1)
        sr = JCx * JCz * (J_r * cos(k0 * dot(rr, rp6)) - J_i * sin(k0 * dot(rr, rp6)))
        si = JCx * JCz * (J_i * cos(k0 * dot(rr, rp6)) + J_r * sin(k0 * dot(rr, rp6)))
        N_r.append(sr)
        N_i.append(si)
        qr = MCx * MCz * (M_r * cos(k0 * dot(rr, rp6)) - M_i * sin(k0 * dot(rr, rp6)))
        qi = MCx * MCz * (M_i * cos(k0 * dot(rr, rp6)) + M_r * sin(k0 * dot(rr, rp6)))
        L_r.append(qr)
        L_i.append(qi)
# seventh octant x < 0, y < 0, z < 0
        rp7 = Expression(('-x[0]', '-x[1]-1.25', '-x[2]'), degree = 1)
        sr = JCy * JCx * JCz * (J_r * cos(k0 * dot(rr, rp7)) - J_i * sin(k0 * dot(rr, rp7)))
        si = JCy * JCx * JCz * (J_i * cos(k0 * dot(rr, rp7)) + J_r * sin(k0 * dot(rr, rp7)))
        N_r.append(sr)
        N_i.append(si)
        qr = MCy * MCx * MCz * (M_r * cos(k0 * dot(rr, rp7)) - M_i * sin(k0 * dot(rr, rp7)))
        qi = MCy * MCx * MCz * (M_i * cos(k0 * dot(rr, rp7)) + M_r * sin(k0 * dot(rr, rp7)))
        L_r.append(qr)
        L_i.append(qi)
# Eighth octant x > 0, y < 0, z < 0
        rp8 = Expression(('x[0]', '-x[1]-1.25', '-x[2]'), degree = 1)
        sr = JCy * JCz * (J_r * cos(k0 * dot(rr, rp8)) - J_i * sin(k0 * dot(rr, rp8)))
        si = JCy * JCz * (J_i * cos(k0 * dot(rr, rp8)) + J_r * sin(k0 * dot(rr, rp8)))
        N_r.append(sr)
        N_i.append(si)
        qr = MCy * MCz * (M_r * cos(k0 * dot(rr, rp8)) - M_i * sin(k0 * dot(rr, rp8)))
        qi = MCy * MCz * (M_i * cos(k0 * dot(rr, rp8)) + M_r * sin(k0 * dot(rr, rp8)))
        L_r.append(qr)
        L_i.append(qi)

# Compute E_ff
        Et_i = -K * assemble((dot(sum(L_r), rphi) + eta * dot(sum(N_r), rtheta)) * dsm(2))
        Et_r = K * assemble((dot(sum(L_i), rphi) + eta * dot(sum(N_i), rtheta)) * dsm(2))
        Ep_i = K * assemble((dot(sum(L_r), rtheta) - eta * dot(sum(N_r), rphi)) * dsm(2))
        Ep_r = -K * assemble((dot(sum(L_i), rtheta) - eta * dot(sum(N_i), rphi)) * dsm(2))

# Compute magnitudes
        Gvert = (Et_r * Et_r + Et_i * Et_i) / (2.0 * TwoPi * eta * (P * 8.0 / (2.0 * K * eta)))
        Ghoriz = (Ep_r * Ep_r + Ep_i * Ep_i) / (2.0 * TwoPi * eta * (P * 8.0 / (2.0 * K * eta)))

        print(" {0:f} {1:f} {2:f} {3:f}".format(theta, phi, Gvert, Ghoriz))
        print(" {0:f} {1:f} {2:f} {3:f}".format(theta, phi, Gvert, Ghoriz), file = fp)

fp.close()
sys.exit(0)

