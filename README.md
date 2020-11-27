# FenicsElectromag
These files have been successfully run on a Dell Optiplex 3060 with 4 cores and 32GB of RAM under Ubuntu 18.04

Monopole2 - Fenics simulation files for a single monopole antenna\
 *.geo = GMSH mesh generator command file\
 *.msh = Finite element mesh file generated by GMSH\
 Monopole2Conv.py = Mesh converter file:  Run "python3 Monopole2Conv.py" before running solver\
 Monopole2.py = Simulation file:  Run "python3 Monopole2.py" to execute Fenics finite element solver\
 
 DualMonopole* - Set of files for simulating monopole linear array\
 *.geo = GMSH generator command file
 *.msh = Finite element mesh file generated by GMSH\
 DualMonopoleConv.py = Mesh converter\
 DualMonopoleInPhase.py = Simulation where even excitation mode is generated\
 DualMonopoleRvePhase.py = Simulation where odd excitation mode is computed\

CoaxWGLaunch - Set of files for simulating coaxial to rectangular waveguide transition\
CoaxWaveguideTransition.geo = GMSH mesh generator command file\
CoaxWaveguideTransition.msh = GMSH generated mesh file\
CoaxWaveguideTransitionConv.py = Mesh converter routine\
CoaxWGLaunch.py = Fenics simulation of waveguide transition\
