// Patch antenna model
h = 0.15; // substrate height
r = 2.0; // Ground plane radius
a = 0.60; // Monopole length
b = 0.8;
t = 0.1;
lc = 1.0; // length coax
rc = 0.11; // Coax shield radius
cc = 0.0325; // Coax center conductor radius
s = 0.0; // Coax probe offset
ndr = 0.05; // Refinement near patch edge

SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, r};
Box(2) = {-r, -r, 0, 2.0*r, 2.0*r, r};
BooleanIntersection(3) = {Volume{1}; Delete;}{Volume{2}; Delete;};
Cylinder(4) = {0, 0, -h, 0, 0, h, r}; // Substrate
Cylinder(5) = {0, s, -h-lc, 0, 0, lc, rc}; // Coax shield
Cylinder(6) = {0, s, -h-lc, 0, 0, h+lc+t/2, cc}; // coax center
Cylinder(7) = {0, s, 0, 0, 0, a, cc}; // Antenna element
//Box(7) = {-a/2, -b/2, 0, a, b, t}; // patch
BooleanUnion(8) = {Volume{6}; Delete;}{Volume{7}; Delete;}; //to be subtracted out
BooleanDifference(9) = {Volume{3}; Delete;}{Volume{8};};
BooleanDifference(10) = {Volume{4}; Delete;}{Volume{8};};
BooleanDifference(11) = {Volume{5}; Delete;}{Volume{8}; Delete;};
BooleanFragments{Volume{9, 10, 11}; Delete;}{}
//Characteristic Length{6, 7, 8, 9, 10, 11, 12, 13} = ndr;
Box(12) = {-r, -r, -h-lc, 2*r, r, r+h+lc}; // slice out half because of symmetry
BooleanUnion(14) = {Volume{9}; Delete;}{Volume{10}; Delete;};
BooleanUnion(15) = {Volume{14}; Delete;}{Volume{11}; Delete;};
BooleanDifference(13) = {Volume{15}; Delete;}{Volume{12}; Delete;};
Box(16) = {-r, -r, -h-lc,r,2*r,r+h+lc}; // Slice out another half to model 1/4 of geom
BooleanDifference(17) = {Volume{13}; Delete;}{Volume{16}; Delete;};
Physical Volume("Patch") = {17};
Characteristic Length{8, 9, 20} = 0.003;
Mesh.CharacteristicLengthMax = 0.20;
Mesh.MinimumCirclePoints = 25;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.Algorithm3D = 4;

