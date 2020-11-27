// Gmsh project created on Thu Sep 05 11:10:00 2019
// All dimensions in centimeters
lc = 0.5;
lcb = 0.05;
pi = 3.1415926536;
a = 2.286; // WG width
b = 1.016; // WG height
l = 4.0;   // WG length
h = 0.648;  // Probe height
r = 0.065; // Probe radius
w = 0.554;  // Probe distance from WG short
rc = 0.22; // Coax outer diameter
eps_c = 2.2; // Coax dielectric
l_c = 1.0; // Coax length
//+
SetFactory("OpenCASCADE"); 
Box(1) = {0, 0, 0, a, b, l};  //  WG geom
Cylinder(2) = {a/2, -l_c, w, 0.0, l_c, 0.0, rc, 2*Pi};
Cylinder(3) = {a/2, -l_c, w, 0.0, l_c + h, 0.0, r, 2*Pi};
//+
//+
BooleanUnion(4) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
BooleanDifference(5) = { Volume{4}; Delete; }{ Volume{3}; Delete; };
//+
Characteristic Length{ 7, 12} = lcb;
Characteristic Length{1, 2, 3, 4, 5, 6, 8, 9, 10, 11} = lc;

Physical Volume("WG_Launch") = {5};

Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.3;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.MinimumCirclePoints = 18;
Mesh.Algorithm3D = 4;


