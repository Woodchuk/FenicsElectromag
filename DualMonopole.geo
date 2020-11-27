// Wideband patch model
h = 3.0; // monopole height (cm)
r = 10.0; // radiation sphere radius
roffset = -1.25; // Rad sphere offset
lc = 1.0; // Coax length
rc = 0.25; // Coax outer radius
cc = 0.107; // Coax cntr cond radius 50 ohm with air diel
s = 0.0; // coax probe  x offset
q = 1.75; // coax probe y offset

SetFactory("OpenCASCADE");
//Radiation sphere
Sphere(1) = {0, roffset, 0, r};
Box(2) = {0, roffset, 0, r,  r, r};
BooleanIntersection(3) = {Volume{1}; Delete;}{Volume{2}; Delete;};

// Coax center
Cylinder(8) = {s, q, -lc, 0, 0, lc+h, cc};
// Coax shield
Cylinder(9) = {s, q, -lc, 0, 0, lc, rc};
Box(10) = {s, q-rc, -lc, rc, 2*rc, lc};
BooleanIntersection(11) = {Volume{9}; Delete;}{Volume{10}; Delete;};

// Add shield to rad volume
BooleanUnion(12) = {Volume{3}; Delete;}{Volume{11}; Delete;};
// Delete patches, etc from volume
BooleanDifference(14) = {Volume{12}; Delete;}{Volume{8}; Delete;};
Characteristic Length{7, 8, 14} = 0.05;
Characteristic Length{5, 10, 11} = 0.05;
Physical Volume("FullVol") = {14};
Mesh.CharacteristicLengthMax = 0.7;
Mesh.CharacteristicLengthMin = 0.02;
Mesh.MinimumCirclePoints = 18;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.Algorithm3D = 4;


