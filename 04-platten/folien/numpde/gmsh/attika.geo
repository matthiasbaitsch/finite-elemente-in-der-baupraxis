SetFactory("OpenCASCADE");

// rechtecke
Rectangle(1) = {0, 0, 0, 0.3, 1, 0};
Rectangle(2) = {-0.05, 0, 0, 0.05, 1, 0};
Rectangle(3) = {-0.05, 1, 0, 0.4, 0.05, 0};
Rectangle(4) = {0.3, 0.6, 0, 1, 0.05, 0};
Rectangle(5) = {0.3, 0.65, 0, 0.05, 0.35, 0};
Rectangle(6) = {0.3, 0.3, 0, 1, 0.3, 0};

// keil
Point(25) = {0.45, 0.65, 0, 1.0};
Point(26) = {0.35, 0.75, 0, 1.0};
Line(25) = {18, 25};
Line(26) = {25, 26};
Line(27) = {26, 18};
Curve Loop(7) = {25, 26, 27};
Plane Surface(7) = {7};

// verschmelzen
Coherence;

// namen
Physical Curve("wa") = {37, 41, 40, 39, 46, 26, 44};
Physical Curve("wi") = {29, 47};
Physical Surface("sd") = {2, 3, 5, 7, 4};
Physical Surface("sc") = {1, 6};

// netz
MeshSize{:} = 0.00625;
Mesh 2;
Save "attika.msh";


