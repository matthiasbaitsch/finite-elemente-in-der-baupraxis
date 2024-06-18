SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 0.015, 0.1, 0};
Rectangle(2) = {0.015, 0, 0, 0.3, 0.1, 0};
Rectangle(3) = {0.315, 0, 0, 0.02, 0.1, 0};
Coherence;

Physical Curve("wi") = {4};
Physical Curve("wa") = {9};
Physical Surface("s1") = {1};
Physical Surface("s2") = {2};
Physical Surface("s3") = {3};

MeshSize{:} = 0.015/2;
Mesh 2;

Save "wall.msh";
