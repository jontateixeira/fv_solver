cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 0.05, 0, cl1};
Point(4) = {0, 0.05, 0, cl1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(8) = {2, 3, 4, 1};
Plane Surface(8) = {8};

// MATERIAL/REGIONS FLAG: < 100
// DIRICHLET BOUNDARY CONDITION FLAG: > 100
// NEUMANN BOUNDARY CONDITION FLAG: > 200

Physical Point(202) = {1, 4};
Physical Point(201) = {2, 3};
Physical Line(201) = {1, 3};
Physical Line(101) = {2};
Physical Line(202) = {4};
Physical Surface(1) = {8};

Transfinite Line {1,3} = 21 Using Progression 1.0;
Transfinite Line {2,4} = 6 Using Progression 1.0;
Transfinite Surface {8} = {1,2,3,4};

Recombine Surface {8};
//+
View "comments" {
  T2(10,-10,0){"Created by JCT © 2020"};
};
Mesh.MshFileVersion = 2.2;
