// Gmsh project created on Sat Nov 11 16:29:57 2006
Point(1) = {0,0,0,0.02};
Point(2) = {2,0,0,0.02};
Point(3) = {2,1,0,0.02};
Point(4) = {0,1,0,0.02};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};
Physical Line(3) = {3,1};
Physical Line(2) = {2,4};
Physical Surface(10) = {6};
