// structured quadrilateral mesh
n = 50;

Point(1) = {-1, -1, 0, 2.0/n};
Point(2) = {1, -1, 0, 2.0/n};
Point(3) = {1, 1, 0, 2.0/n};
Point(5) = {-1, 1, 0, 2.0/n};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 1};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};
Physical Curve(1) = {4, 3, 2, 1};
Physical Surface(2) = {1};

// Periodic Curve{2} = {-4}; // periodic along x
// Periodic Curve{3} = {-1}; // periodic along y

Recombine Surface {1};
Transfinite Surface {1};
