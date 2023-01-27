// Unstuctured mixed grid
n = 50;

Point(1) = {-1, -1, 0, 1.0};
Point(2) = {1, -1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {-1, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {4, 2};

Curve Loop(1) = {4, 1, -5};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 2, 3};
Plane Surface(2) = {2};

Transfinite Curve {4, 1, 2, 3} = n Using Progression 1;
Transfinite Curve {5} = n*1.5 Using Progression 1;

// Periodic Curve{2} = {-4}; // periodic along x
// Periodic Curve{3} = {-1}; // periodic along y

Physical Curve(1) = {4, 3, 2, 1};
Physical Surface(2) = {1, 2};
Recombine Surface {2};
