// 2D channel with a circular obstacle (DFG 2D benchmark / paper Sec. 5)
//   Omega = (0,2.2) x (0,0.41) \ B_r(0.2,0.2),  r = 0.05
// Same geometry and physical tags as redbKIT's Problems/FEM_CFD_Dfg2D/mesh/obstacleDFG.geo
// Mesh sizes are passed from the command line:
//   gmsh cylinder.geo -2 -format msh22 -setnumber s1 <h_wall> -setnumber s3 <h_cyl> -o out.msh
//   LF (coarse): s1 = 0.02,  s3 = 0.005
//   HF (fine)  : s1 = 0.015, s3 = 0.002
If (!Exists(s1)) s1 = 0.02;  EndIf    // size at the 4 channel corners
If (!Exists(s3)) s3 = 0.005; EndIf    // size on the cylinder

L  = 2.2;  H = 0.41;
xc = 0.2;  yc = 0.2;  r = 0.05;

Point(1) = {0, 0, 0, s1};
Point(2) = {L, 0, 0, s1};
Point(3) = {L, H, 0, s1};
Point(4) = {0, H, 0, s1};

Point(5) = {xc,   yc,   0, 1.0};      // centre (not meshed)
Point(6) = {xc+r, yc,   0, s3};
Point(7) = {xc,   yc+r, 0, s3};
Point(8) = {xc-r, yc,   0, s3};
Point(9) = {xc,   yc-r, 0, s3};

Circle(1) = {6, 5, 7};
Circle(2) = {7, 5, 8};
Circle(3) = {8, 5, 9};
Circle(4) = {9, 5, 6};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9)  = {7, 8, 5, 6};
Line Loop(10) = {1, 2, 3, 4};
Plane Surface(11) = {9, 10};

// Boundary flags used by the solver
Physical Line(1) = {8};            // inlet  x = 0      (Gamma_D2, parabolic inflow)
Physical Line(2) = {6};            // outlet x = 2.2    (Gamma_N, sigma n = 0)
Physical Line(3) = {7, 5};         // walls  y = 0.41, y = 0 (Gamma_D1, no-slip)
Physical Line(4) = {1, 2, 3, 4};   // cylinder          (Gamma_D1, no-slip; drag/lift)
Physical Surface(1) = {11};
