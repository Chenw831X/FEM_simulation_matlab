%
% fem3d_nonlinear_hex.m
%
% Created by Wei Chen on 9/3/21
%

addpath util util_simulation;

% user-defined material properties
E = 1.0;
nu = 0.3;
Lv = 0.01; % load value
nelx = 10;
nely = 10;
nelz = 10;
dx = 0.1;

% initialize
u = E / 2 / (1 + nu);
l = E * nu / (1 + nu) / (1 - 2 * nu);

[vet0, ele, nvet, nele] = GenerateMesh(nelx*dx, nely*dx, nelz*dx, nelx, nely, nelz);
vet = vet0;
[freeDofs, DBCDofs, Fext] = designDomain(nelx, nely, nelz, Lv);
[edofMat, iK, jK] = forAssembly(ele);

disp('fem simulation start');