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
disp('fem simulation start');
u = E / 2 / (1 + nu);
l = E * nu / (1 + nu) / (1 - 2 * nu);

[vet0, ele, nvet, nele] = GenerateMesh(nelx*dx, nely*dx, nelz*dx, nelx, nely, nelz);
vet = vet0;
[freeDofs, DBCDofs, Fext] = designDomain(nelx, nely, nelz, Lv);
[edofMat, iK, jK] = forAssembly(ele);
dN = dNdX(dx/2, dx/2, dx/2); % dNdX, cell(8, 1), dNdX{i, 1} is (8, 3)

disp('Newton method start');
maxloop = 10;
loop = 0;
converged = 0;
while loop < maxloop && ~converged
    loop = loop + 1;
    % compute deformation gradient
    F = deformGradient(ele, vet, vet0, dN); % cell{8, 1}, F{i, 1} is (3, 3, nele)
    % compute stiffness matrix

    K = computeK(F, dN, iK, jK, dx, u, l);
    % boundary stiffness matrix
    K = boundaryK(K, DBCDofs);

    % compute elastic force
    Fint = elasticForce(F, dN, nvet, edofMat, dx, u, l);
    rhs = Fint + Fext;
    rhs(DBCDofs) = 0;
    % solve linear system to compute search direction
    p = K \ rhs;

    [vet, converged] = linesearch(p, dN, rhs, Fext, vet, vet0, ele, edofMat, ...
                                  DBCDofs, dx, u, l);
    fprintf("Iteration %i\n", loop);
end
figure;
scatter3(vet(:, 1), vet(:, 3), vet(:, 2), 'filled');
axis equal;
