%
% fem2d_hex_nonlinear.m
%
% Created by Wei Chen on 8/27/21
%

addpath util util_simulation

% user-defined material properties
E = 1e3;
nu = 0.3;
optDesign = 'cantilever';
vF = 100; % load value
nelx = 20;
nely = 10;
dx = 1;

% initialize
u = E / 2 / (1 + nu);
l = E * nu / (1 + nu) / (1 - 2 * nu);

disp('fem simulation start');
[ref_nodes, eles, eleNum, nodeNum] = GenerateMesh(nelx*dx, nely*dx, nelx, nely);
dN = dNdX(dx/2, dx/2); % dNdX, cell(4, 1), dN{i, 1} is 4*2

deform_nodes = ref_nodes;
[freedofs, fixeddofs, Fext] = designDomain(optDesign, nelx, nely, vF);
[iK, jK, edofMat] = forAssembly(nelx, nely);

disp('Newton method start');
maxloop = 10;
loop = 0;
converged = 0;
while loop < maxloop && ~converged
    loop = loop + 1;
    % compute deformation gradient
    F = deformGradient(eles, deform_nodes, ref_nodes, dN); % cell(4, 1), F{i, 1} is 2*2*eleNum
    % compute stiffness matrix
    K = computeK(F, dN, iK, jK, dx, u, l);
    % boundary stiffness matrix
    K = boundaryK(K, fixeddofs);

    % compute elastic force
    Fint = elasticForce(F, dN, nodeNum, edofMat, dx, u, l);
    rhs = Fint + Fext;
    rhs(fixeddofs) = 0;
    % compute search direction
    p = K \ rhs;

    % line search
    [deform_nodes, converged] = linesearch(p, dN, rhs, Fext, deform_nodes, ref_nodes, ...
                                            eles, edofMat, fixeddofs, dx, u, l);
    fprintf("Iteration %i\n", loop);
end
figure;
scatter(deform_nodes(:, 1), deform_nodes(:, 2), 'filled');
axis equal;