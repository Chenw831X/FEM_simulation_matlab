%
% fem2d_linear_quad.m
%
% Created by Wei Chen on 11/7/21
%

function fem2d_linear_quad()
    disp('fem simulation start');

    % quadrilateral element size is a*b
    a = 0.1;
    b = 0.1;
    % prepare finite element analysis
    nelx = 10;
    nely = 10;
    nele = nelx * nely;
    ndof = 2 * (nelx+1) * (nely+1);
    F = zeros(ndof, 1);
    U = zeros(ndof, 1);

    % user-defined load dofs
    disp('compute load');
    loadnid = nelx * (nely + 1) + 1;
    loaddof = 2 * loadnid - 1;
    F(loaddof, 1) = 0.01;

    % user-defined support fixed dofs
    disp('compute DBC');
    fixednid = (nely+1):(nely+1):(nelx+1)*(nely+1);
    fixeddof = [2*fixednid(:); 2*fixednid(:)-1];
    freedof = setdiff(1:ndof, fixeddof);
    
    nodegrd = reshape(1:(nely+1)*(nelx+1), nely+1, nelx+1);
    nodeids = reshape(nodegrd(1:end-1, 1:end-1), nely*nelx, 1);
    edofVec = 2*nodeids(:) + 1;
    edofMat = repmat(edofVec, 1, 8) + ...
        repmat([0 1 2*nely+[2 3 0 1] -2 -1], nele, 1);
    iK = reshape(kron(edofMat, ones(8, 1))', 8*8*nele, 1);
    jK = reshape(kron(edofMat, ones(1, 8))', 8*8*nele, 1);
    KE = lk;

    % FE-ANALYSIS
    sK = repmat(KE(:), nele, 1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    disp('solve');
    U(freedof, :) = K(freedof, freedof) \ F(freedof, :);

    % display result
    cor = zeros((nelx+1)*(nely+1), 2);
    for i = 0:nelx
        ii = i * (nely + 1);
        for j = 0:nely
            tmp = ii + nely + 1 -j;
            cor(tmp, 1) = i*a;
            cor(tmp, 2) = j*b;
        end
    end
    cor1 = cor + reshape(U, 2, (nelx+1)*(nely+1))';
    figure;
    scatter(cor1(:, 1), cor1(:, 2));

    disp('simulation end');
end

function [KE]=lk
    E = 1.; 
    nu = 0.3;
    k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
    -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
    KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end