%
% elasticForce.m
%
% Created by Wei Chen on 9/5/21
%

function Fint = elasticForce(F, dN, nvet, edofMat, dx, u, l)
% compute elastic force
%
% Syntax: Fint = elasticForce(F, dN, nvet, edofMat, dx, u, l)
%
% @Input:
%   F: deformation gradient, cell(8, 1), F{i, 1} is (3, 3, nele)
%   dN: dNdX, cell(8, 1), dN{i, 1} is (8, 3)
%   nvet: number of vertices
%   edofMat: 24 edofs in each element, (nele, 24)
%   dx: size of hexahedron
%   u, l: Lame coefficient
% @Output:
%   Fint: elastic force, (3*nvet, 1)
    Fint = zeros(3*nvet, 1);
    nele = size(edofMat, 1);
    jac = dx * dx * dx / 8;

    for eleI = 1:nele
        edof = edofMat(eleI, :);
        Finte = zeros(24, 1);
        for gp = 1:8
            Fe = F{gp, 1}(:, :, eleI); % (3, 3)
            B = dN{gp, 1}; % (8, 3)
            P = PK1(Fe, u, l); % (3, 3)
            tmp = P * B';
            Finte = Finte + tmp(:) * jac;
        end
        Fint(edof, 1) = Fint(edof, 1) - Finte; % NOTE: MINUS
    end
end