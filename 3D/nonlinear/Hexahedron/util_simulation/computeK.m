%
% computeK.m
%
% Created by Wei Chen on 9/5/21
%

function K = computeK(F, dN, iK, jK, dx, u, l)
% compute stiffness matrix
%
% Syntax: K = computeK(F, dN, iK, jK, dx, u, l)
%
% @Input:
%   F: deformation gradient, cell(8, 1), F{i, 1} is (3, 3, nele)
%   dN: dNdX, cell(8, 1), dN{i, 1} is (8, 3)
%   iK, jK: row-index/col-index of each element stiffness matrix
%   dx: size of hexahedron element
%   u, l: Lame coefficient
% @Output:
%   K: stiffness matrix
    nele = size(F{1, 1}, 3);
    sK = zeros(24*24, nele);
    jac = dx * dx * dx / 8;
    dFe = zeros(3, 3);

    for eleI = 1:nele
        Ke = zeros(24, 24);
        for gp = 1:8
            Fe = F{gp, 1}(:, :, eleI); % (3, 3)
            B = dN{gp, 1}; % (8, 3)

            for vI = 1:8
                for dof = 1:3
                    col = (vI - 1) * 3 + dof;
                    dFe(:) = 0;
                    dFe(dof, :) = B(vI, :);
                    dP = dPdF(Fe, dFe, u, l); % (3, 3)

                    tmp = dP * B';
                    Ke(:, col) = Ke(:, col) + tmp(:) * jac;
                end
            end
        end
        sK(:, eleI) = Ke(:);
    end

    K = sparse(iK, jK, sK);
    K = (K + K') / 2;
end