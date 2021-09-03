%
% computeK.m
%
% Created by Wei Chen on 8/28/21
%

% @Function:
%   compute and assemble stiffness matrix
% @Input:
%   F: deformation gradient
%   dN: dNdX
%   iK, jK: i-index and j-index of nonzero element in stiffness matrix
%   dx: size of element
%   u, l: Lame coefficients
function K = computeK(F, dN, iK, jK, dx, u, l)
    eleNum = size(F{1, 1}, 3);
    jac = dx * dx / 4;
    dFe = zeros(2, 2);
    sK = zeros(64, eleNum);

    for ele = 1:eleNum
        Ke = zeros(8, 8);
        for gp = 1:4
            Fe = F{gp, 1}(:, :, ele);
            Be = dN{gp, 1};

            for nI = 1:4
                for dof = 1:2
                    col = 2 * (nI - 1) + dof;
                    dFe(:) = 0;
                    dFe(dof, :) = Be(nI, :);
                    dP = dPdF(Fe, dFe, u, l);
                    
                    tmp = dP * Be';
                    Ke(:, col) = Ke(:, col) + tmp(:) * jac;
                end
            end
        end
        sK(:, ele) = Ke(:);
    end

    K = sparse(iK, jK, sK);
    K = (K + K') / 2;
end