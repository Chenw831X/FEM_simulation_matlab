%
% elasticForce.m
%
% Created by Wei Chen on 8/28/21
%

% @Function:
%   compute elastic force of Neo-Hookean material
% @Input:
%   F: deforamation gradient, cell(4, 1), F{i, 1} is 2*2*eleNum
%   dN: dNdX, cell(4, 1), dN{i, 1} is 4*2
%   nodeNum: number of nodes
%   edofMat: eleNum*8
%   dx: size of element
%   u, l: Lame cofficients
% @Output:
%   Fint: elastic force, (2*nodeNum, 1)
function Fint = elasticForce(F, dN, nodeNum, edofMat, dx, u, l)
    Fint = zeros(2*nodeNum, 1);
    eleNum = size(edofMat, 1);
    jac = dx * dx / 4;

    for ele = 1:eleNum
        edof = edofMat(ele, :);
        Finte = zeros(8, 1);
        
        for gp = 1:4
            Fe = F{gp, 1}(:, :, ele);
            Be = dN{gp, 1};
            P = PK1(Fe, u, l);

            tmp = P * Be';
            Finte = Finte + tmp(:) * jac;
        end
        Fint(edof, 1) = Fint(edof, 1) - Finte; % NOTE: MINUS
    end
end