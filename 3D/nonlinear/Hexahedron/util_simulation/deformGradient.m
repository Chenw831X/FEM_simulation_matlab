%
% deformGradient.m
%
% Created by Wei Chen on 9/5/21
%

function F = deformGradient(ele, vet, vet0, dN)
% compute deformation gradient of 8 Gauss points of each element
%
% Syntax: F = deformGradient(ele, vet, vet0, dN)
%
% @Input:
%   ele: node ids of each element
%   vet: physical coordinates
%   vet0: reference coordinates
%   dN: dNdX
% @Output:
%   F: deformation gradient, cell(8, 1), F{i, 1} is (3, 3, nele)
%      F{i, 1}(:, :, eleI) is deformation gradient of i-th Gauss nodes of eleI-th element
    nele = size(ele, 1);
    F = cell(8, 1);

    for gp = 1:8
        Fgp = zeros(3, 3, nele);
        B = dN{gp, 1}; % (8, 3)
        for eleI = 1:nele
            edof = ele(eleI, :);
            x0 = vet0(edof, :);
            x = vet(edof, :);
            u = x - x0; % (8, 3)
            Fgp(:, :, eleI) = eye(3) + u.' * B;
        end
        F{gp, 1} = Fgp;
    end
end