%
% deformGradient.m
%
% Created by Wei Chen on 8/28/31
%

% @Function:
%   compute deformation gradient of 4 Gauss Points in each element
% @Input:
%   eles: nodal index of each elements
%   deform_nodes: current nodal coordinates
%   ref_nodes: reference nodal coordiantes
%   dN: dNdX
% @Output:
%   F: cell(4, 1), F{i, 1} is 2*2*eleNum, 
%      F{i, 1}(:, :, ele) is deformation gradient of i-th Gauss Points in ele-th element
function F = deformGradient(eles, deform_nodes, ref_nodes, dN)
    eleNum = size(eles, 1);
    F = cell(4, 1);
    for gp = 1:4
        B = dN{gp, 1};
        Fm = zeros(2, 2, eleNum);

        for ele = 1:eleNum
            nodeIdx = eles(ele, :);
            x = deform_nodes(nodeIdx, :);
            X = ref_nodes(nodeIdx, :);
            Fm(:, :, ele) = eye(2) + (x - X)' * B;
        end
        F{gp, 1} = Fm;
    end
end