%
% GenerateMesh.m
%
% Created by Wei Chen on 8/27/21
%

% @Function:
%   To mesh a rectangle membrane to use in FEA
% @Input:
%   lx: length of membrane along X-axis
%   ly: length of membrane along Y-axis
%   elex: number of elements along X-axis
%   eley: number of elements along Y-axis
% @Output:
%   nodes: the nodal coordinates of the mesh
%   eles: the nodal index of the elements
%   eleNum: number of elements
%   nodeNum: number of nodes
%    _______
%   | 4   3 |
%   |       |
%   |_1___2_|
function [nodes, eles, eleNum, nodeNum] = GenerateMesh(lx, ly, elex, eley)
    eleNum = elex * eley;
    nx = elex + 1;
    ny = eley + 1;
    nodeNum = nx * ny;

    x = linspace(0, lx, nx);
    y = linspace(0, ly, ny);
    [xx, yy] = meshgrid(x, fliplr(y));
    nodes = [xx(:), yy(:)];

    eles = zeros(eleNum, 4);
    nodeNo = reshape(1:nodeNum, ny, nx);
    eles(:, 4) = reshape(nodeNo(1:ny-1, 1:nx-1), eleNum, 1);
    eles(:, 3) = reshape(nodeNo(1:ny-1, 2:nx), eleNum, 1);
    eles(:, 2) = reshape(nodeNo(2:ny, 2:nx), eleNum, 1);
    eles(:, 1) = reshape(nodeNo(2:ny, 1:nx-1), eleNum, 1);
end