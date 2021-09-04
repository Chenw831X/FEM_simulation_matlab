%
% forAssembly.m
%
% Created by Wei Chen on 9/4/21
%

function [edofMat, iK, jK] = forAssembly(ele)
% compute edofMat, iK, jK for stiffness matrix assembly
%
% Syntax: [edofMat, iK, jK] = forAssembly(ele)
%
% @Input:
%   ele: node ids of each element
% @Output:
%   edofMat: (nele, 24), 24 dofs of each element
%   iK, jK: (nele*24*24, 1), element stiffness matrix row-index/col-index of each element
    nele = size(ele, 1);
    edofMat = kron(3*ele, ones(1, 3));
    edofMat(:, 1:3:24) = edofMat(:, 1:3:24) - 2;
    edofMat(:, 2:3:24) = edofMat(:, 2:3:24) - 1;
    iK = reshape(kron(edofMat, ones(24, 1)).', nele*24*24, 1);
    jK = reshape(kron(edofMat, ones(1, 24)).', nele*24*24, 1);
end