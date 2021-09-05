%
% boundaryK.m
%
% Created by Wei Chen on 9/5/21
%

function K = boundaryK(K, DBCDofs)
% boundary K based on DBCDofs
%
% Syntax: K = boundaryK(K, DBCDofs)
%
% @Input:
%   K: stiffness matrix
%   DBCDofs: dofs in DBC
% @Output:
%   K: stiffness matrix after boundary
    K(DBCDofs, :) = 0;
    K(:, DBCDofs) = 0;

    for i = 1:size(DBCDofs, 2)
        dofI = DBCDofs(1, i);
        K(dofI, dofI) = 1;
    end
end