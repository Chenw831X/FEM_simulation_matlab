%
% boundaryK.m
%
% Created by Wei Chen on 8/28/21
%

% @Function:
%   add DBC to K, K(fixeddofs) should be 1 in diag and 0 otherwise
% @Input:
%   K: stiffness matrix
%   fixeddofs: DBC, displacements are 0 in fixeddofs
% @Output:
%   K: stiffness matrix after boundary
function K = boundaryK(K, fixeddofs)
    K(fixeddofs, :) = 0;
    K(:, fixeddofs) = 0;

    for i = 1:size(fixeddofs, 2)
        dof = fixeddofs(1, i);
        K(dof, dof) = 1;
    end
end