%
% designDomain.m
%
% Created by Wei Chen on 8/28/21
%

% @Function:
%   set boundary conditions
% @Input:
%   opt: mode of design domain
%   nelx: element number along X-axis
%   nely: element number along Y-axis
%   vF: force value
% @Output:
%   freedofs: dofs not in DBC
%   fixeddofs: dofs in DBC
%   Fext: load in NBC
function [freedofs, fixeddofs, Fext] = designDomain(opt, nelx, nely, vF)
    ndofs = 2 * (nelx + 1) * (nely + 1);

    switch opt
    case 'cantilever'
        fixeddofs = 1:2*(nely+1);
        loadnid = 2*(nely+1)*(nelx+1);
        vF = -vF;
    end

    freedofs = setdiff(1:ndofs, fixeddofs);
    Fext = sparse(loadnid, 1, vF, ndofs, 1);
end