%
% designDomain.m
%
% Created by Wei Chen on 9/4/21
%

function [freeDofs, DBCDofs, Fext] = designDomain(nelx, nely, nelz, Lv)
% set boundary condition
%
% Syntax: [freeDofs, DBCDofs, Fext] = designDomain(nelx, nely, nelz, Lv)
%
% @Input:
%   nelx, nely, nelz: number of elements along x, y, z axis
%   Lv: load value
% @Output:
%   freeDofs: free dofs
%   FixedDofs: fixed dofs (DBC)
%   Fext: load force (NBC)
    nx = 1 + nelx;
    ny = 1 + nely;
    nz = 1 + nelz;
    ndofs = 3 * nx * ny * nz;
    [Dx, Dy, Dz] = meshgrid(0:nelx, 0, 0:nelz);
    DBCVInds = (Dz(:) * nx * ny) + (Dx(:) * ny) + (nely + 1 - Dy(:));
    DBCDofs = [DBCVInds*3-2; DBCVInds*3-1; DBCVInds*3]';
    freeDofs = setdiff(1:ndofs, DBCDofs);

    [Nx, Ny, Nz] = meshgrid(0, nely, 0:nelz);
    NBCVInds = (Nz(:) * nx * ny) + (Nx(:) * ny) + (nely + 1 -Ny(:));
    NBCDofs = NBCVInds * 3 - 2;
    Fext = sparse(NBCDofs, 1, -Lv, ndofs, 1);
end