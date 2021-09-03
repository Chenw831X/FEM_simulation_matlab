%
% GenerateMesh.m
%
% Created by Wei Chen on 9/3/21
%

function [vet, ele, nvet, nele] = GenerateMesh(lx, ly, lz, nelx, nely, nelz)
% To mesh a hexahedron cube to use in FEA
%
% Syntax: [vet, ele, nvet, nele] = GenerateMesh(lx, ly, lz, nelx, nely, nelz)
%
% @Input:
%   lx, ly, lz: length of mesh along x, y, z axis
%   nelx, nely, nelz: number of element along x, y, z axis
% @Output:
%   vet: coordinates of the mesh
%   ele: nodal index of the elements
%   nvet: number of vertices
%   nele: number of elements
    nele = nelx * nely * nelz;
    nx = nelx + 1;
    ny = nely + 1;
    nz = nelz + 1;
    nvet = nx * ny * nz;

    x = linspace(0, lx, nx);
    y = linspace(0, ly, ny);
    z = linspace(0, lz, nz);
    [xx, yy] = meshgrid(x, fliplr(y));
    xx = repmat(xx(:), nz, 1);
    yy = repmat(yy(:), nz, 1);
    zz = kron(z(:), ones(nx*ny, 1));
    vet = [xx, yy, zz];

    vetgrd = reshape(1:nx*ny, ny, nx);
    vetids = reshape(vetgrd(1:end-1, 1:end-1), nelx*nely, 1);
    vetidz = 0:nx*ny:(nelz-1)*nx*ny;
    vetids = repmat(vetids, 1, nelz) + repmat(vetidz, nelx*nely, 1);
    ele = vetids(:) + 1;
    ele = repmat(ele, 1, 8) + repmat([0 ny ny-1 -1 nx*ny+[0 ny ny-1 -1]], nele, 1);
end