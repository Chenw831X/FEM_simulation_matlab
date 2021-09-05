%
% dNdX.m
%
% Created by Wei Chen on 9/5/21
%

function dN = dNdX(a, b, c)
% compute dNdX of 8 Gauss nodes
% where X is the reference coordinates, x is the deformed coordinates
%
% Syntax dN = dNdX(a, b, c)
%
% @Input:
%   a, b, c: element size along x, y, z axis
% @Output
%   dN: cell(8, 1), dN{i, 1} is (8, 3)
    GN_x = [-1/sqrt(3); 1/sqrt(3)];
    GN_y = GN_x;
    GN_z = GN_x;
    dN = cell(8, 1);

    for i = 1:length(GN_x)
        for j = 1:length(GN_y)
            for k = 1:length(GN_z)
                x = GN_x(i);
                y = GN_y(j);
                z = GN_z(k);
                dNx = 1/8 * [-(1-y)*(1-z)  (1-y)*(1-z)  (1+y)*(1-z) -(1+y)*(1-z) ...
                             -(1-y)*(1+z)  (1-y)*(1+z)  (1+y)*(1+z) -(1+y)*(1+z)];
                dNy = 1/8 * [-(1-x)*(1-z) -(1+x)*(1-z)  (1+x)*(1-z)  (1-x)*(1-z) ...
                             -(1-x)*(1+z) -(1+x)*(1+z)  (1+x)*(1+z)  (1-x)*(1+z)];
                dNz = 1/8 * [-(1-x)*(1-y) -(1+x)*(1-y) -(1+x)*(1+y) -(1-x)*(1+y) ...
                              (1-x)*(1-y)  (1+x)*(1-y)  (1+x)*(1+y)  (1-x)*(1+y)];
                J = [dNx; dNy; dNz] * [-a a a -a -a a a -a; -b -b b b -b -b b b; ...
                                       -c -c -c -c c c c c].';
                B = J \ [dNx; dNy; dNz];
                dN{(i-1)*4+(j-1)*2+k, 1} = B.';
            end
        end
    end
end