%
% dNdX.m
%
% Created by Wei Chen on 8/27/21
%

% @Function:
%   compute dNdX of 4 Gauss nodes
%   where X is reference coordinates, x is deformed coordinates
% @Input:
%   a, b is the length along x, y axes of each element
% @Output:
%   dN: cell(4, 1), dN{i, 1} is 4*2
function dN = dNdX(a, b)
    GaussNodes = [-1/sqrt(3); 1/sqrt(3)];
    dN = cell(4, 1);

    for i = 1:2
        for j = 1:2
            GN_x = GaussNodes(i);
            GN_y = GaussNodes(j);
            
            dN_x = 1/4 * [-(1-GN_y)  (1-GN_y)  (1+GN_y) -(1+GN_y)];
            dN_y = 1/4 * [-(1-GN_x) -(1+GN_x)  (1+GN_x)  (1-GN_x)];

            J = [dN_x; dN_y] * [ -a  a  a -a;
                                 -b -b  b  b]';
            B = J \ [dN_x; dN_y]; 
            dN{2*(i-1)+j, 1} = B';
        end
    end
end