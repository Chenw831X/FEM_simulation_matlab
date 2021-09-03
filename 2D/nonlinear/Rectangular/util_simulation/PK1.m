%
% PK1.m
%
% Created by Wei Chen on 8/28/21
%

% @Function:
%   compute PK1 stress of Neo-Hookean material
% @Input:
%   F: 2*2, deformation gradient
%   u, l: Lame coefficient
% @Output:
%   P: 2*2, PK1 stress
function P = PK1(F, u, l)
    JJ = log(det(F));
    Finv = inv(F);
    FinvT = Finv';

    P = u * (F - FinvT) + l * JJ * FinvT;
end