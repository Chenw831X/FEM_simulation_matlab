%
% PK1.m
%
% Created by Wei Chen on 9/5/21
%

function P = PK1(F, u, l)
% compute PK1 stress tensor
%
% Syntax: P = PK1(F, u, l)
%
% @Input:
%   F: deformation gradient, (3, 3)
%   u, l: Lame coefficients
% @Output:
%   P: PK1 stress tensor, (3, 3)
    JJ = log(det(F));
    Finv = inv(F);
    FinvT = Finv';

    P = u * (F - FinvT) + l * JJ * FinvT;
end