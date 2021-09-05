%
% dPdF.m
%
% Created by Wei Chen on 9/5/21
%

function dP = dPdF(F, dF, u, l)
% compute PK1 stress differentials for Neohookean materials
%
% Syntax: dP = dPdF(F, dF, u, l)
%
% @Input:
%   F: deformation gradient, (3, 3)
%   dF: differentials of F, (3, 3)
%   u, l: Lame coefficients
% @Output:
%   dP: PK1 stress differentials, (3, 3)
    JJ = log(det(F));
    Finv = inv(F);
    FinvT = Finv';

    dP = u * dF;
    dP = dP + (u - l * JJ) * FinvT * dF' * FinvT;
    dP = dP + l * trace(F \ dF) * FinvT;
end