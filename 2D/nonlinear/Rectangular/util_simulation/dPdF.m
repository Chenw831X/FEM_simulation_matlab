%
% dPdF.m
%
% Created by Wei Chen on 8/28/21
%

% @Function:
%   compute PK1 stress differentials of Neo-Hookean material
% @Input:
%   F: deformation gradient
%   dF: differential of F
%   u, l: Lame coefficients
% @Output:
%   dP: PK1 stress differentials
function dP = dPdF(F, dF, u, l)
    JJ = log(det(F));
    Finv = inv(F);
    FinvT = Finv';

    dP = u * dF;
    dP = dP + (u - l * JJ) * FinvT * dF' * FinvT;
    dP = dP + l * trace(F \ dF) * FinvT;
end