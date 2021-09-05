%
% linesearch.m
%
% Created by Wei Chen on 9/5/21
%

function [vetNew, converged] = linesearch(p, dN, rhs, Fext, vet, vet0, ele, ...
                                          edofMat, DBCDofs, dx, u, l)
% line search to compute step
%
% Syntax: [vetNew, converged] = linesearch(p, dN, rhs, Fext, vet, vet0, ele, ...
%                                          edofMat, DBCDofs, dx, u, l)
    nvet = size(vet, 1);
    deform0 = reshape(vet', [], 1);
    ref = reshape(vet0', [], 1);
    approxE0 = -sum(dot(rhs, deform0-ref));
    norm0 = max(max(abs(rhs)));
    step = 1;
    foundStep = 0;

    maxloop = 10;
    loop = 0;
    while loop < maxloop
        loop = loop + 1;

        deform1 = deform0 + step * p;
        F = deformGradient(ele, reshape(deform1, 3, [])', vet0, dN);

        inverted = checkInversion(F);
        if (inverted)
            step = step * 0.5;
            continue;
        end

        Fint = elasticForce(F, dN, nvet, edofMat, dx, u, l);
        rhs = Fint + Fext;
        rhs(DBCDofs) = 0;

        approxE1 = -sum(dot(rhs, deform1-ref));
        norm1 = max(max(abs(rhs)));
        if ((norm1 < norm0) || (approxE1 < approxE0))
            foundStep = 1;
            break;
        end
        step = step * 0.5;
    end
    
    converged = 0;
    if(~foundStep || abs(norm1 - norm0) < 1e-5)
        converged = 1;
    end

    deform1(DBCDofs) = deform0(DBCDofs);
    vetNew = reshape(deform1, 3, [])';
end

function inverted = checkInversion(F)
    inverted = false;
    nele = size(F{1, 1}, 3);
    eeps = 1e-5;
    for gp = 1:8
        Fgp = F{gp, 1}; % (3, 3, nele)
        for eleI = 1:nele
            Fe = Fgp(:, :, eleI);
            if (det(Fe) < eeps)
                inverted = true;
                break;
            end
        end
        if (inverted)
            break;
        end
    end
end