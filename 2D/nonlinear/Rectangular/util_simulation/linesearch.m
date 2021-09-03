%
% linesearch.m
%
% Created by Wei Chen on 8/28/21
%

function [deform_nodes1, converged] = linesearch(p, dN, rhs, Fext, deform_nodes0, ...
                                                 ref_nodes, eles, edofMat, fixeddofs, ...
                                                 dx, u, l)
    nodeNum = size(deform_nodes0, 1);
    deform0 = reshape(deform_nodes0', [], 1);
    ref = reshape(ref_nodes', [], 1);
    approxE0 = -sum(dot(rhs, deform0-ref));
    norm0 = max(max(abs(rhs)));
    step = 1;
    foundStep = 0;

    maxloop = 10;
    loop = 0;
    while loop < maxloop
        loop = loop + 1;

        deform1 = deform0 + step * p;
        F = deformGradient(eles, reshape(deform1, 2, [])', ref_nodes, dN);
        
        inverted = checkInversion(F);
        if (inverted)
            step = step * 0.5;
            continue;
        end

        Fint = elasticForce(F, dN, nodeNum, edofMat, dx, u, l);
        rhs = Fint + Fext;
        rhs(fixeddofs) = 0;

        approxE1 = -sum(dot(rhs, deform1 - ref));
        norm1 = max(max(abs(rhs)));
        if ((norm1 < norm0) || (approxE1 < approxE0))
            foundStep = 1;
            break;
        end
        step = step * 0.5;
    end

    converged = 0;
    if(~foundStep || (abs(norm1 - norm0) < 1e-5))
        converged = 1;
    end

    deform1(fixeddofs) = deform0(fixeddofs);
    deform_nodes1 = reshape(deform1, 2, [])';
end

function inverted = checkInversion(F)
    inverted =false;
    eleNum = size(F{1, 1}, 3);
    eeps = 1e-5;
    for gp = 1:4
        Fgp = F{gp, 1}; % 2*2*eleNum
        for ele = 1:eleNum
            Fe = Fgp(:, :, ele);
            if det(Fe) <= eeps
                inverted = true;
                break;
            end
        end
        if (inverted)
            break;
        end
    end
end