%
% forAssembly.m
%
% Created by Wei Chen on 8/28/21
%

function [iK, jK, edofMat] = forAssembly(nelx, nely)
    m = nelx * nely;
    nodeId = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
    nodeVec = reshape(2*nodeId(1:end-1, 1:end-1)+1, m, 1);
    edofMat = repmat(nodeVec, 1, 8) + repmat([0 1 2*nely+[2 3 0 1] -2 -1], m, 1);
    iK = reshape(kron(edofMat, ones(8, 1)).', 64*m, 1);
    jK = reshape(kron(edofMat, ones(1, 8)).', 64*m, 1);
end