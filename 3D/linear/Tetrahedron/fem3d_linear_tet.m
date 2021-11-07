%
% fem3d_tet_linear.m
% 
% Created by Wei Chen on 8/4/21
%

function fem3d_linear_tet() % nvet: number of vertices, nele: number of elements
nvet = 3168;
nele = 16563;

disp('fem simulation start')

eeps = 1e-8;
% USER-DEFINED MATERIAL PROPERTIES
E = 1.0;          % Young's modulus of solid material
nu = 0.3;         % Poisson's ratio 
% constitutive matrix
D = E / (1 + nu) / (1 - 2 * nu) * ...
[ 1-nu   nu   nu     0          0          0     ;
    nu 1-nu   nu     0          0          0     ;
    nu   nu 1-nu     0          0          0     ;
     0    0    0 (1-2*nu)/2     0          0     ;
     0    0    0     0      (1-2*nu)/2     0     ;
     0    0    0     0          0      (1-2*nu)/2];

% read coordinates of vertices, tetrahedron information
cor = load('input/cor.txt');
elem = load('input/elem.txt');

% PREPARE FINITE ELEMENT ANALYSIS
ndof = nvet * 3;
U = zeros(ndof,1);
F = zeros(ndof,1);

% USER-DEFINED LOAD DOFs
for vI = 1:nvet
    if (abs(cor(vI, 3) - 1.0) < eeps) && (abs(cor(vI, 1)) < eeps)
        F(vI*3-2, 1) = -0.01; 
    end
end
disp('compute load');

% USER-DEFINED SUPPORT FIXED DOFs
fixednids = find(abs(cor(:, 3)) < eeps);
fixeddofs = [3*fixednids-2; 3*fixednids-1; 3*fixednids];
freedofs = setdiff(1:ndof,fixeddofs);
disp('compute DBC')

% assembly the stiffness matrix
edofMat = kron(elem, 3*ones(1, 3));
edofMat(:, 1:3:12) = edofMat(:, 1:3:12) - 2;
edofMat(:, 2:3:12) = edofMat(:, 2:3:12) - 1;
iK = reshape(kron(edofMat,ones(12,1))',12*12*nele,1);
jK = reshape(kron(edofMat,ones(1,12))',12*12*nele,1);

sK = zeros(12*12, nele);
for eleI = 1:nele
    X = cor([elem(eleI, 1), elem(eleI, 2), elem(eleI, 3), elem(eleI, 4)], :);
    KE = initKE(D, X);
    sK(:, eleI) = KE(:);
end
sK = reshape(sK, 12*12*nele, 1);

% FE-ANALYSIS
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
disp('solve')

save U.mat U;
% display result
cor1 = cor + reshape(U, 3, nvet)';
scatter3(cor1(:, 1), cor1(:, 2), cor1(:, 3));

disp('simulation end')

end

% === GENERATE ELEMENT STIFFNESS MATRIX ===
% X: 4*3 matrix, which is corridinates of element vertices
function [KE] = initKE(D, X)
H = ...
[ 1.0 X(1, 1) X(1, 2) X(1, 3);
  1.0 X(2, 1) X(2, 2) X(2, 3);
  1.0 X(3, 1) X(3, 2) X(3, 3);
  1.0 X(4, 1) X(4, 2) X(4, 3)];
V6 = abs(det(H)); % 6*V, V is volumn of tetrahedron element

% a1 =  det(H([2, 3, 4], [2, 3, 4]));
% a2 = -det(H([1, 3, 4], [2, 3, 4]));
% a3 =  det(H([1, 2, 4], [2, 3, 4]));
% a4 = -det(H([1, 2, 3], [2, 3, 4]));

b1 = -det(H([2, 3, 4], [1, 3, 4]));
b2 =  det(H([1, 3, 4], [1, 3, 4]));
b3 = -det(H([1, 2, 4], [1, 3, 4]));
b4 =  det(H([1, 2, 3], [1, 3, 4]));

c1 =  det(H([2, 3, 4], [1, 2, 4]));
c2 = -det(H([1, 3, 4], [1, 2, 4]));
c3 =  det(H([1, 2, 4], [1, 2, 4]));
c4 = -det(H([1, 2, 3], [1, 2, 4]));

d1 = -det(H([2, 3, 4], [1, 2, 3]));
d2 =  det(H([1, 3, 4], [1, 2, 3]));
d3 = -det(H([1, 2, 4], [1, 2, 3]));
d4 =  det(H([1, 2, 3], [1, 2, 3]));

B = zeros(6, 12);

B(:, 1:3) = ...
[ b1  0  0;
   0 c1  0;
   0  0 d1;
  c1 b1  0;
   0 d1 c1;
  d1  0 b1];

B(:, 4:6) = ...
[ b2  0  0;
   0 c2  0;
   0  0 d2;
  c2 b2  0;
   0 d2 c2;
  d2  0 b2];

B(:, 7:9) = ...
[ b3  0  0;
   0 c3  0;
   0  0 d3;
  c3 b3  0;
   0 d3 c3;
  d3  0 b3];

B(:, 10:12) = ...
[ b4  0  0;
   0 c4  0;
   0  0 d4;
  c4 b4  0;
   0 d4 c4;
  d4  0 b4];

B = B / V6;
KE = V6 / 6.0 * (B' * D * B);
end