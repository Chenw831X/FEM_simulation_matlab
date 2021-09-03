%
% fem3d_hex_linear.m
% 
% Created by Wei Chen on 8/4/21
%

function fem3d_hex_linear()
disp(['fem simulation start'])

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

% hexahedron element size is a*b*c
a = 0.1;
b = 0.1;
c = 0.1;
% PREPARE FINITE ELEMENT ANALYSIS
nelx = 10;
nely = 10;
nelz = 10;
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
U = zeros(ndof,1);

% USER-DEFINED LOAD DOFs
disp(['compute load']);
[il,jl,kl] = meshgrid(0, 0:nely, nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 2;                             % DOFs
F = sparse(loaddof, 1, -0.01, ndof, 1);

% USER-DEFINED SUPPORT FIXED DOFs
disp(['compute DBC'])
[iif,jf,kf] = meshgrid(0:nelx,0:nely,0);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
freedofs = setdiff(1:ndof,fixeddof);

nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
KE = intKE(a/2, b/2, c/2, D);

% FE-ANALYSIS
sK = repmat(KE(:), nele, 1);
K = sparse(iK,jK,sK); K = (K+K')/2;
disp(['solve'])
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);

save U.mat U;
% display result
cor = zeros((nelx+1)*(nely+1)*(nelz+1), 3);
for k = 0:nelz
    kk = k * (nelx + 1) * (nely + 1);
    for i = 0:nelx
        ii = i * (nely + 1);
        for j = 0:nely
            tmp = kk + ii + nely + 1 - j;
            cor(tmp, 1) = i*a;
            cor(tmp, 2) = j*b;
            cor(tmp, 3) = k*c;
        end
    end
end
cor1 = cor + reshape(U, 3, (nelx+1)*(nely+1)*(nelz+1))';
figure;
scatter3(cor1(:, 1), cor1(:, 2), cor1(:, 3));

disp(['simulation end'])

end


% generate element stiffness matrix
% element size is a*b*c
function Ke = intKE(a, b, c, DH)
GN_x=[-1/sqrt(3),1/sqrt(3)]; GN_y=GN_x; GN_z=GN_x; GaussWeigh=[1,1];
Ke = zeros(24,24); L = zeros(6,9);
L(1,1) = 1; L(2,5) = 1; L(3,9) = 1;
L(4,2) = 1; L(4,4) = 1; L(5,6) = 1;
L(5,8) = 1; L(6,3) = 1; L(6,7) = 1;
for i=1:length(GN_x)
    for j=1:length(GN_y)
        for k=1:length(GN_z)
            x = GN_x(i);y = GN_y(j);z = GN_z(k);
            dNx = 1/8*[-(1-y)*(1-z)  (1-y)*(1-z)  (1+y)*(1-z) -(1+y)*(1-z) -(1-y)*(1+z)  (1-y)*(1+z)  (1+y)*(1+z) -(1+y)*(1+z)];
            dNy = 1/8*[-(1-x)*(1-z) -(1+x)*(1-z)  (1+x)*(1-z)  (1-x)*(1-z) -(1-x)*(1+z) -(1+x)*(1+z)  (1+x)*(1+z)  (1-x)*(1+z)];
            dNz = 1/8*[-(1-x)*(1-y) -(1+x)*(1-y) -(1+x)*(1+y) -(1-x)*(1+y)  (1-x)*(1-y)  (1+x)*(1-y)  (1+x)*(1+y)  (1-x)*(1+y)];
            J = [dNx;dNy;dNz]*[ -a  a  a  -a  -a  a  a  -a ;  -b  -b  b  b  -b  -b  b  b; -c -c -c -c  c  c  c  c]';
            G = [inv(J) zeros(3) zeros(3);zeros(3) inv(J) zeros(3);zeros(3) zeros(3) inv(J)];
            dN(1,1:3:24) = dNx; dN(2,1:3:24) = dNy; dN(3,1:3:24) = dNz;
            dN(4,2:3:24) = dNx; dN(5,2:3:24) = dNy; dN(6,2:3:24) = dNz;
            dN(7,3:3:24) = dNx; dN(8,3:3:24) = dNy; dN(9,3:3:24) = dNz;
            Be = L*G*dN;
            Ke = Ke + GaussWeigh(i)*GaussWeigh(j)*GaussWeigh(k)*det(J)*(Be'*DH*Be);
        end
    end
end
end

% === GENERATE ELEMENT STIFFNESS MATRIX ===
% element size is 1*1*1
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
