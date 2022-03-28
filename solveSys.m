function [u,R] = solveSys(vL,vR,uR,KG,Fext)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - vL      Free degree of freedom vector
%   - vR      Prescribed degree of freedom vector
%   - uR      Prescribed displacement vector
%   - KG      Global stiffness matrix [n_dof x n_dof]
%              KG(I,J) - Term in (I,J) position of global stiffness matrix
%   - Fext    Global force vector [n_dof x 1]
%              Fext(I) - Total external force acting on DOF I
%--------------------------------------------------------------------------
% It must provide as output:
%   - u       Global displacement vector [n_dof x 1]
%              u(I) - Total displacement on global DOF I
%   - R       Global reactions vector [n_dof x 1]
%              R(I) - Total reaction acting on global DOF I
%--------------------------------------------------------------------------

F_ext_L = zeros(size(vL,1), 1);
F_ext_R = zeros(size(vR,1), 1);

for i = 1:size(vR,1)
    F_ext_R(i, 1) = Fext(vR(i,1),1);   
end

for i = 1:size(vL,1)
    F_ext_L(i, 1) = Fext(vL(i,1),1);
end

KLL = zeros(size(vL, 1));
KRR = zeros(size(vR, 1));
KRL = zeros(size(vR, 1), size(vL, 1));
KLR = zeros(size(vL, 1), size(vR, 1));

for i = 1:size(KLL, 1)
   for j = 1:size(KLL, 2) 
   KLL(i,j) = KG(vL(i),vL(j));
   end
end

for i = 1:size(KRR, 1)
   for j = 1:size(KRR, 2) 
   KRR(i,j) = KG(vR(i),vR(j));
   end
end

for i = 1:size(KRL, 1)
   for j = 1:size(KRL, 2) 
   KRL(i,j) = KG(vR(i),vL(j));
   end
end

for i = 1:size(KLR, 1)
   for j = 1:size(KLR, 2) 
   KLR(i,j) = KG(vL(i),vR(j));
   end
end

uL = KLL\(F_ext_L - KLR*uR);
R = KRR*uR + KRL*uL - F_ext_R;

z = size(vR, 1) + size(vL, 1);
u = zeros(z, 1);

for i = 1:size(vL, 1)
   u(vL(i), 1) = uL(i,1);
end

for i = 1:size(vR, 1)
   u(vR(i), 1) = uR(i,1);
end

