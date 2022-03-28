function [vL,vR,uR] = applyCond(n_i,n_dof,fixNod)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_i      Number of DOFs per node
%                  n_dof    Total number of DOFs
%   - fixNod  Prescribed displacements data [Npresc x 3]
%              fixNod(k,1) - Node at which the some DOF is prescribed
%              fixNod(k,2) - DOF (direction) at which the prescription is applied
%              fixNod(k,3) - Prescribed displacement magnitude in the corresponding DOF
%--------------------------------------------------------------------------
% It must provide as output:
%   - vL      Free degree of freedom vector
%   - vR      Prescribed degree of freedom vector
%   - uR      Prescribed displacement vector
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering to
% determine at which DOF in the global system each displacement is prescribed.

vR = zeros(size(fixNod,1), 1);
vL = zeros(n_dof-size(fixNod,1), 1);
for i = 1:size(fixNod,1)
    vR(i, 1)=fixNod(i,2);
end

g = 1;
for j = 1:n_dof  
a = 0;
        for i=1:size(vR, 1)
            if j == vR(i,1)
            a = 1;
            end
        end
        if a == 0
            vL(g,1) = j;
            g=g+1;
        end    
 end

uR = zeros(size(fixNod,1), 1);
%Fer loop per casos diferents de 0