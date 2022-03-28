function Td = connectDOFs(n_el,n_nod,n_i,Tn)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_el     Total number of elements
%                  n_nod    Number of nodes per element
%                  n_i      Number of DOFs per node
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering.
Td = zeros(n_el, n_i*n_nod);

for i=1:n_el
    Td(i,:) = [Tn(i,1)*2-1 Tn(i,1)*2 Tn(i,2)*2-1 Tn(i,2)*2];
end

% Td = [1 2 3 4;
%     3 4 5 6;
%     5 6 7 8;
%     9 10 11 12;
%     11 12 13 14;
%     13 14 15 16;
%     1 2 9 10;
%     1 2 11 12;
%     3 4 9 10;
%     3 4 11 12;
%     3 4 13 14;
%     5 6 11 12;
%     5 6 13 14;
%     5 6 15 16;
%     7 8 13 14;
%     7 8 15 16
% ];

end