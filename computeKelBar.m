function [Kel, l_e_vector] = computeKelBar(n_i,n_el,x,Tn,E, Iz)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - Kel   Elemental stiffness matrices [n_el_dof x n_el_dof x n_el]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------
Kel = zeros(2*n_i, 2*n_i, n_el);
l_e_vector = zeros(n_el,1);

for i=1:n_el
    x_1_e= x(1,Tn(i,1));
    x_2_e= x(1,Tn(i,2));
    l_e= x_2_e-x_1_e;
    mat_aux = [12 6*l_e -12 6*l_e;
        6*l_e 4*l_e^2 -6*l_e 2*l_e^2
        -12 -6*l_e 12 -6*l_e;
        6*l_e 2*l_e^2 -6*l_e 4*l_e^2
        ];
    Kel(:,:,i)=Iz*E/(l_e^3) * mat_aux;
    l_e_vector(i,1) = l_e;
end



