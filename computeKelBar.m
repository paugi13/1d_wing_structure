function Kel = computeKelBar(n_d,n_el,x,Tn,mat,Tmat)
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
Kel = zeros(2*n_d, 2*n_d, n_el);

for i=1:n_el
    x_1_e= x(Tn(i,1),1);
    x_2_e= x(Tn(i,2),1);
    y_1_e= x(Tn(i,1),2);
    y_2_e= x(Tn(i,2),2);
    l_e= sqrt((x_2_e-x_1_e)^2+(y_2_e-y_1_e)^2);
    s_e=(y_2_e-y_1_e)/l_e;
    c_e=(x_2_e-x_1_e)/l_e;
    mat_aux = [c_e^2 c_e*s_e -(c_e)^2 -c_e*s_e;
        c_e*s_e (s_e)^2 -c_e*s_e -(s_e)^2;
        -(c_e)^2 -c_e*s_e (c_e)^2 c_e*s_e;
        -c_e*s_e -(s_e)^2 c_e*s_e (s_e)^2;
        ];
    Kel(:,:,i)=Tmat(i, 2)*Tmat(i, 1)/l_e * mat_aux;
end



