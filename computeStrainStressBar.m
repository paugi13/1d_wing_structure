function [eps,sig] = computeStrainStressBar(n_d,n_el,u,Td,x,Tn,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
%   - u     Global displacement vector [n_dof x 1]
%            u(I) - Total displacement on global DOF I
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
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
%   - eps   Strain vector [n_el x 1]
%            eps(e) - Strain of bar e
%   - sig   Stress vector [n_el x 1]
%            sig(e) - Stress of bar e
%--------------------------------------------------------------------------
R_e = zeros(2*n_d, 2*n_d, n_el);
l_e = zeros(n_el,1);
E_e = zeros(n_el,1);
for i=1:n_el
    x_1_e= x(Tn(i,1),1);
    x_2_e= x(Tn(i,2),1);
    y_1_e= x(Tn(i,1),2);
    y_2_e= x(Tn(i,2),2);
    l_e(i)= sqrt((x_2_e-x_1_e)^2+(y_2_e-y_1_e)^2);
    s_e=(y_2_e-y_1_e)/l_e(i);
    c_e=(x_2_e-x_1_e)/l_e(i);
    mat_aux = [c_e s_e 0 0;
        -s_e c_e 0 0;
        0 0 c_e s_e;
        0 0 -s_e c_e
        ];
   R_e(:,:,i)= mat_aux;
   E_e(i) = Tmat(i, 1);
end

u_e_global = zeros(size(Td,2), 1, n_el);

%Assigning every node its displacement in element analysis. 
%The elements share nodes, so its displacement is applies equally to both
%elements.
%With Td we select the degree of liberty whose displacement must be applied
%to the analysis of the desired element. 

for i = 1:n_el      %Line selector in Td
    for j = 1:size(Td,2)    %Column selector in Td
        u_e_global(j, 1, i) = u(Td(i, j), 1);
    end
end

u_e_local = zeros(size(Td,2), 1, n_el);

for i = 1:n_el
    u_e_local(:, 1, i) = R_e(:,:,i)*u_e_global(:,1,i);
end

%eps must contanin every element's epsilon. 
%Same for sig. 

eps = zeros(n_el, 1);
sig = zeros(n_el, 1);

for i = 1:n_el
    eps(i,1) = 1/l_e(i)*[-1 0 1 0]*u_e_local(:, 1, i);
    sig(i,1) = E_e(i)*eps(i,1);
end


end