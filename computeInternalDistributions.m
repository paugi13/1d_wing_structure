function [Fe_y, Me_z, coeff_u, coeff_m] = computeInternalDistributions(l_e_vector,u, Nel, n_ne, n_i, Td, Kel)
% Compute internal distributions.  

Fe_y = zeros(Nel, 2);
Me_z = zeros(Nel, 2);
coeff_u = zeros(Nel, 4);
coeff_m = zeros(Nel, 3);

for e = 1:Nel
u_e = zeros(n_ne*ni,1);

    for i = 1:n_ne*n_i
        I = Td(e,i);
        u_e(i,1) = u(I);
    end
    
    F_e_int = Kel(:,:,e)*u_e;
    Fe_y(e,1) = -F_e_int(1,1);
    Fe_y(e,2) = F_e_int(3,1);
    Me_z(e,1) = -F_e_int(2,1);
    Me_z(e,2) = F_e_int(4,1);
    ord3_matrix =   [2 l_e_vector(e,1) -2 l_e_vector(e,1);
        -3*l_e_vector(e,1) -2*l_e_vector(e,1)^2 3*l_e_vector(e,1) -l_e_vector(e,1)^2;
        0 l_e_vector(e,1)^3 0 0;
        l_e_vector(e,1)^3];
    pol_coeff = (1/l_e_vector^3)*ord3_matrix*u_e;
    coeff_u(e, :) = [pol_coeff(1) pol_coeff(2) pol_coeff(3) pol_coeff(4)];
    coeff_m(e, :) = [3*pol_coeff(1) 2*pol_coeff(2) pol_coeff(3)];
end

