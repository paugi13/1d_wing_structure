function F_el = computeElementForceVector(n_i, n_nod, Nel, x, l_e_vector, L1, L2, l, M)
% Cumpte every force vector from every element.
% Lift and density distribution is needed 

F_el = zeros(2*n_i, Nel);

% To pick the left side -> i = 1:(Nel-1)
for i = 2:n_nod     % Right side
    q_e = l*(0.8-0.2*cos(pi*x(i)/L1));
    lambda_e = 9.81*(M/(4*(L1+L2)) + 3*M*(L1-x(i))/(2*L2^2));
    if x(i)>5
        q_e = l*(1-((x(i)-L1)/L2))*(1+((x(i)-L1)/L2));
        lambda_e = 9.81*M/(4*(L1+L2));
    end
    q_e_res = q_e - lambda_e;
    F_el(:,i-1) = q_e_res*l_e_vector(i-1,1)/2*[1; l_e_vector(i-1,1)/6; 1;...
        -l_e_vector(i-1,1)/6];
end

