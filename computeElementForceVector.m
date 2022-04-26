function F_el = computeElementForceVector(n_ne, n_i, n_el, Tn, x, M, L1, L2, g, l)
% Compute every force vector from every element.
% Lift and density distribution is needed. 

 syms a
lambda_1 = M/(4*(L1+L2))+3*M/(2*L2^2)*(L1-a);
lambda_2 = M/(4*(L1+L2));
q_1 = l*(0.8-0.2*cos(pi*a/L1));
q_2 = l*(1-(a-L1)/L2)*(1+(a-L1)/L2);

% The force matrix has the same number of rows as DOFs has an element.
F_el=zeros(n_ne*n_i, n_el);

for i = 1:n_el

    x_e_1=x(1,Tn(i,1));
    x_e_2=x(1,Tn(i,2));
    l_e=abs(x_e_2-x_e_1);
% When both nodes are in the fist section of the wing.
    if x_e_1<L1 && x_e_2<=L1
        q_mean_e=(int(q_1,a,x_e_1,x_e_2)-int(lambda_1,a,x_e_1,x_e_2)*g)/l_e;
    end 
% When both nodes are in the second section of the wing.    
    if x_e_1>=L1 && x_e_2>L1
        q_mean_e=(int(q_2,a,x_e_1,x_e_2)-int(lambda_2,a,x_e_1,x_e_2)*g)/l_e;
    end 
% When there is one node in each section of the wing.    
    if x_e_1<L1 && x_e_2>L1
        q_mean_e=(int(q_1,a,x_e_1,L1)+int(q_2,a,L1,x_e_2)-int(lambda_1,a,x_e_1,L1)*g-int(lambda_2,a,L1,x_e_2)*g)/l_e;
    end
% The force element is stored.
    F_el(:,i)=q_mean_e*l_e/2*[1; l_e/6; 1; -l_e/6;];
end



