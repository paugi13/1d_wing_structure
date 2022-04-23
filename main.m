%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date: 
% Author/s:
%

clear;

%% INPUT DATA

% Material properties
E = 85e9;

% Cross-section parameters
t1 = 1.5e-3;
t2 = 4e-3;
h1 = 0.5;
h2 = 0.25;
b = 0.775;


% Other data
g = 9.81;
L1 = 5;
L2 = 10;
L = L1+L2;
Me = 2550;
M = 35000;

% Number of elements for each part
nel = [3,6, 12, 24, 48, 96];

%%  Essential parameters. 
% Everything is built so that it works by changing Nel's value.

Nel = 96; % Number of elements for the "exact" solution
n_d = 1;    %Problem dimension
n_ne = 2;   % Nodes in a beam
n_i = 2;    % DOF per node.

% Fixed degrees of freedom
% First node can't move upwards or rotate.
fixNod = [1 1 0;
    1 2 0;
];

%% PRECOMPUTATIONS

% Compute section: 
% A  - Section area 
h_aux = h1/2-h2/2;
longitud = sqrt(b^2+h_aux^2);
A1 = t2*h1;
A2 = t2*h2;
A3 = 2*longitud*t1;
A =  A1 + A2 + A3;

% Iz - Section inertia
z_cg = (A2*b + A3*b/2)/A;  

theta = atan((h1/2-h2/2)/b);

% Iz inertia moments -> z is the horizontal axis.
Izz1 = 1/12*t2*h1^3;
Izz2 = 1/12*t2*h2^3;
Izz3 = 2*(1/12*t1*longitud^3*sin(theta)^2 + A3/2*(h1/2-longitud/2*sin(theta))^2);
Iz = Izz1 + Izz2 + Izz3;

%Iy inertia moments -> Vertical axis.
Iyy1 = 1/12*h1*t2^3;
Iyy2 = 1/12*h2*t2^3;
Iyy3 = 1/12*t1*longitud^3*cos(theta)^2;
Iyy = (Iyy1 + A1*z_cg^2) + (Iyy2 + A2*(z_cg-b)^2) + 2*(Iyy3 + A3/2*(z_cg-b/2)^2);

% Compute parameter l:
% l - Equilibrium parameter
aux_M = (1.7863e+04);
aux_Q = 32/3; %[N]
l = (aux_M*g)/aux_Q;

u_ext = zeros(size(nel,2),1);

% This plot is empty at first
fig = plotBeamsInitialize(L1+L2);

for k = 1:length(nel)
    
 n_nod = nel(k)+1;  %Total number of nodes.
 n_dof = n_nod*n_i;  %Total number of degrees of freedom.
 
 %% PREPROCESS
 x = 0:L/nel(k):L;
    % Material properties matrix
    mat = [% Young M.        Section A.    Inertia 
              E,                A,         Iz;  % Material (1)
    ];
    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    
Tn = zeros(nel(k), n_ne);
for i = 1:nel(k)
    Tn(i,:) = [i, i+1];
end

%% SOLVER
    
    % Compute:
    % u  - Displacements and rotations vector [ndof x 1]
    % pu - Polynomial coefficients for displacements for each element [nel x 4]
    % pt - Polynomial coefficients for rotations for each element [nel x 3]
    % Fy - Internal shear force at each elements's nodes [nel x nne]
    % Mz - Internal bending moment at each elements's nodes [nel x nne]
    
% Compute Td: Only 2 degrees of freedom per node (deflection and rotation).
Td = connectDOFs(nel(k),n_ne,n_i,Tn);

% Independent element stiffness matrix
[Kel, l_e_vector] = computeKelBar(n_i,nel(k),x,Tn,E, Iz);

% Element force vector
F_el = computeElementForceVector(n_i, n_nod, nel(k), x, l_e_vector, L1, L2, l, M);

% Global matrix
[KG, Fext] = assemblyKG(nel(k),n_ne*n_i,n_dof,Td,Kel, F_el);

% Apply conditions 
[vL,vR,uR] = applyCond(n_i,n_dof,fixNod);

% System resolution
[u,R] = solveSys(vL,vR,uR,KG,Fext);

%Internal distributions
[Fy, Mz, pu, pt] = computeInternalDistributions(l_e_vector,u, nel(k), n_ne, n_i, Td, Kel);


%% POSTPROCESS

u_an = zeros(size(u,1)/2,1);
theta_an = zeros(size(u,1)/2,1);

j = 1;
for i = 1:2:size(u,1)
    u_an(j) = u(i,1);
    theta_an(j) = u(i+1,1);
    j = j+1;
end

Fy_an = zeros(nel(k)+1,1);
Mz_an = zeros(nel(k)+1,1);
Fy_an(2:nel(k)+1,1) = Fy(:,1);
Mz_an(2:nel(k)+1,1) = Mz(:,1);
Fy_an(1,1) = Fy(1,1);
Mz_an(1,1) = Mz(1,1);

% Number of subdivisions and plots
    nsub = Nel/nel(k);
    plotBeams1D(fig,x,Tn,nsub,pu,pt,Fy,Mz)
    drawnow;

% Vector acumulating displacement at x = L1 + L2
% 2*nel(k)+2 contains the rotation of the last node
u_ext(k, 1) = u(2*nel(k)+1,1);

% End of general loop
end

% Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel))),'location','northeast');


%%  Postprocessing
% Error calculus
error_vector = zeros(size(u_ext,1),1);

% Taking into account the deflection at the wing tip
% The last value is the one for Nel = 96 -> 'Exact solution'

for i = 1:size(u_ext,1)
   error_vector(i,1) = abs((u_ext(i,1) - u_ext(length(u_ext)))/u_ext(length(u_ext)));
end

figure
semilogx(nel, error_vector);
xlabel('log_{10}(n_{el})');
ylabel('Relative error');
title('Relative error along the number of elements');

%% VON MISES CRITERION
% Bending moment in z axis. 
% Axial stress -> h1/2 is the maximum value of y. 
sig_max = (h1/2)*Mz(1,1)/Iz;
Sy = Fy(1,1);
const = Sy/Iz;

% Shear stress
% Calculate Qs0
% Normal distance from x_cg to upper and lower bars.
u1 = 1/(sqrt(z_cg^2 + (h1/2)^2)) * [z_cg, h1/2];
u2 = [cos(theta) sin(theta)];
mu = acos(dot(u1,u2));
dn = sqrt(z_cg^2 + (h1/2)^2)*sin(mu);

% Qopen definite integrals
% Most solicitated position expected at 2-3. 
q01f = @(x) -const*t2*x;
q12f = @(x) -const*t1*(h2/2+x*sin(theta));
q23f = @(x) -const*t2*(h1/2-x);
q34f = @(x) -const*t1*(-h1/2+x*sin(theta));
q40f = @(x) -const*t2*(-h2/2+x);

q1op = integral(q01f, 0, h2/2);
q2op = q1op + integral(q12f, 0, longitud);
q3op = q2op + integral(q23f, 0, h1);
q4op = q3op + integral(q34f, 0, longitud);

% Qopen shear flow moment differentials in indefinite integrals (function_handle)
syms x
q01op_m_ind = matlabFunction(-(b-z_cg)*int(-const*t2*x));
q12op_m_ind = matlabFunction(-dn*(int(-const*t1*(h2/2+x*sin(theta)))+q1op));
q23op_m_ind = matlabFunction(-(-z_cg)*(int(-const*t2*(h1/2-x))+q2op));
q34op_m_ind = matlabFunction(-dn*(int(-const*t1*(-h1/2+x*sin(theta)))+q3op));
q40op_m_ind = matlabFunction(-(b-z_cg)*(int(-const*t2*(-h2/2+x))+q4op));

% q01op_m_ind = matlabFunction((b-z_cg)*int(const*t2*x));
% q12op_m_ind = matlabFunction(dn*(int(const*t1*(h2/2+x*sin(theta)))+q1op));
% q23op_m_ind = matlabFunction(-z_cg*(int(const*t2*(h1/2-x))+q2op));
% q34op_m_ind = matlabFunction(dn*(int(const*t1*(-h1/2+x*sin(theta)))+q3op));
% q40op_m_ind = matlabFunction((b-z_cg)*(int(const*t2*(-h2/2+x))+q4op));

% Shear flow (open) moment definite integrals
q01op_m = integral(q01op_m_ind, 0, h2/2);
q12op_m = integral(q12op_m_ind, 0, longitud);
q23op_m = integral(q23op_m_ind, 0, h1);
q34op_m = integral(q34op_m_ind, 0, longitud);
q40op_m = integral(q40op_m_ind, 0, h2/2);

% Final calculus.
M_total = q01op_m + q12op_m + q23op_m + q34op_m +q40op_m;
A_in = b*((h1-h2)/2) + b*h2;
Qs0 = M_total/(2*A_in);

q1op = q1op + Qs0;
q2op = q2op + Qs0;
q3op = q3op + Qs0;
q4op = q4op + Qs0;

tau_max = q2op/t1;

sigma_eq = sqrt(sig_max^2 + 3*tau_max^2);


% Sigma_eq (root) = 4.8281e+08



