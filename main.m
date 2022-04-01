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

angle_aux = atan((h1/2-h2/2)/b);

% Iz inertia moments
Izz1 = 1/12*t2*h1^3;
Izz2 = 1/12*t2*h2^3;
Izz3 = 2*(1/12*t1*longitud^3*sin(angle_aux)^2 + A3/2*(h1/2-longitud/2*sin(angle_aux))^2);
Iz = Izz1 + Izz2 + Izz3;

%Iy inertia moments
Iyy1 = 1/12*h1*t2^3;
Iyy2 = 1/12*h2*t2^3;
Iyy3 = 1/12*t1*longitud^3*cos(angle_aux)^2;
Iyy = (Iyy1 + A1*z_cg^2) + (Iyy2 + A2*(z_cg-b)^2) + 2*(Iyy3 + A3/2*(z_cg-b/2)^2);

% Compute parameter l:
% l - Equilibrium parameter
aux_M = 1.7863e+04;
aux_Q = 32/3; %[N]
l = (aux_M*g)/aux_Q;

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


end

%%  Postprocessing

% Error calculus
%Taking into account the deflection at the wing tip
err = (u(193,1) - u(193,1))/u(193,1);
    
%Von Mises criterion
sig = y_max*Mz/Iz;

sigma = sqrt(sig^2 + 3*tau^2);
    
% for k = 1:length(nel)le
% end

% Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel))),'location','northeast');