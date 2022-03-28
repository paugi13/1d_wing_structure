%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date:
% Author/s:
%

clear;
close all;

%% INPUT DATA

% Material properties
E = 0;

% Cross-section parameters
t1 = 0;
t2 = 0;
h1 = 0;
h2 = 0;
b = 0;


% Other data
g = 9.81;
L1 = 5;
L2 = 10;
L = L1+L2;
Me = 0;
M = 0;

% Number of elements for each part
nel = [3,6];
Nel = 96; % Number of elements for the "exact" solution

%% PRECOMPUTATIONS

% Compute section: 
% A  - Section area 
% Iz - Section inertia

% Compute parameter l:
% l - Equilibrium parameter

% Plot analytical solution
fig = plotBeamsInitialize(L1+L2);

% Loop through each of the number of elements
for k = 1:length(nel)

    %% PREPROCESS
    
    % Nodal coordinates
    %  x(a,j) = coordinate of node a in the dimension j
    % Complete the coordinates
    
    % Nodal connectivities  
    %  Tnod(e,a) = global nodal number associated to node a of element e
    
    % Material properties matrix
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
    %  mat(m,3) = Section inertia of material m
    mat = [% Young M.        Section A.    Inertia 
              E,                A,         Iz;  % Material (1)
    ];

    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
        
    %% SOLVER
    
    % Compute:
    % u  - Displacements and rotations vector [ndof x 1]
    % pu - Polynomial coefficients for displacements for each element [nel x 4]
    % pt - Polynomial coefficients for rotations for each element [nel x 3]
    % Fy - Internal shear force at each elements's nodes [nel x nne]
    % Mz - Internal bending moment at each elements's nodes [nel x nne]
    
    %% POSTPROCESS
    
    % Number of subdivisions and plots
    nsub = Nel/nel(k);
    plotBeams1D(fig,x,Tnod,nsub,pu,pt,Fy,Mz)
    drawnow;
    
end

% Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel))),'location','northeast');