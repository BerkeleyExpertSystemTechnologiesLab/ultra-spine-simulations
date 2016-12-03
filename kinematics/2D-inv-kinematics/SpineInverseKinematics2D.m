%% 2D Spine Inverse Kinematics
% The goal of this program is to solve for the force density values of
% cables in a 2D spine consisting of two stacked tetrahedra. The bottom
% tetra is fixed, and the top tetra is translated. If the configuration is
% feasible, then infinitely many solutions exist; the minimum-effort
% solution will be selected.
%
% Author: Mallory Daly
% Created: 12/1/16
% Modified: 12/2/16

clear; close all

%% Spine Parameters
% Specific to case: 2D spine with two tetrahedra, stacked vertically.

% Load spine data
load('spine_geometric_parameters_2D.mat')

% Geometric parameters
length = spine_geometric_parameters.l; % length of long bar
h = spine_geometric_parameters.h; % height from top to bottom of tetra
w = sqrt(length^2-(h/2)^2); % width from center of tetra

% Mass and force parameters
g = spine_geometric_parameters.g; % m/s^2
m = spine_geometric_parameters.m; % kg

% Number of bars, cables, and nodes
r = 6; % bars
s = 4; % cables
n = 8; % nodes

%% Connectivity Matrices
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."
% NOTE: Schek's work finds equilbrium positions given load and position of
% fixed nodes. We want to find loads given positions of nodes.

% Rows 1-6 are bars
% Rows 7-10 are cables
% Columns 1-4 are "free" nodes (they are fixed later)
% Columns 5-8 are fixed nodes

%         1  2  3  4    
C_free = [1 -1  0  0;  %  1
          0  1 -1  0;  %  2
          1  0  0 -1;  %  3
          0  0  0  0;  %  4
          0  0  0  0;  %  5
          0  0  0  0;  %  6
          0  1  0  0;  %  7
          0  1  0  0;  %  8
          0  0  0  1;  %  9
          0  0  0  1]; % 10
     
%          5  6  7  8    
C_fixed = [0  0  0  0;  %  1
           0  0  0  0;  %  2
           0  0  0  0;  %  3
           1 -1  0  0;  %  4
           1  0 -1  0;  %  5
           1  0  0 -1;  %  6
           0 -1  0  0;  %  7
           0  0 -1  0;  %  8
           0  0 -1  0;  %  9
           0  0  0 -1]; % 10
      
C = [C_free C_fixed];

%% VARIABLES: Translation and rotation of free tetrahedron
xTrans = 0*w; % horizontal translation
zTrans = 6/3*h; % vertical translation
theta = pi/20; % rotation (degrees)

%% Nodal Positions
% Coordinate system such that nodes 6 and 8 are on the x axis and nodes 5
% and 7 are centered on the z axis.

% Fixed nodal positions
%            5   6   7   8
x_fixed = [  0  -w   0   w]';
z_fixed = [h/2   0   h   0]';

% Rotation matrix for given angle of theta
rot = [cos(theta), -sin(theta);
       sin(theta),  cos(theta)];
   
% Need to multiply each node by the rotation matrix and add the (x,z)
% offset from the xi state vector. Copy out (x,z) offset to a 2x4 matrix to
% make this easier.
xz_offset = repmat([xTrans; zTrans], 1, 4);

% Translate and rotate position of fixed nodes to get coordinates of free
% nodes
xz_free = rot*[x_fixed'; z_fixed'] + xz_offset;

% Extract x and z of free nodes
x_free = xz_free(1,:)';
z_free = xz_free(2,:)';

% Combined nodal positions
x = [x_free; x_fixed];
z = [z_free; z_fixed];

% Plot nodal positions
% figure
% plot(x,z,'.','MarkerSize',10)
figure
plot(x_fixed,z_fixed,'k.','MarkerSize',10)
hold on
plot(x_free,z_free,'r.','MarkerSize',10)
% figure
% geometry.l = length;
% geometry.h = h;
% plot_2d_spine([xTrans zTrans theta],geometry)

%% Lengths of Bars and Cables
% Don't actually need this, but started making it because it was in Schek's
% paper.

% Rows 1-6 are bars
% Rows 7-10 are cables
l = [ length;  % 1
      h/2;  % 2
      length;  % 3
      length;  % 4
      h/2;  % 5
      length;  % 6
      norm([x(2),z(2)]-[x(6),z(6)]);  % 7
      norm([x(2),z(2)]-[x(7),z(7)]);  % 8
      norm([x(4),z(4)]-[x(7),z(7)]);  % 9
      norm([x(4),z(4)]-[x(8),z(8)])]; %10
% The forces in the cables will depend on F = k*(l-l0)

% Diagonal length matrix
L = diag(l);
L_cables = diag(l(7:10));

%% Loading of the tetrahedrons
% Assume equal distribution of mass at point nodes. Bottom four (fixed)
% nodes support weight of top four (free) nodes.
px = zeros(n,1);
pz = [-m*g*ones(n/2,1); m*g*ones(n/2,1)];
p = [px; pz];

%% Set Up Optimization Problem
% See J. Friesen's "DuCTT: a Tensegrity Robot for Exploring Duct Systems."
% He adapts Schek's work to solve for force density in Section III.

% Define A matrix
A = [C'*diag(C*x);
     C'*diag(C*z)];

% Take last s columns of the A matrix, which correspond to the cables
As = A(:,end-(s-1):end);

% VARIABLE: Minimum force densities for each cable
c = zeros(s,1);

% Moore-Penrose Pseudoinverse
As_pinv = pinv(As);
V = eye(size(As_pinv*As))-As_pinv*As;

% Set up QUADPROG
H = 2*(V'*V);
f = (2*V'*As_pinv*p)';
Aineq = -V;
bineq = -(c - As_pinv*p);

% Call QUADPROG
wOpt = quadprog(H,f,Aineq,bineq);

% Find q
qOpt = As_pinv*p + V*wOpt;

%% Find the Forces in the Cables
f = L_cables*qOpt

%% Set Up Optimization Problem (Bars and Cables)
% See J. Friesen's "DuCTT: a Tensegrity Robot for Exploring Duct Systems."
% He adapts Schek's work to solve for force density in Section III.
% 
% % Moore-Penrose Pseudoinverse
% A_pinv = pinv(A);
% V = eye(size(A_pinv*A))-A_pinv*A;
% 
% % VARIABLE: Minimum force densities for each cable
% c = zeros(r+s,1);
% 
% % Set up QUADPROG
% H = 2*(V'*V);
% f = (2*V'*A_pinv*p)';
% A = -V;
% b = -(c - A_pinv*p);
% 
% % Call QUADPROG
% wOpt = quadprog(H,f,A,b)
% 
% % Find q
% qOpt = A_pinv*p + V*wOpt
% 
% % Solve with YALMIP
% yalmip('clear')
% w = sdpvar(size(V,2),1);
% options = sdpsettings('solver','quadprog','verbose',2);
% obj = w'*(V'*V)*w + 2*w'*V'*As_pinv*p;
% constr = As_pinv*p + V*w - c >= 0;
% optimize(constr,obj,options)
% % optimize(constr,obj)
% wOpt_Y = value(w)
% qOpt_Y = As_pinv*p + V*wOpt_Y