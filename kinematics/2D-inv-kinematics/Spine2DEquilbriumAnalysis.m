%% Find Forces in 2D Spine Using Equilibrium Force Balance
% The goal of this program is to solve for the force density values of
% cables in a 2D spine consisting of two stacked tetrahedra. The bottom
% tetra is fixed, and the top tetra is translated. If the configuration is
% feasible, then infinitely many solutions exist; the minimum-effort
% solution will be selected.
%
% This script solves a force balance, in which the tetrahedron are used as
% solid bodies. We assume that the spine is sitting on a surface, so that
% there are vertical reaction forces at each of the two bottom nodes in
% contact.
%
% Authors: Mallory Daly and Ellande Tang
% Created: 12/8/16
% Modified: 12/9/16

% TO DO:
% - Make into function with inputs zi, min. cable tension, and spring
% constant, and outputs cable tensions and rest lengths
% - Change coordinate system to match Drew's (translate bottom tetra so
% middle is at (0,0)
% - Change tetra numbering system to match Drew's

clear; close all

%% Spine Parameters

% Number of bars, cables, and nodes
r = 6; % bars
s = 4; % cables
n = 8; % nodes

% Load spine data
load('spine_geometric_parameters_2D.mat')

% Geometric parameters
ll = spine_geometric_parameters.l; % m, length of long bars
h = spine_geometric_parameters.h; % m, height from top to bottom of tetra
ls = h/2; % m, length of short bar
w = sqrt(ll^2-(h/2)^2); % m, width from center of tetra

% Mass and force parameters
g = spine_geometric_parameters.g; % m/s^2, acceleration due to gravity
m = spine_geometric_parameters.m; % kg/node
M = m*4; % kg/tetra

%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."
% NOTE: Schek's work finds equilbrium positions given load and position of
% fixed nodes. We want to find loads given positions of nodes.

% Full connectivity matrix
% Rows 1-6 are bars
% Rows 7-10 are cables
% Columns 1-4 are "free" nodes (they are fixed later)
% Columns 5-8 are fixed nodes
%    1  2  3  4  5  6  7  8  
C = [1 -1  0  0  0  0  0  0;  %  1
     1  0 -1  0  0  0  0  0;  %  2
     1  0  0 -1  0  0  0  0;  %  3
     0  0  0  0  1 -1  0  0;  %  4
     0  0  0  0  1  0 -1  0;  %  5
     0  0  0  0  1  0  0 -1;  %  6
     0  0  0  1  0  0  0 -1;  %  7
     0  0  1  0  0  0 -1  0;  %  8
     0  1  0  0  0  0  0 -1;  %  9
     0  1  0  0  0  0 -1  0]; % 10

% Connection matrix of cables
Cs = C(end-(s-1):end,:);

%% VARIABLES: Translation and rotation of free tetrahedron
xTrans = 1/2*w; % horizontal translation
zTrans = 3/4*h; % vertical translation
theta = deg2rad(5);      % rotation (radians)

%% Nodal Positions
% Coordinate system such that nodes 3 and 4 are on the x axis and nodes 1
% and 2 are centered on the z axis.

% Nodal positions of bottom tetrahedra
%          1   2   3   4
x_bot = [  0   0   w  -w]';
z_bot = [  0 h/2 -h/2 -h/2]';

% Rotation matrix for given angle of theta
rot = [cos(theta), -sin(theta);
       sin(theta),  cos(theta)];
   
% Need to multiply each node by the rotation matrix and add the (x,z)
% offset from the xi state vector. Copy out (x,z) offset to a 2x4 matrix to
% make this easier.
xz_offset = repmat([xTrans; zTrans], 1, 4);

% Translate and rotate position of fixed nodes to get coordinates of free
% nodes
xz_top = rot*[x_bot'; z_bot'] + xz_offset;

% Extract x and z of free nodes
x_top = xz_top(1,:)';
z_top = xz_top(2,:)';

% Combined nodal positions
x = [x_bot; x_top];
z = [z_bot; z_top];

% Plot nodal positions
figure
plot(x_bot,z_bot,'k.','MarkerSize',10)
hold on
plot(x_top,z_top,'r.','MarkerSize',10)
% figure
% geometry.l = length;
% geometry.h = h;
% plot_2d_spine([xTrans zTrans theta],geometry)

%% Lengths of Bars and Cables

% Rows 1-6 are bars
% Rows 7-10 are cables
l = [ls;                             % 1
     ll;                             % 2
     ll;                             % 3
     ls;                             % 4
     ll;                             % 5
     ll;                             % 6
     norm([x(4),z(4)]-[x(8),z(8)]);  % 7
     norm([x(3),z(3)]-[x(7),z(7)]);  % 8
     norm([x(2),z(2)]-[x(8),z(8)]);  % 9
     norm([x(2),z(2)]-[x(7),z(7)])]; %10

% Diagonal length matrices
L = diag(l);
L_cables = diag(l(end-(s-1):end));

%% Solve for Reaction Forces
% Assume spine is sitting on a surface. Then there are vertical reaction
% forces at nodes 3 and 4.

% Solve AR*[R3; R4] = bR, where AR will always be invertible
AR = [1 1; 0 (x(3)-x(4))];
bR = [2*M*g; M*g*(x(1)-x(4))+M*g*(x(5)-x(4))];
R = AR\bR;
R3 = R(1);
R4 = R(2);

% Note R3 and R4 cannot be negative, but this constraint is not imposed
% here. However, a condition under which R3 or R4 becomes negative creates
% an infeasible problem for the cable tensions, so the issue solves itself.
% It's possible that a more generalized problem would need to solve for
% reaction forces using a solver in order to impose this constraint.

%% Equilibrium Force Equations

% Create vector of distance differences
dx = Cs*x;
dz = Cs*z;

% Define A*q = p, where q is a vector of the cable force densities and p is
% a vector of the external forces. Note that A is not a full rank matrix
% (not invertible).
A = [ -dx(1) -dx(2) -dx(3) -dx(4);  % horizontal forces, bottom tetra
       dx(1)  dx(2)  dx(3)  dx(4);  % horizontal forces, top tetra
      -dz(1) -dz(2) -dz(3) -dz(4);  % vertical forces, bottom tetra
       dz(1)  dz(2)  dz(3)  dz(4)]; % vertical forces, top tetra
p = [ 0; 0; M*g-R3-R4; M*g;];

%% Solve Problem for Minimized Cable Tension

% Minimum cable tension
minTension = 0; % N

% Solve with YALMIP
yalmip('clear')
q = sdpvar(s,1);
obj = q'*q;
constr = [A*q == p, L_cables*q >= minTension*ones(s,1)];
options = sdpsettings('solver','quadprog','verbose',0);
sol = optimize(constr,obj,options);
% optimize(constr,obj)

if sol.problem == 0
    % Display optimal force densities
    qOpt = value(q) % N/m

    % Display cable tensions
    tensionsOpt = L_cables*qOpt % N
else
    display(sol.info)
end


%% Check Distance Vectors Symbollically
% 
% x1 = sym('x1','real');
% x2 = sym('x2','real');
% x3 = sym('x3','real');
% x4 = sym('x4','real');
% x5 = sym('x5','real');
% x6 = sym('x6','real');
% x7 = sym('x7','real');
% x8 = sym('x8','real');
% x = [x1 x2 x3 x4 x5 x6 x7 x8]';
% 
% z1 = sym('z1','real');
% z2 = sym('z2','real');
% z3 = sym('z3','real');
% z4 = sym('z4','real');
% z5 = sym('z5','real');
% z6 = sym('z6','real');
% z7 = sym('z7','real');
% z8 = sym('z8','real');
% z = [z1 z2 z3 z4 z5 z6 z7 z8]';
% 
% Cs*x
% Cs*z
% 
% q7 = sym('q7','real');
% q8 = sym('q8','real');
% q9 = sym('q9','real');
% q10 = sym('q10','real');
% qs = [q7 q8 q9 q10]'