%% 2D Spine Inverse Kinematics
% The goal of this program is to solve for the force density values of
% cables in a 2D spine consisting of two stacked tetrahedra. The bottom
% tetra is fixed, and the top tetra is translated. If the configuration is
% feasible, then infinitely many solutions exist; the minimum-effort
% solution will be selected.
%
% Author: Mallory Daly
% Created: 12/1/16
% Modified: 12/1/16

clear; close all

%% Load Spine Parameters
load('spine_geometric_parameters_2D.mat')

% Geometric parameters
length = spine_geometric_parameters.l; % length of long bar
h = spine_geometric_parameters.h; % height from top to bottom of tetra
w = sqrt(length^2-(h/2)^2); % width from center of tetra

% Mass and force parameters
g = spine_geometric_parameters.g; % m/s^2
m = spine_geometric_parameters.m; % kg

%% Connectivity Matrices
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."
% NOTE: Schek's work finds equilbrium positions given load and position of
% fixed nodes. We want to find loads given positions of nodes.

% Rows 1-6 are bars
% Rows 7-10 are cables
% Columns 1-4 are free nodes
% Columns 5-8 are fixed nodes

%         1  2  3  4    
C_free = [1 -1  0  0;  %  1
          0  1 -1  0;  %  2
          0  0  0 -1;  %  3
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

%% Nodal Positions
% Coordinate system such that nodes 6 and 8 are on the x axis and nodes 5
% and 7 are centered on the z axis.

% Fixed nodal positions
%            5   6   7   8
x_fixed = [  0  -w   0   w]';
z_fixed = [h/2   0   h   0]';

% VARIABLES: Values to translate center
hTrans = 0; % horizontal translation
vTrans = h/2; % vertical translation

% Translate center
xCenter = x_fixed(1)+hTrans;
zCenter = z_fixed(1)+vTrans;

% Free nodal positions relative to base
%               1           2           3           4
x_free = [xCenter   xCenter-w     xCenter   xCenter+w]';
z_free = [zCenter zCenter-h/2 zCenter+h/2 zCenter-h/2]';

% Combined nodal positions
x = [x_free; x_fixed];
z = [z_free; z_fixed];

%% Lengths of Bars and Cables

% Rows 1-6 are bars
% Rows 7-10 are cables
l = [           length;  % 1
                   h/2;  % 2
                length;  % 3
                length;  % 4
                   h/2;  % 5
                length;  % 6
       norm(x(2)-x(6));  % 7
       norm(x(2)-x(7));  % 8
       norm(x(4)-x(7));  % 9
       norm(x(4)-x(8))]; %10
% The forces in the cables will depend on F = k*(l-l0)

% Diagonal length matrix
L = diag(l);

%% Set Up Optimization Problem
% See J. Friesen's "DuCTT: a Tensegrity Robot for Exploring Duct Systems."
% He adapts Schek's work to solve for force density in Section III.

% Define number of bars, cables, and nodes
r = 6; % bars
s = 4; % cables
n = 8; % nodes

% Define A matrix
A = [C'*diag(C*x);
     C'*diag(C*z)];
 
% Define p vector

% Take last s rows/columns A matrices, which correspond to the cables
As = A(:,end-(s-1):end);

% Moore-Penrose Pseudoinverse (1/2)
As_pinv = pinv(As);
V = eye(size(As_pinv*As))-As_pinv*As;

% Loading of the tetrahedrons (assume equal distribution of mass at point
% nodes)
px = zeros(n,1);
pz = -m*g*ones(n,1);
p = [px; pz];

% VARIABLE: Minimum minimum force densities for each cable
c = zeros(4,1);

% % Set up YALMIP
% yalmip('clear')
% w = sdpvar(size(V,2),1);

% % Define objective and constraint
% obj = w'*(V'*V)*w + 2*w'*V'*As_pinv*p;
% constr = As_pinv*p + V*w - c >= 0;

% Set up QUADPROG
H = 2*(V'*V);
f = (2*V'*As_pinv*p)';
A = -V;
b = -(c - As_pinv*p);

% Call QUADPROG
wOpt = quadprog(H,f,A,b)

% Find q
q = As_pinv*p + V*wOpt