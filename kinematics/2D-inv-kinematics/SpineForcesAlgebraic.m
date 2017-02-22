%% Find Forces in 2D Spine Algebraically
% This is a case where the second spine is translated vertically upwards.

% forceVec = [Rz F2 F10 F8 F9]

% Load spine data
load('spine_geometric_parameters_2D.mat')

% Geometric parameters
ll = spine_geometric_parameters.l; % length of long bars
h = spine_geometric_parameters.h; % height from top to bottom of tetra
ls = h/2; % length of short bar
w = sqrt(ll^2-(h/2)^2); % width from center of tetra

% Angle between bottom nodes of top tetra and top node of bottom tetra
theta = deg2rad(45);

% Mass and force parameters
g = spine_geometric_parameters.g; % m/s^2
m = spine_geometric_parameters.m; % kg/node
M = m*4; % kg/tetra

% Equilibrium equations (Ax = b)
A = [1 1 1 -sin(theta) -sin(theta);
     0 0 0 -cos(theta) cos(theta);
     0 -1 -1 sin(theta) sin(theta);
     0 -w w cos(theta)*ls -cos(theta)*ls;
     0 w -w (cos(theta)*h/2-sin(theta)*w) (sin(theta)*w-cos(theta)*h/2)];
b = [M*g; 0; M*g; 0; 0];

% Minimum cable tension
minTension = 1; % N

% Solve with YALMIP
yalmip('clear')
forceVec = sdpvar(5,1);
options = sdpsettings('solver','quadprog','verbose',2);
obj = forceVec'*forceVec;
constr = [A*forceVec == b, forceVec >= minTension];
optimize(constr,obj,options)
% optimize(constr,obj)
value(forceVec)