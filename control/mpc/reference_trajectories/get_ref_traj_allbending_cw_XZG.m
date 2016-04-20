% get_ref_traj_allbending_cw_XZG.m
% Copyright 2015 Andrew P. Sabelhaus
% This function returns a trajectory for all three vertebra of a 4-vertebra spine that
% bends clockwise around the Y+ axis.
% It includes full position state information for all three rigid bodies. No velocities though, those are zero-padded.

function [traj, num_points] = get_ref_traj_allbending_cw_XZG(tetra_vertical_spacing)
% Inputs:
%   tetra_vertical_spacing = the distance between successive vertebrae. On 2016-04-18, was 0.1 meters.
% Outputs:
%   traj = the output trajectory of the whole 3-vertebra system. Will have 36 states.
%   num_points = number of waypoints in the trajectory

% Define this trajectory by determining the total angle of bending to perform. 
% From prior work, we found that a reasonable pattern of movement among the tetrahedra was 
% to have additive movement between successive spines: e.g., if vertebra 2 moved in Z by dZ and in G by dG,
% then vertebra 3 would move by 2*dZ and 2*dG, and vertebra 3 by 3*dZ and 3*dG.
% So, define two sets of angles: first, the total angle of translational movement of the top vertebra, with respect to the origin.
% Understanding: if a line was drawn between the center of the top vertebra and the origin, this would be the angle of change of that line.

% Call the translational movement XZ, since it's movement in that plane, and "3" for top tetrahedron.
start_deg = 0;
max_deg_XZ3 = pi/3;

% The second angle: the total rotation about its own axis for the top tetrahedron.
% To go clockwise, the system moves in +X, -Z, and +G.
max_deg_G3 = pi/12;

% Number of points to have in this trajectory. 
% Note that it's been estimated that timesteps should only put the top tetras about 0.0014 units distance away from each other (in sequential
% timesteps) for the optimization to work.
num_points = 30;

% Create a sequence of successive angles between min and max, for each of the vertebrae, for both angles.
% Here, theta_XZ1 refers to the angle for movement in the XZ plane for the first moving vertebra.
theta_XZ3 = linspace(start_deg, max_deg_XZ3, num_points);
theta_XZ2 = theta_XZ3 .* (2/3);
theta_XZ1 = theta_XZ3 .* (1/3);
theta_G3 = linspace(start_deg, max_deg_G3, num_points);
theta_G2 = theta_G3 .* (2/3);
theta_G1 = theta_G3 .* (1/3);

% Convert these into coordinates in x and z
x_3 = (3 * tetra_vertical_spacing) * sin(theta_XZ3);
z_3 = (3 * tetra_vertical_spacing) * cos(theta_XZ3);
x_2 = (2 * tetra_vertical_spacing) * sin(theta_XZ2);
z_2 = (2 * tetra_vertical_spacing) * cos(theta_XZ2);
x_1 = (1 * tetra_vertical_spacing) * sin(theta_XZ1);
z_1 = (1 * tetra_vertical_spacing) * cos(theta_XZ1);

% Finally, place all of these points into a big array of all points in the trajectory.
% We need this to be a concatenation of row vectors: traj is 36 rows by num_points columns.
    
% This trajectory rotates around the Z+ axis (e.g., in the X,Y plane).
traj = [ x_1; ...             % X
    zeros(1, num_points); ...   % Y
    z_1; ...                  % Z
    zeros(1, num_points); ...   % T
    theta_G1; ...               % G
    zeros(1, num_points); ...   % P
    zeros(1, num_points); ...   % dX
    zeros(1, num_points); ...   % dY
    zeros(1, num_points); ...   % dZ
    zeros(1, num_points); ...   % dT
    zeros(1, num_points); ...   % dG
    zeros(1, num_points); ...   % dP
    x_2; ...                  % X
    zeros(1, num_points); ...   % Y
    z_2; ...                  % Z
    zeros(1, num_points); ...   % T
    theta_G2; ...               % G
    zeros(1, num_points); ...   % P
    zeros(1, num_points); ...   % dX
    zeros(1, num_points); ...   % dY
    zeros(1, num_points); ...   % dZ
    zeros(1, num_points); ...   % dT
    zeros(1, num_points); ...   % dG
    zeros(1, num_points); ...   % dP
    x_3; ...                  % X
    zeros(1, num_points); ...   % Y
    z_3; ...                  % Z
    zeros(1, num_points); ...   % T
    theta_G3; ...               % G
    zeros(1, num_points); ...   % P
    zeros(1, num_points); ...   % dX
    zeros(1, num_points); ...   % dY
    zeros(1, num_points); ...   % dZ
    zeros(1, num_points); ...   % dT
    zeros(1, num_points); ...   % dG
    zeros(1, num_points)];      % dP
    

% end function.
    
    
    
    