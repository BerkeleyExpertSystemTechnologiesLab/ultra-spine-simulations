% get_ref_traj_topbending2.m
% Copyright 2015 Andrew P. Sabelhaus
% This function returns a trajectory for the topmost spine vertebra that is along a constant-bending radius for the whole spine.
% Trajectory 2 is much shorter than trajectory 1.

function [traj, num_points] = get_ref_traj_topbending2()
% No inputs.
% Outputs:
%   traj = the output trajectory of the topmost tetrahedron
%   num_points = number of waypoints in the trajectory

% The whole spine will bend along a circle of radius (total spine length).
% Like the other trajectories (as of 2016-02-27), this assumes 4 spine vertebra.
r = 3 * 0.1; % this is tetra_vertical_spacing as of 2016-02-27

% how many degrees of movement should the top tetrahedron have?
% Note, this circle starts at theta = pi/2, and we are bending forwards here.
start_deg = pi/2;
max_deg = 4*pi / 9;

% Number of points to have in this trajectory. 
% Note that it's been estimated that timesteps should only be about 0.0014 units distance away from each other for the optimization to work.
% For 4*pi/9: (0.3 * pi/18) / 0.0014 = 38. Round it up to 50.
num_points = 50;

% Create a sequence of successive angles between min and max
theta = linspace(start_deg, max_deg, num_points);

% This trajectory bends in y and z.
traj = [zeros(1, num_points); ...
        -r*cos(theta); ...
        r*sin(theta); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ... 
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points)];