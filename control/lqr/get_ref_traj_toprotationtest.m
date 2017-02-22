% get_ref_traj_toprotationtest.m
% Copyright 2015 Andrew P. Sabelhaus
% This function returns a trajectory for the topmost spine vertebra that rotates along one axis.

function [traj, num_points] = get_ref_traj_toprotationtest()
% No inputs.
% Outputs:
%   traj = the output trajectory of the topmost tetrahedron
%   num_points = number of waypoints in the trajectory

% The top tetra will rotate around one axis.
% Let's have it go 0 to pi/3.
start_deg = 0;
max_deg = pi/3;

% Number of points to have in this trajectory. 
% Note that it's been estimated that timesteps should only put the top tetras about 0.0014 units distance away from each other (in sequential
% timesteps) for the optimization to work.
% For pi/3 rotations only: take a guess. 
num_points = 5;

% Create a sequence of successive angles between min and max
theta = linspace(start_deg, max_deg, num_points);

% This trajectory bends in y and z.
traj = [zeros(1, num_points); ...   % X
        zeros(1, num_points); ...   % Y
        zeros(1, num_points); ...   % Z
        zeros(1, num_points); ...   % T
        theta; ...                  % G
        zeros(1, num_points); ...   % P
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points)];