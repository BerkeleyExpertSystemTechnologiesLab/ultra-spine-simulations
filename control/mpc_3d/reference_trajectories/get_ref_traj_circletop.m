% get_ref_traj_circletop.m
% Copyright 2015 Andrew P. Sabelhaus, Abishek Akella
% This function returns a trajectory for the topmost spine vertebra that is an circle in (x, z), offset by some z amount.
% Trajectory is ONLY 12 STATES.

function [traj, num_points] = get_ref_traj_circletop(tetra_vertical_spacing, num_points, direction)
% Inputs:
%   tetra_vertical_spacing = the distance between successive vertebrae. On 2016-04-18, was 0.1 meters.
%   num_points = the number of timesteps/waypoints in this trajectory. On 2016-04-18, was 30 or 300.
%   direction = either 1 or -1, for clockwise (1) or counterclockwise (-1) rotation. NOT USED AS OF 2016-04-25
% Outputs:
%   traj = the output trajectory of the topmost tetrahedron
%   L = number of waypoints in the trajectory

% Number of vertebrae
num_vertebrae = 3;

% Radius of the circle
r = 0.04; % meters
% z-height offset for the circle
circle_height = tetra_vertical_spacing * num_vertebrae; % meters

% Multiplier here: used to be just 180
%num_points = 180;
theta = linspace(-pi, pi, num_points);

traj_top = [r*cos(theta) + r; ...
        zeros(1, num_points); ...
        r*sin(theta) + circle_height; ...
        zeros(1, num_points); ...
        zeros(1, num_points); ... 
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points); ...
        zeros(1, num_points)];

% full trajectory is zeros for bottom two tetras plus circle for top tetra.
% All tetras have 12 states (rigid bodies).
%traj = [ zeros(12, num_points); zeros(12, num_points); traj_top];

traj = traj_top;