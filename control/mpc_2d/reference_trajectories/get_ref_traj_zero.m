% get_ref_traj_zero.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This function returns a trajectory of zeroes - i.e., the controller's 
% goal will be to have the system stay still.

function [traj, num_points] = get_ref_traj_zero(num_points, horizon, num_states)
% Inputs:
%   num_points = number of waypoints desired in the trajectory
%   num_states = number of states in system
% Outputs:
%   traj = the output trajectory of the topmost vertebra
%       The shape of traj is a num_states by num_points matrix where each
%       column corresponds to the reference state at that step
%   num_points = passing the variable through

% This trajectory is only zeroes! The length accounts for 0 indexing
traj = zeros(num_states,num_points+horizon);