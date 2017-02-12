% get_ref_traj_eq.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This function returns a trajectory of system equillibrium values

function [xi_traj, u_traj, num_points] = get_ref_traj_eq(num_points, horizon)
% Inputs:
%   num_points = number of waypoints desired in the trajectory
%   num_states = number of states in system
% Outputs:
%   traj = the output trajectory of the topmost vertebra
%       The shape of traj is a num_states by num_points matrix where each
%       column corresponds to the reference state at that step
%   num_points = passing the variable through


% System equillibrium values (From Drew's dynamics simulation)
xi_eq = [-0.001; 0.1237; 0.0657; 0; 0; 0];
u_eq = [0.1; 0.12; 0.05; 0.06];

xi_traj = repmat(xi_eq,1,num_points+horizon);
u_traj = repmat(u_eq,1,num_points+horizon);