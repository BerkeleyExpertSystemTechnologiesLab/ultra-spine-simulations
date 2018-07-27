% get_ref_traj_line.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This function returns a straight line trajectory which holds the
% equilibrium z and theta, but moves horizontally

function [xi_ref, num_pts] = get_ref_traj_line(num_pts, horizon)
% Inputs:
%   num_points = number of waypoints desired in the trajectory
%   num_states = number of states in system
% Outputs:
%   xi_ref = the state reference trajectory of the topmost vertebra
%       The shape of xi_ref is a num_states by num_points matrix where each
%       column corresponds to the reference state at that step
%   num_points = passing the variable through

xi_eq = [-0.001; 0.1237; 0.0657; 0; 0; 0];

x_ref = linspace(xi_eq(1)-0.01,xi_eq(1)+0.01,num_pts+horizon);
z_ref = xi_eq(2)*ones(1,num_pts+horizon);
T_ref = xi_eq(3)*ones(1,num_pts+horizon);

xi_ref = [x_ref; z_ref; T_ref];