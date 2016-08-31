% analyze_ultra_spine_mpc_results.m
% Copyright 2016 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function analyzes the data from a single run of MPC,
% including graphs etc.

function [ mpc_error ] = analyze_ultra_spine_mpc_results( path, make_graphs )
% Inputs:
%   path = the file name of the .mat from a single run of MPC to read in and analyze.
%   make_graphs = a boolean flag. If set to true, graphs will be plotted.
% Outputs:
%   error_analysis = a matrix of ???

% Some notes on the data that's loaded in:
%
% traj is the trajectory to track, as loaded in by ultra_spine_mpc. It should 
%       be (12 * num_vertebrae) x (num_points_ref_traj), where num_points_ref_traj is also just the length of the simulation
%       TO-DO: is this true now even with regulation appended to the end of a trajectory??
%
% x_ref is actually the output of the MPC simulation. This was used by Abishek when he
%       used a series of LQR controllers after running MPC (this didn't make much sense,
%       so Drew removed it, but kept the naming.)
% But! It's easier to use refx, which is just x_ref converted into a matrix. See
% the code in ultra_spine_mpc_single_simulation.
% Also note that refx and x_ref only include the trajectories as they evolve, e.g, do not include the initial condition.
% So, refx(:,1) corresponds to traj(:,2), for example.

% first, load in the data.
load(path);

% The easiest way to do this calculation is to strip traj of its first column (which is the initial condition).
% That data isn't interesting (since the controller also started from there, the error would be zero by definition.)
traj = traj(:, 2:end);

% Then, subtract the reference from the state at each step, x_ref - traj (this is like x - reference in usual controls terminology.)
mpc_error = refx - traj;

% We'll use Q_track to determine which states were tracked during the MPC.
% Where Q_track is zero, that state was not tracked. Otherwise, it was part of the objective function of the MPC.
% When 


end

