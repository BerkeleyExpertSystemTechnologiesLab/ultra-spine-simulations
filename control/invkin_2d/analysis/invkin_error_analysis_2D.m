% invkin_error_analysis_2D
% Copyright 2017-18 Andrew P. Sabelhaus and Berkeley Emergent Space Tensegrities Lab
% This script does the error analysis of the 2D Model-Predictive Controller for the spine.
% May 2018: this script adapted to make the same plots as the ACC2017
% paper but for the 2D MPC data.

% July 2018 - for the inverse kinematics from mpc.

function [ errors ] = invkin_error_analysis_2D( file_name, path_to_data_folder)
% Inputs:
%   file_name = name of the data file, needs to include '.mat'
%   path_to_data_folder = location of the data file to read in.
%   plots_flag = create plots (1) or do not make plots (0).
% Outputs:
%   errors = a struct containing all the types of errors that were calculated 

% Read in the data file.
data_path = strcat( path_to_data_folder, file_name );
data = load(data_path);

% A quick guide to Shirley's variable naming conventions...
% xi_traj = reference state trajectory
% xi_cl = closed loop results (after running MPC)
% opt_params = struct with all parameters, see the MPC script for a def'n.

% ...not used here, but could be meaningful?
% u_traj = reference inv kin inputs
% u_cl = closed loop calculated controls from the MPC algorithm

% The algorithm that's currently used calculates a trajectory out past the
% last timestep: e.g., there are timesteps + horizon number of datapoints.
% Thus, remove the last N points from the reference state traj.
% Rows are states, columns are timesteps.
N = 4; %horizon length

xi_traj = data.xi_traj(:, 1 : end-N);
xi_cl = data.xi_cl; % for convenience, so we don't keep having to do data.whatever

% Append these to a struct to return.
% BUT, rename according to the 3D naming convention, that way we can re-use
% the 3D analysis code.
errors.ref_traj = xi_traj;
errors.result_traj = xi_cl;
errors.opt_params = data.opt_params;

% calculate what's needed for the combined plots. specifically, that's:

% Relative errors for the states:
tracking_error = xi_cl - xi_traj;
errors.tracking_error = tracking_error;

% save some other misc things for easier use
errors.dt = data.opt_params.dt;
% the number of points is really num_pts + 1, since this traj includes the
% 0-th state.
errors.num_pts = data.opt_params.num_pts + 1; % e.g., 399 + 1 = 400

% the 'combined' script will do all the plotting for us.
end


