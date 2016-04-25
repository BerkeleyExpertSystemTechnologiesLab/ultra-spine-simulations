% mpc_error_analysis.m
% Copyright 2016 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function loads in a set of saved data from an individual MPC run, calculates
% the error in trajectory tracking, and ...

function [ errors ] = mpc_error_analysis( start_time_string, path_to_data_folder )
% Inputs:
%   start_time_string = the prefix of the data file to read in, for analysis.
%       note that the remainder of the string is hard-coded below. TO-DO: change this.
%   path_to_data_folder = location of the data file to read in.
%       NOTE that this script currently ONLY WORKS for three-vertebrae spines.
% Outputs:
%   errors = a struct containing all the types of errors that were calculated 

% Read in the data file.
data_path = strcat( path_to_data_folder, 'ultra-spine-mpc_data_', start_time_string );
data = load(data_path);

% A quick sanity check: only 3-vertebra trajectories allowed as of 2016-04-25.
assert( size(data.traj, 1) == 36, 'MPC error script only works for three-vertebra spines. Use a different data file.');

% Pick out the reference trajectory and the resulting trajectory.
% Each of these should be num_states x length_trajectory
ref_traj = data.traj;
result_traj = data.refx;

% ref_traj is off by one. Prune out the first point in the reference trajectory (this point, the intial state, is not controlled for.)
ref_traj = ref_traj(:, 2:end);

% Take the difference between these two trajectories, and square each point: this is the squared tracking error (all states)
tracking_error_squared = (result_traj - ref_traj).^2;
% what does just the value of the error look like? These units will be in meters, at least, for comparison.
tracking_error = result_traj - ref_traj;
% Save this into the output struct
errors.tracking_error_squared = tracking_error_squared;
errors.tracking_error = tracking_error;

% Take the total tracking error, as a function of timestep, for the trajectories that were actually tracked
% using MPC. (ie., remove the tracking for velocities.)
% tracking_XYZ will be for the positions of each tetrahedron. 3 tetras, 3 states, = 9 trajectories.
% tracking_angle will be for rotations
tracking_XYZ_squared = zeros(9, size(tracking_error_squared, 2));
tracking_angle_squared = zeros(9, size(tracking_error_squared, 2));
tracking_XYZ = zeros(9, size(tracking_error, 2));
tracking_angle = zeros(9, size(tracking_error, 2));

for i=1:3
    % Pick out the XYZ and angle states
    tracking_XYZ_squared( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error_squared( (i-1)*12 + 1: (i-1)*12 + 3, :);
    tracking_angle_squared( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error_squared( (i-1)*12 + 4: (i-1)*12 + 6, :);
    tracking_XYZ( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error( (i-1)*12 + 1: (i-1)*12 + 3, :);
    tracking_angle( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error( (i-1)*12 + 4: (i-1)*12 + 6, :);
end

% save these into the output struct
errors.tracking_XYZ_squared = tracking_XYZ_squared;
errors.tracking_angle_squared = tracking_angle_squared;
errors.tracking_XYZ = tracking_XYZ;
errors.tracking_angle = tracking_angle;

% Smash each of these down to a single time series of scalars
tracking_XYZ_squared_sum = sum(tracking_XYZ_squared);
tracking_angle_squared_sum = sum(tracking_angle_squared);
tracking_XYZ_sum = sum(tracking_XYZ);
tracking_angle_sum = sum(tracking_angle);

% save these into the output struct
errors.tracking_XYZ_squared_sum = tracking_XYZ_squared_sum;
errors.tracking_angle_squared_sum = tracking_angle_squared_sum;
errors.tracking_XYZ_sum = tracking_XYZ_sum;
errors.tracking_angle_sum = tracking_angle_sum;

% add these two to get the total error for all position states
tracking_squared_total = tracking_XYZ_squared_sum + tracking_angle_squared_sum;
errors.tracking_squared_total = tracking_squared_total;
tracking_total = tracking_XYZ_sum + tracking_angle_sum;
errors.tracking_abs_total = tracking_total;


end






