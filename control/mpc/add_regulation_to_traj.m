% add_regulation_to_traj.m
% Copyright 2016 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function augments a reference trajectory with a 'regulation' trajectory,
% i.e., it copies the final state of a trajectory a specified number of times.

function [ full_traj, num_points_ref_traj ] = add_regulation_to_traj( ref_traj, num_points_ref_traj_regulation )
% Inputs:
%   ref_traj = the part of the trajectpory when the robot moves. This is, for example, from one of 
%       the get_ref_traj functions. Should be 
%   num_points_ref_traj_regulation = the number of times to copy the final state of ref_traj.
% Outputs:
%   full_traj = ref_traj plus the final state of ref_traj copied many times. Should be size (tracking + regulation).
%   num_points_ref_traj = the final size of full_traj, passed back for debugging purposes.

% If the trajectory should NOT be augmented, then num_points_ref_traj_regulation == 0, so just return ref_traj.
if( num_points_ref_traj_regulation == 0)
    full_traj = ref_traj;
else
    % Recall that traj is (num system states) by (num timesteps), e.g., 36 by something large.
    num_states = size(ref_traj, 1);
    num_points_ref_traj_tracking = size(ref_traj, 2);
    
    % Create a new matrix with the correct sizes
    full_traj = zeros(num_states, num_points_ref_traj_tracking + num_points_ref_traj_regulation);

    % Insert the tracking portion of the trajectory
    full_traj(:, 1:num_points_ref_traj_tracking) = ref_traj;

    % Copy the final state of the tracking portion up to the end of the array
    last_state = ref_traj(:,end);
    full_traj(:, num_points_ref_traj_tracking+1 : end) = repmat(last_state, 1, num_points_ref_traj_regulation);
    
    % end if-else
end

num_points_ref_traj = size(full_traj, 2);

% end function
end