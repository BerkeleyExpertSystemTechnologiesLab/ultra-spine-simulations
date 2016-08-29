% get_ref_traj_upright.m
% Copyright 2015 Andrew P. Sabelhaus
% This function returns a trajectory for all three vertebra of a 4-vertebra spine that
% stays fully upright in an equilibrium position.
% It includes full position state information for all three rigid bodies. No velocities though, those are zero-padded.

function [traj, num_points] = get_ref_traj_upright(tetra_vertical_spacing, num_points, direction)
% Inputs:
%   tetra_vertical_spacing = the distance between successive vertebrae. On 2016-04-18, was 0.1 meters.
%   num_points = the number of timesteps/waypoints in this trajectory. On 2016-04-18, was 30 or 300.
%       On 2016-04-29: This will now be split in 1/2 between movement and regulation at the end state.
%   direction = NOT USED IN THIS SCRIPT
% Outputs:
%   traj = the output trajectory of the whole 3-vertebra system. Will have 36 states.
%   num_points = number of waypoints in the trajectory

% Hardcode the number of moving vertebrae here
% TO-DO: make this a parameter.
num_vertebrae = 3; 

% Create the single state that defines upright equilibrium for the spine
traj = zeros(num_vertebrae * 12, num_points);

for i=1:num_vertebrae
    % z points are third in each rigid body's state vector
    % z is at 3, 15, 27
    traj( 12*(i-1) + 3, :) = tetra_vertical_spacing * i;
end

% end function.
    
    
    
    