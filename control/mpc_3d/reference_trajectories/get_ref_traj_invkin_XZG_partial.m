% get_ref_traj_invkin_XZG_partial.m
% Copyright 2015 Andrew P. Sabelhaus
% This function returns a trajectory for all three vertebra of a 4-vertebra spine that
% bends around the Y+ axis, according to the inverse kinematics script, in either direction.
% This script is a version of the trajectory with a smaller bending, only a fraction of the full one.
% It includes full position state information for all three rigid bodies. No velocities though, those are zero-padded.

function [traj, num_points] = get_ref_traj_invkin_XZG_partial(tetra_vertical_spacing, num_points, direction)
% Inputs:
%   tetra_vertical_spacing = the distance between successive vertebrae. On 2016-04-18, was 0.1 meters.
%   num_points = the number of timesteps/waypoints in this trajectory. On 2016-04-18, was 30 or 300.
%   direction = either 1 or -1, for clockwise (1) or counterclockwise (-1) rotation.
% Outputs:
%   traj = the output trajectory of the whole 3-vertebra system. Will have 36 states.
%   num_points = number of waypoints in the trajectory

% Hardcode the number of moving vertebrae here
% TO-DO: make this a parameter.
num_vertebrae = 3; 

% Check to be sure "direction" is only 1 or -1, no scaling allowed.
assert( (direction == 1) | (direction == -1), 'Direction can only be 1 (clockwise) or -1 (counterclockwise)');

% The inverse kinematics script currently translates the tetrahedra according to their angle
% from the origin. The radius of the arc drawn out by the tetra's center point varies with angle.
% Call that angle beta.
% For example, these curves are z = r cos(beta), where r = const * beta + r0.
% Here, r0 is the initial radius of this curve, which is the z height of each tetra.
% Define this constant for the curve. I called it c1 in my notes.
% From looking at the inverse kinematics results, it seemed like 1e-4 worked well enough.
% Anything much larger than this would define a curve that increases upwards in z (we want downward movement.)
c1 = 1e-4;

% Measure this trajectory by beta, from initial to final.
% Define it here for the top trajectory.
% Note that we include the "direction" flag here.
% pi/8 is approximately the max angle for vertebra 4 (3rd moving vertebra) at timestep 20 (max) from inv-kin on 2016-04-23.
% On 2016-06-04, let's do 'partial' as one-tenth of the full.
beta_0 = 0;
beta_f = direction * pi/8 * (1/10); 

% The sweep angles for the other tetrahedra are as the following, considering vertebra 1 to be the first moving vertebra:
% beta2 = 1.5 * beta1, beta3 = 2 * beta1.
% This equation looks like: vertebra_number_multiplier = 1 + 0.5 * (vertebra_number - 1).
% Equivalently, solving this for beta_i = (1/vertebra_number_multipler) * beta_3, going from vertebra i=1 to i=3, 
% gives (1/vertebra_number_multipler(i)) = 1 / (1 + (1/2) * (3 - i) ).
% So we can calculate a full set of beta vectors for the angles for all vertebrae.

beta = zeros(num_points, num_vertebrae);
% beta = zeros(num_points/2, num_vertebrae);
for i=1:num_vertebrae
    % For the i-th moving vertebra: create points from beta_0 to beta_f adjusted by the multiplier:
    % (remember that we're using this multipler here to "make the higher-up vertebrae move further")
    beta_f_current = beta_f * 1/( 1 + (1/2) * (3-i));
    beta(:,i) = linspace( beta_0, beta_f_current, num_points)';
%     beta(:,i) = linspace( beta_0, beta_f_current, num_points/2)';
end

% Then, create the longitudinal displacements for each tetrahedron.
% These are curves swept out with varying radius.
% Call these trajectories _ref.
x_ref = zeros(num_points, num_vertebrae);
z_ref = zeros(num_points, num_vertebrae);
% x_ref = zeros(num_points/2, num_vertebrae);
% z_ref = zeros(num_points/2, num_vertebrae);

for i=1:num_vertebrae
    % Use the equations defined above, for this varying radius curve.
    x_ref(:,i) = c1 .* beta(:,i) .* sin(beta(:,i)) + (tetra_vertical_spacing * i) .* sin(beta(:,i));
    z_ref(:,i) = c1 .* beta(:,i) .* cos(beta(:,i)) + (tetra_vertical_spacing * i) .* cos(beta(:,i));
end

% Also, define the rotation of each tetra around its own axis. This is variable gamma, state 5 out of 12 for each vertebra.
% Call this rate of change c2.
% In the inverse kinematics script, each vertebra rotated around its axis at evenly-space rates.
% For example, the rate of rotation in the 3rd vertebra was -0.03 rad/timestep.
% But, we must convert that to be a function of beta, not of timestep (see Drew's notes.)
% Luckily, beta is a linear function of timestep, and since gamma is also a linear function of timestep,
% we get the following linear relationships between beta and gamma, for each of the 3 tetras:
% TO-DO: find some reasonable relationship between these numbers. They look a bit like a power law?

% Note, no need to adjust these by clockwise or counterclockwise, since beta is changed directly above.
c2 = [1.06, 1.39, 1.54];

g_ref = zeros(num_points, num_vertebrae);
% g_ref = zeros(num_points/2, num_vertebrae);
for i=1:num_vertebrae
    % We already have our betas, just convert to gammas.
    g_ref(:,i) = c2(i) .* beta(:,i);
end

% Finally, place all of these points into a big array of all points in the trajectory.
% We need this to be a concatenation of row vectors: traj is 36 rows by num_points columns.

% For each of the three tetrahedra:
traj = zeros(num_vertebrae * 12, num_points);

% Now, for the first half of the trajectory, plug in x_ref, z_ref, and g_ref.
for i=1:num_vertebrae
    % Plug in the x, z, and g references
    % The x position will be at 1, 13, 25
    traj( 12*(i-1) + 1, :) = x_ref(:,i)';
    % z is at 3, 15, 27
    traj( 12*(i-1) + 3, :) = z_ref(:,i)';
    % g is at 5, 17, 29
    traj( 12*(i-1) + 5, :) = g_ref(:,i)';
end


% end function.
    
    
    
    