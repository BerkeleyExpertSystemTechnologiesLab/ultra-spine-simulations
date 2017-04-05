% get_ref_traj_invkin_XZG.m
% Revised from Drew's code 2015 for 3-dimensional to 2-dimensional version,
% Added derivatives too
% This function returns a trajectory for all three vertebra of a 4-vertebra spine that
% bends around the Y+ axis, according to the inverse kinematics script, in either direction.
% It includes full position state information for all three rigid bodies. 
% No velocities though, those are zero-padded.

function [traj, num_points] = get_ref_traj_invkin_XZG(tetra_vertical_spacing, num_points, direction,dt)
% Inputs:
%   tetra_vertical_spacing = the distance between successive vertebrae. On 2016-04-18, was 0.1 meters.
%   num_points = the number of timesteps/waypoints in this trajectory. On 2016-04-18, was 30 or 300.
%       On 2016-04-29: This will now be split in 1/2 between movement and regulation at the end state.
%   direction = either 1 or -1, for clockwise (1) or counterclockwise (-1) rotation.
% Outputs:
%   traj = the output trajectory of the whole 3-vertebra system. Will have 36 states.
%   num_points = number of waypoints in the trajectory

% Hardcode the number of moving vertebrae here
% TO-DO: make this a parameter.
num_vertebrae = 1; 

% Check to be sure "direction" is only 1 or -1, no scaling allowed.
assert( (direction == 1) | (direction == -1), 'Direction can only be 1 (clockwise) or -1 (counterclockwise)');

% On 2016-04-29, have this script only accept even numbers, just to avoid errors when splitting
% the trajectory in half between movement and regulation at the end state.
% On 2016-05-02, requirement removed since add_regulation_to_traj is a separate function
% and no division is required.
% assert( mod(num_points,2) == 0, 'Error, num_points is odd. Only even-numbered trajectory lengths are allowed
% for now.');

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
% pi/8 is approximately the max angle for vertebra 4 (3rd moving vertebra)
% at timestep 20 (max) from inv-kin on 2016-04-23.
beta_0 = 0;
beta_f = -direction * pi/8; 
% On 2016-09-18: made beta larger for illustrating the trajectory in a figure for the ACC 2017 paper.
% beta_f = direction * pi/4;
%beta_f = direction * pi/16;
%beta_f = direction * pi/32;
% beta_f = direction * pi/256;

% Number of points to have in this trajectory. 
% Note that it's been estimated that timesteps should only put the top tetras about 0.0014 units distance
% away from each other (in sequential
% timesteps) for the optimization to work. (Feb. 2016)
% num_points = 300;

% The sweep angles for the other tetrahedra are as the following, considering vertebra 1 to be the first 
% moving vertebra:
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
dx_ref = zeros(num_points, num_vertebrae);
dz_ref = zeros(num_points, num_vertebrae);
% x_ref = zeros(num_points/2, num_vertebrae);
% z_ref = zeros(num_points/2, num_vertebrae);

for i=1:num_vertebrae
    % Use the equations defined above, for this varying radius curve.
    %x_ref(:,i) = -0.014+c1 .* beta(:,i) .* sin(beta(:,i)) + (tetra_vertical_spacing * i) .* sin(beta(:,i));
    % Need to have a very slight offset so that Mallory's inverse kinematics code
    % does not give -inf (which occurs when the system is symmetric.)
    % x_ref(:,i) = 1e-3 + c1 .* beta(:,i) .* sin(beta(:,i)) + (tetra_vertical_spacing * i) .* sin(beta(:,i));
    x_ref(:,i) = c1 .* beta(:,i) .* sin(beta(:,i)) + (tetra_vertical_spacing * i) .* sin(beta(:,i));
    z_ref(:,i) = c1 .* beta(:,i) .* cos(beta(:,i)) + (tetra_vertical_spacing * i) .* cos(beta(:,i));

end

for i=1:num_points-1
    % Use the equations defined above, for this varying radius curve.
    dx_ref(i,:) = (x_ref(i+1,:)-x_ref(i,:))/dt;
    dz_ref(i,:) = (z_ref(i+1,:)-z_ref(i,:))/dt;
end
% For the last columns of each trajectory:
dx_ref(num_points,:) = dx_ref(num_points-1,:);
dz_ref(num_points,:) = dz_ref(num_points-1,:);

% Also, define the rotation of each tetra around its own axis. 
% This is variable gamma, state 5 out of 12 for each vertebra.
% Call this rate of change c2.
% In the inverse kinematics script, each vertebra rotated around its axis at evenly-space rates.
% For example, the rate of rotation in the 3rd vertebra was -0.03 rad/timestep.
% But, we must convert that to be a function of beta, not of timestep (see Drew's notes.)
% Luckily, beta is a linear function of timestep, and since gamma is also a linear function of timestep,
% we get the following linear relationships between beta and gamma, for each of the 3 tetras:
% TO-DO: find some reasonable relationship between these numbers. They look a bit like a power law?

% Note, no need to adjust these by clockwise or counterclockwise, since beta is changed directly above.
% c2 = [1.06, 1.39, 1.54, 2.1, 2.5];
c2 = -1.06;

g_ref = zeros(num_points, num_vertebrae);
dg_ref = zeros(num_points, num_vertebrae);
% g_ref = zeros(num_points/2, num_vertebrae);
for i=1:num_vertebrae
    % We already have our betas, just convert to gammas.
    g_ref(:,i) = c2(i) .* beta(:,i);
end

for i=1:num_points-1
    % Use the equations defined above, for this varying radius curve.
    dg_ref(i,:) = (g_ref(i+1,:)-g_ref(i,:))/dt;
end
% For the last columns of each trajectory:
dg_ref(num_points,:) = dg_ref(num_points-1,:);

% Finally, place all of these points into a big array of all points in the trajectory.
% We need this to be a concatenation of row vectors: traj is 36 rows by num_points columns.

% For each of the three tetrahedra:
traj = zeros(num_vertebrae * 6, num_points);

% Now, for the first half of the trajectory, plug in x_ref, z_ref, and g_ref.
for i=1:num_vertebrae
    % Plug in the x, z, and g references
    % The x position will be at 1, 13, 25
    traj( 6*(i-1) + 1, :) = x_ref(:,i)';
%     traj( 12*(i-1) + 1, 1:(num_points/2)) = x_ref(:,i)';
    % z is at 3, 15, 27
    traj( 6*(i-1) + 2, :) = z_ref(:,i)';
%     traj( 12*(i-1) + 3, 1:(num_points/2)) = z_ref(:,i)';
    % g is at 5, 17, 29
    traj( 6*(i-1) + 3, :) = g_ref(:,i)';
%     traj( 12*(i-1) + 5, 1:(num_points/2)) = g_ref(:,i)';
    traj( 6*(i-1) + 4, :) = dx_ref(:,i)';
    traj( 6*(i-1) + 5, :) = dz_ref(:,i)';
    traj( 6*(i-1) + 6, :) = dg_ref(:,i)';
end

% Then, copy the last state through to the end of the trajectory.
 %last_state = traj(:,num_points/2);
% traj(:, (num_points/2)+1 : end) = repmat(last_state, 1, num_points/2);

% end function.
    
    
    
    