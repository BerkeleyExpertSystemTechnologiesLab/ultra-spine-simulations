function [ xi_kp1 ] = simulate_2d_spine_dynamics( xi, inputs, dt, num_steps, dyn_type )
%simulate_d2_spine_dynamics.m
%   Forward-simulate the dynamics of the 2D tensegrity spine model, 
%   the "inverted-Y" shape, using the basic Euler integration method.
%   This function uses the dynamics functions, solved symbolically,
%   which is called two_d_spine_xi_dot.
%   Inputs to this function:
%       xi, the system states, \in R^6, which are: [x, z, theta, dx, dz, dTheta]
%       inputs, the rest lengths of each of the four cables in the system
%           (note that this is NOT "change in rest length", but the entire rest length.
%            Thus, equilibrium is not at inputs == 0.)
%       dt, the amount of time to forward-simulate the system.
%       num_steps, the number of steps to use in the Euler integration.
%       dyn_type, a parameter that controls which type of dynamics simulation to run.
%           See the dynamics derivation script, or below, for more information.
%   Outputs of this function:
%       xi_kp1, the system states xi after the simulation.
%           (note that I use "kp1" to mean k-plus-1.)

% Recommended use of this function is to use the system in the form of
% xi_kp1 = g( xi_k, u_k)
% with num_steps = 1 and a small dt,
% ...even though these dynamics are CONTINUOUS-TIME for now,
% this discretizes them reasonably well enough.

% Start: check the dimensions on the inputs.
% assert( all(size(xi) == [6 1]), 'xi is not the proper size: must be 6x1 column vector.');
%assert( all(size(inputs) == [4 1]), 'inputs is not the proper size: must be 4x1 column vector.');
assert( num_steps > 0, 'num_steps must be positive.');

% The time interval for each step is:
dt_step = dt / num_steps;

% 2018-05-17: Adding noise, for comparison with 3D MPC results.
% The 3D results used the following:
%noise_mag_pos = 0.0005;
%noise_mag_vel = 0.0002;
% To keep things consistent, compare with the bottom vertebra, which for
% the 3D case moved the same distance. Thus, with the same total length,
% just divide by the number of points, which here is 80 /400 as in 3D vs
% 2D. Since 80/400 = 0.2,
% ALSO! Since we're comparing top vertebra in 3D versus bottom vertebra in
% 2D, scale the noise again by the ratio of size travelled on the top
% vertebra in 3D to the bottom vertebra. With a bit of math from the 2017
% ACC paper, the top vertebra of the 3D spine moves 0.0015 m each step in
% x, and the bottom vertebra moves 0.00025 m each each step, so that's
% another conversion of 0.00025 / 0.0015 = 0.16 repeating.
% Multiplying 0.2 * 0.16 = 0.033. Wow! That's small!
%noise_mag_pos = 0.0005 * 0.033;
%noise_mag_vel = 0.0002 * 0.033;
% For the slightly larger spine vertebra, here's how we'll scale.
% The larger robot still has the 0.2 conversion due to datapoints. The step
% size is now, however, only about 1/2 that of the top vertebra of the
% larger spine. I estimate this based on looking at the max position (was
% about -12 in X for the top vert 3D, -2 for the lowest vertebra in 2D, and
% -6 for the larger single vertebra.) This makes sense since the 12/2 = 1/6
% = 0.16, what I used previously.
% So, now we've got 0.2 * 0.5 = 0.1, not 0.033.
% noise_mag_pos = 0.0005 * 0.1;
% noise_mag_vel = 0.0002 * 0.1;

% 2018-07-27: for the IROS workshop paper, more noise, just for
% illustration.
noise_mag_pos = 0.0005 * 0.2;
noise_mag_vel = 0.0002 * 0.2;

% Turn noise on or off.
noise_flag = 1;

% Forward simulate for the given number of steps
for i = 1:num_steps
    % At this timestep, we need to calculate xi_dot.
    % There are multiple ways to do this, depending on what type of dynamics
    % to run. Use a case-switch statement, based on the dynamics type passed in.
    switch dyn_type
        case 1
            % For dynamics approach 1, which
            % calculates the tensions and acceleration separately,
            % and does not constrain the tensions:
            % (a) calculate the tension
            tensions = two_d_spine_tensions(xi,inputs);
            % (b) calculate the accelerations
            accel = two_d_spine_accel(xi, tensions);
            % (c) re-form the state vector derivative
            xi_dot = zeros(size(xi));
            xi_dot(1:3) = xi(4:6);
            xi_dot(4:6) = accel;
        case 2
            % Approach 2 is like approach 1, but includes
            % a rectification step for the cable tensions
            % (this makes it so cables cannot "push".)
            % (a) calculate the tension
            tensions = two_d_spine_tensions(xi,inputs);
            % (b) rectify the tensions: they must be nonnegative
            tensions = (tensions >= 0).*tensions;
            % (c) calculate the accelerations
            accel = two_d_spine_accel(xi, tensions);
            % (d) re-form the state vector derivative
            xi_dot = zeros(size(xi));
            xi_dot(1:3) = xi(4:6);
            xi_dot(4:6) = accel;
        case 3
            % Approach 3 has a combined xi_dot function which
            % calculates tensions and accelerations automatically.
            % Note that like approach 1, this does NOT account for
            % negative cable tensions, e.g., cables can "push" using 
            % this method.
            xi_dot = two_d_spine_xi_dot(xi, inputs);
        case 4
            % Approach 4 uses the xi_dot function with the logistic
            % barrier on the tensions included.
            xi_dot = two_d_spine_xi_dot_barrier(xi, inputs);
    end
    % Next, now that xi_dot is calculated, 
    % Forward-Euler integration:
    xi_kp1 = xi + dt_step*xi_dot;
    % 5) update xi for the next iteration
    xi = xi_kp1;
    % 6) Let's also check to see if we've violated any constraints
    %    on the system: specifically, did any of the cables "push"?
    %    All the cable tensions must be positive.
    tensions = two_d_spine_tensions(xi, inputs);
%     for j=1:length(tensions)
%         % If this tension was negative, print out a statement.
%         if( tensions(j) < 0 )
%             disp(strcat('ERROR: TENSION IS NEGATIVE. Cable number #', num2str(j), ...
%                 ', had tension :', num2str(tensions(j))));
%         end
%     end
end

% Added in May 2018: add noise for a meaningful comparison to the 3D model.
% TO-DO: the magnitude may need to be less here, since the distance changes
% are MUCH smaller between timesteps.
if( noise_flag )
   % Add to the positions
   xi_kp1(1:3) = xi_kp1(1:3) + noise_mag_pos*randn(1);
   % add to the velocities
   xi_kp1(4:6) = xi_kp1(4:6) + noise_mag_vel*randn(1);
end

% done! xi_kp1 is returned.

end







