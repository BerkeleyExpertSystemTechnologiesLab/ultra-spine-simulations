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
assert( all(size(xi) == [6 1]), 'xi is not the proper size: must be 6x1 column vector.');
%assert( all(size(inputs) == [4 1]), 'inputs is not the proper size: must be 4x1 column vector.');
assert( num_steps > 0, 'num_steps must be positive.');

% The time interval for each step is:
dt_step = dt / num_steps;

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
    for j=1:length(tensions)
        % If this tension was negative, print out a statement.
        if( tensions(j) < 0 )
            disp(strcat('ERROR: TENSION IS NEGATIVE. Cable number #', num2str(j), ...
                ', had tension :', num2str(tensions(j))));
        end
    end
end

% done! xi_kp1 is returned.

end







