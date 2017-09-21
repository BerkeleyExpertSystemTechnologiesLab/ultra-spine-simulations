function [ xi_kp1 ] = simulate_2d_spine_dynamics( xi, inputs, dt, num_steps )
%simulate_2d_spine_dynamics.m
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
assert( all(size(inputs) == [4 1]), 'inputs is not the proper size: must be 4x1 column vector.');
assert( num_steps > 0, 'num_steps must be positive.');

% The time interval for each step is:
dt_step = dt / num_steps;

% Forward simulate for the given number of steps
for i = 1:num_steps
    % At this timestep:
    % 1) Calculate the derivative of the system state, xi_dot,
    %    at this timestep
    xi_dot = two_d_spine_xi_dot(xi, inputs)
    % 4) Forward-integrate using Euler's method
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







