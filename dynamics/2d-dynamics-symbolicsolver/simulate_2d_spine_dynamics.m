function [ xi_kp1 ] = simulate_2d_spine_dynamics( xi, inputs, dt, num_steps )
%simulate_d2_spine_dynamics.m
%   Forward-simulate the dynamics of the 2D tensegrity spine model, 
%   the "inverted-Y" shape, using the basic Euler integration method.
%   This function uses the dynamics functions, solved symbolically,
%   which are:
%       two_d_spine_getTensions
%       two_d_spine_accel
%   Note that the two_d_spine_getTensions function calls the following:
%       two_d_spine_lengths
%       two_d_spine_dlengths_dt
%   All these m-files must be in the same folder or in the MATLAB path.
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
% ...even though these dynamics are CONTINUOUS-TIME for now.
% To-Do: continuous -> discrete dynamics change. 
%   Is there a way to have Lagrange's equations do that automatically for us??

% Start: check the dimensions on the inputs.
assert( all(size(xi) == [6 1]), 'xi is not the proper size: must be 6x1 column vector.');
assert( all(size(inputs) == [4 1]), 'inputs is not the proper size: must be 4x1 column vector.');
assert( num_steps > 0, 'num_steps must be positive.');

% Pull out the individual states from xi
x = xi(1);
z = xi(2);
theta = xi(3);
dx = xi(4);
dz = xi(5);
dtheta = xi(6);

% The time interval for each step is:
dt_step = dt / num_steps;

% Forward simulate for the given number of steps
for i = 1:num_steps
    % At this timestep:
    % 1) Calculate the cable tensions
    %Tensions = two_d_spine_getTensions(x, z, theta, dx, dz, dtheta, ...
    %    inputs(1), inputs(2), inputs(3), inputs(4))
    % 2) Calculate the accelerations at this timestep
    %accel = two_d_spine_accel(x, z, theta, dx, dz, dtheta, ...
    %    Tensions(1), Tensions(2), Tensions(3), Tensions(4));
    % 3) Form the xi_dot vector out of these 3 accelerations
    %    as well as the velocities from the current state
    %xi_dot = [dx; dz; dtheta; accel(1); accel(2); accel(3)];
    % 4) Forward-integrate using Euler's method
    xi_dot = two_d_spine_xi_dot(xi, inputs);
    xi_kp1 = xi + dt_step*xi_dot;
    % 5) update xi for the next iteration
    xi = xi_kp1;
    % 5) un-roll the system states for the next iteration
    x = xi(1);
    z = xi(2);
    theta = xi(3);
    dx = xi(4);
    dz = xi(5);
    dtheta = xi(6);
end

% done! xi_kp1 is returned.

end







