% linearize_dynamics_2d.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This function calculates the linearized AND DISCRETIZED 
% approximation to the dynamics of
% a 2D 2 vertebra spine. Steps of this script are:
%   1) Linearization is performed numerically by calculating the Jacobians 
%       of the nonlinear dynamics: x_dot = f(x, u)
%   2) Then, the affine linear dynamics are converted into a linear system
%       by augmenting the state space, x = [x; 1]
%   3) MATLAB's zero order hold discretization is used, c2d
%   4) The system is re-expanded into \dot x = A x + B u + C, where
%       each of the A, B, C are taken as blocks from the affine-augmented
%       system that pops out of c2d. This makes the MPC formulation work
%       appropriately: the affine-augmented discretized system seemed to 
%       have some trouble.
% 
% Inputs:
%   x_bar = nx by 1 vector of states, this is the state about which
%       the linearized A is calculated
%   u_bar = nu by 1 vector of imputs, this is the input about which
%       the linearized B is calculated
%   dt = timestep of the simulation
%
% Outputs:
%   A, B = Linearized and discretized state matricies of the system
%       x(k+1) = A*x(k) + B*u(k)
%   c = 0th order term in the Taylor series expansion due to the
%   linearization about a trajectory point instead of an equillibrium point
%
% This function calls the following functions related to the spine
% dynamics:
%   simulate_2d_spine_dynamics (used to be Shirley's "new" but that didn't
%   seem to ensure that the integration period was correct.)

function [A_k, B_k,c_k] = linearize_dynamics_2d_affine(x_bar,xp1_bar, u_bar, dt, dyn_type)

% Number of states and inputs
nx = length(x_bar);
nu = length(u_bar);

% Quick check: if any not-a-numbers are passed in, quit.
assert( ~(any(isnan(x_bar))), 'Error! NaN was passed in to linearize_dynamics as x_bar.');
assert( ~(any(isnan(u_bar))), 'Error! NaN was passed in to linearize_dynamics as u_bar.');

% Small constant for taking the numerical derivative (finite difference approx.)
% eps = 1e-5;
% Shirley used the following, summer 2017:
%eps = 1e-3;
% We should be using the same linearization/Jacobian step as the dt for
% the discretization.
eps = dt;

% The forward simulation of the nonlinear dynamics uses Euler integration,
% which gives a better result when is performed in smaller "chunks."
% So, let's declare a number of steps for the integration. Maybe 4? Seems
% reasonable.
num_integration_steps = 4;

% Initialize linearized state space matricies
A = zeros(nx);
B = zeros(nx,nu);

% 1) Linearize the system, calculate Jacobians / taylor series expansion.

% Linearize state matrix A
for i = 1:nx
    % Calculate the "upper" and "lower" limits for states and inputs,
    % since we'll be iterating over 
    x_bar_U = x_bar;
    x_bar_L = x_bar;
    x_bar_U(i) = x_bar(i) + eps*(xp1_bar(i) - x_bar(i));
    x_bar_L(i) = x_bar(i) - eps*(xp1_bar(i) - x_bar(i));
    % Now we've got a perturbed state and input vector, let's take the
    % finite difference between the newly perturbed states:
    % (state perturbed upper - state perturbed lower) / (2*timestep)
    A(:,i) = (simulate_2d_spine_dynamics(x_bar_U, u_bar, dt, ...
                    num_integration_steps, dyn_type) - ...
              simulate_2d_spine_dynamics(x_bar_L, u_bar, dt, ...
                    num_integration_steps, dyn_type)) / (2*eps);
    
    % Shirley used the following:
    %A(:,i) = (simulate_2d_spine_dynamics_new(x_bar_U,u_bar,dt,1,dyn_type) - ...
    %    simulate_2d_spine_dynamics_new(x_bar_L,u_bar,dt,1,dyn_type))/(2*eps);
end

% Linearize input matrix B
for i = 1:nu
    % Same setup as with A/states, but now we're perturbing the inputs.
    u_bar_U = u_bar;
    u_bar_L = u_bar;
    u_bar_U(i) = u_bar_U(i) + eps;
    u_bar_L(i) = u_bar_L(i) - eps;
    B(:,i) = (simulate_2d_spine_dynamics(x_bar, u_bar_U, dt, ...
                num_integration_steps, dyn_type) - ...
              simulate_2d_spine_dynamics(x_bar, u_bar_L, dt, ...
                num_integration_steps, dyn_type))/(2*eps);
    % Shirley used:
    %B(:,i) = (simulate_2d_spine_dynamics_new(x_bar,u_bar_U,dt,1,dyn_type) - ...
    %    simulate_2d_spine_dynamics_new(x_bar,u_bar_L,dt,1,dyn_type))/(2*eps);
end

% Finally, the constant offset between the linearized system and the 
% result of the nonlinear system. See taylor series expansion in paper.
c = simulate_2d_spine_dynamics(x_bar, u_bar, dt, num_integration_steps, ...
        dyn_type) - (A*x_bar+B*u_bar);
% Shirley used:
%c = simulate_2d_spine_dynamics_new(x_bar,u_bar,dt,1,dyn_type) - (A*x_bar+B*u_bar);

% 2) Convert linear system into affine system by concatenating the c
% vector.
A_af = [A, c;  ...
        zeros(1,nx), 1];
% Shirley used: (didn't we need a 1 as the bottom right corner element?)
%A_af = [A, c;zeros(1,nx+1)];

% The affine-augmented B matrix is not affected by the "1" state.
B_af = [B;zeros(1,nu)];

% 3) Covert the linear continuous time system into discrete time.
% The C matrix doesn't matter since we're doing state feedback
C = eye(nx+1);
% Continuous time system is: (no D, no input -> output feedthrough)
sysC = ss(A_af, B_af, C, 0);
% MATLAB's internal zero-order-hold function
sysD = c2d(sysC, eps);
% Pick out the discretized A and B matrices.
A_af_k = sysD.A;
B_af_k = sysD.B;

% 4) Finally, out of the discretized linear result, convert back into
% the affine system that it really is. (To-Do: check A_af_k and confirm
% that the nx+1-th row, column 0 to nx, is truly 0 still.)
A_k = A_af_k(1:nx,1:nx);
B_k = B_af_k(1,1:nu);
c_k = A_af_k(1:nx,nx+1);

end

% some older code...
% A_k = eye(nx) + eps * A;
% B_k = eps * B;
% c_k = simulate_2d_spine_dynamics(x_bar,u_bar,dt,10,dyn_type)-(A_k*x_bar+B_k*u_bar);

% A_k = exp(A * eps);
% B_k = inv(A) * (A_k - eye(nx)) * B;
% %  A_k = eye(nx) + A_t * eps;
% %  B_k = B_t * eps;
% c_k = simulate_2d_spine_dynamics(x_bar,u_bar,dt,1,dyn_type)-(A_k*x_bar+B_k*u_bar);
%      
%      C = eye(nx);
%      sysC = ss(A, B, C, 0);
%      sysD = c2d(sysC, eps);
%      A_k = sysD.A;
%      B_k = sysD.B;
%      c_k = simulate_2d_spine_dynamics(x_bar,u_bar,dt,1,dyn_type)-(A_k*x_bar+B_k*u_bar);
