% linearize_dynamics_2d_new.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% Updated with disceretization methods by Shirley Zhao to work on
% ultra_spine_mpc_2d_alt_DT.m, Feb 2017
%
% This function calculates the linearized approximation to the dynamics of
% a 2D 2 vertebra spine. Linearization is performed numerically by
% calculating the Jacobians of the nonlinear dynamics: x_dot = f(x, u)
% 
% Disceretization is produced upon linearized model, and comes out with
% x_kp1 = A_k *x_k + B_k*u_k + c_k
%
% Inputs:
%   x_bar = nx by 1 vector of states, this is the state about which
%       the linearized A is calculated
%   u_bar = nu by 1 vector of imputs, this is the input about which
%       the linearized B is calculated
%   dt = timestep of the simulation
% Outputs:
%   A_k, B_k = Linearized and disceretized state matricies of the system
%       x(k+1) = A_k*x(k) + B_k*u(k)
%   c_k = 0th order term in the Taylor series expansion due to the
%   linearization about a trajectory point instead of an equillibrium point
%
%   An affine model provides different method to obtain c_k. However, under
%   current simulation for mpc, that one doesn't prove to work be
%
% This function calls the following functions related to the spine
% dynamics:
%   simulate_2d_spine_dynamics

function [A_k, B_k, c_k] = linearize_dynamics_2d_new(x_bar, u_bar, dt, dyn_type)

% Number of states and inputs
nx = length(x_bar);
nu = length(u_bar);

% Quick check: if any not-a-numbers are passed in, quit.
assert( ~(any(isnan(x_bar))), 'Error! NaN was passed in to linearize_dynamics as x_bar.');
assert( ~(any(isnan(u_bar))), 'Error! NaN was passed in to linearize_dynamics as u_bar.');

% Small constant for taking the numerical derivative (finite difference approx.)
% eps = 1e-5;
eps = 1e-10;

% Initialize linearized state space matricies
A = zeros(nx);
B = zeros(nx,nu);

% Linearize state matrix A
for i = 1:nx
    x_bar_U = x_bar;
    x_bar_L = x_bar;
    x_bar_U(i) = x_bar_U(i) + eps;
    x_bar_L(i) = x_bar_L(i) - eps;
    A(:,i) = (simulate_2d_spine_dynamics_new(x_bar_U,u_bar,dt,10,dyn_type) - ...
        simulate_2d_spine_dynamics(x_bar_L,u_bar,dt,1,dyn_type))/(2*eps);
end

% Linearize input matrix B
for i = 1:nu
    u_bar_U = u_bar;
    u_bar_L = u_bar;
    u_bar_U(i) = u_bar_U(i) + eps;
    u_bar_L(i) = u_bar_L(i) - eps;
    B(:,i) = (simulate_2d_spine_dynamics_new(x_bar,u_bar_U,dt,10,dyn_type) - ...
        simulate_2d_spine_dynamics(x_bar,u_bar_L,dt,1,dyn_type))/(2*eps);
end

% Discretization with Euler Method
% A_k = eye(nx) + eps * A;
% B_k = eps * B;
% c_k = simulate_2d_spine_dynamics(x_bar,u_bar,dt,10,dyn_type)-(A_k*x_bar+B_k*u_bar);

% Discretization with Zero-Order-Hold (ZOH) Method, using handwritten
% equations
% A_k = exp(A * eps);
% B_k = inv(A) * (A_k - eye(nx)) * B;
% %  A_k = eye(nx) + A_t * eps;
% %  B_k = B_t * eps;
% c_k = simulate_2d_spine_dynamics(x_bar,u_bar,dt,1,dyn_type)-(A_k*x_bar+B_k*u_bar);

% Discretization with Zero-Order-Hold (ZOH) Method, using embedded c2d command
% with A, B treated as a system     
     C = eye(nx);
     sysC = ss(A, B, C, 0);
     sysD = c2d(sysC, eps);
     A_k = sysD.A;
     B_k = sysD.B;
     c_k = simulate_2d_spine_dynamics(x_bar,u_bar,dt,1,dyn_type)-(A_k*x_bar+B_k*u_bar);