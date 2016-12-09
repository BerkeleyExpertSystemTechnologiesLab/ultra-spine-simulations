% linearize_dynamics_2d.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This function calculates the linearized approximation to the dynamics of
% a 2D 2 vertebra spine. Linearization is performed numerically by
% calculating the Jacobians of the nonlinear dynamics: x_dot = f(x, u)
% 
% Inputs:
%   x_bar = nx by 1 vector of states, this is the state about which
%       the linearized A is calculated
%   u_bar = nu by 1 vector of imputs, this is the input about which
%       the linearized B is calculated
%   dt = timestep of the simulation
% Outputs:
%   A, B = Linearized state matricies of the system
%       x(k+1) = A*x(k) + B*u(k)
%   c = 0th order term in the Taylor series expansion due to the
%   linearization about a trajectory point instead of an equillibrium point
%
% This function calls the following functions related to the spine
% dynamics:
%   simulate_2d_spine_dynamics

function [A, B, c] = linearize_dynamics_2d(x_bar, u_bar, dt)

% Number of states and inputs
nx = length(x_bar);
nu = length(u_bar);

% Small constant for taking the numerical derivative
eps = 1e-3;

% Initialize linearized state space matricies
A = zeros(nx);
B = zeros(nx,nu);

% Dynamics type flag
dyn_type = 2;

% Linearize state matrix A
for i = 1:nx
    x_bar_U = x_bar;
    x_bar_L = x_bar;
    x_bar_U(i) = x_bar_U(i) + eps;
    x_bar_L(i) = x_bar_L(i) - eps;
    A(:,i) = (simulate_2d_spine_dynamics(x_bar_U,u_bar,dt,1,dyn_type) - ...
        simulate_2d_spine_dynamics(x_bar_L,u_bar,dt,1,dyn_type))/(2*eps);
end

% Linearize input matrix B
for i = 1:nu
    u_bar_U = u_bar;
    u_bar_L = u_bar;
    u_bar_U(i) = u_bar_U(i) + eps;
    u_bar_L(i) = u_bar_L(i) - eps;
    B(:,i) = (simulate_2d_spine_dynamics(x_bar,u_bar_U,dt,1,dyn_type) - ...
        simulate_2d_spine_dynamics(x_bar,u_bar_L,dt,1,dyn_type))/(2*eps);
end

c = simulate_2d_spine_dynamics(x_bar,u_bar,dt,1,dyn_type)-(A*x_bar+B*u_bar);