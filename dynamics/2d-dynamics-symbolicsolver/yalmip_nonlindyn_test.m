% yalmip_nonlindyn_test.m
% Andrew P. Sabelhaus
% This is a short script that tests if YALMIP will take
% the nonlinear dynamics of the 2D spine system as a constraint
% directly, instead of having to linearize to system to write
% dynamics constraints.

clear all;
close all;
clc;

% Load in the symbolic solutions to make into inline functions below
sym_vars = load('two_d_dynamics_symbolic_vars.mat');

%%

% Let's do a problem that minimizes the difference between
% the system  state and some desired state,
% for one control input.
% An equilibrium state, with the vertebra rotated slightly, is
xi_eq = [-0.001; 0.1237; 0.0657; 0; 0; 0];
% An initial state to start from is
xi_0 = [-0.01; 0.12; 0; 0; 0; 0];
% The corresponding u input is
%u_eq = [0.1; 0.12; 0.05; 0.06];
xi = sdpvar(6,2);
u = sdpvar(4,1);

tensions = sdpvar(4,1);

% The timestep of integration should be small
dt = 0.001;

% The objective function could be
obj = (xi(:,2)-xi_eq)'*ones(6)*(xi(:,2)-xi_eq);
% We then want to find a control input u
% that makes the next state of xi as close to equilibrium as possible.

% The dynamics are 
% xi_kp1 = xi_k + dt*xi_dot
% which is
% x_kp1 = xi_k + dt*g(xi_k, u_k)
% Where xi_dot is the function g:
%g = two_d_spine_xi_dot_barrier(xi(:,1), u);
%g = two_d_spine_xi_dot(xi(:,1), u);

% Try, instead, to separate out the tensions
% and accel functions (the dynamics split up into
% two parts.)
% tensions = two_d_spine_tensions(xi(:,1), u);
% accel = two_d_spine_accel(xi(:,1), tensions);
xi_dot = sdpvar(6,1);

tensions_fun = inline(sym_vars.two_d_spine_tensions);
accel_fun = inline(sym_vars.two_d_spine_accel);

tensions_soln = tensions_fun(u(1), u(2), u(3), u(4), xi(1), xi(2), xi(3), xi(4), xi(5), xi(6));
accel = accel_fun(tensions(1), tensions(2), tensions(3), tensions(4), xi(1), xi(2), xi(3), xi(6));

%%


% Constraints are 
% 1) assigning the initial state
% 2) nonlinear dynamics
% 3) nonnegative cable tensions
% constraints = [ xi(:,1) == xi_0, ...
%     xi(:,2) == xi(:,1) + dt*g, ...
%     u >= zeros(4,1)];
constraints = [ tensions_soln == tensions, ...
    xi(:,1) == xi_0, ...
    xi_dot(1:3) == xi(4:6,1), ...
    xi_dot(4:6) == accel, ...
    xi(:,2) == xi(:,1) + dt*xi_dot, ...
    u >= zeros(4,1)];


% Try to solve:
diagnostics = optimize(constraints, [])

% ...nope, YALMIP complains. Tried this a bunch of ways, didn't work.
% This motivates a switch to something in Python, with a different solver.

