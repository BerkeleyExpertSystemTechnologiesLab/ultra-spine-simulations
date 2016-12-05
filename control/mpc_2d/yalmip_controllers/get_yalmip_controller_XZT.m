% get_yalmip_controller_XZT.m.m
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
% This function contains the formulation of the 2D spine MPC optimization
% problem (This is for a 2 vertebra system)

function [controller, constraints, objective, parameters_in, solution_out] = get_yalmip_controller_XZT ...
    (N, inputs, states, A_t, B_t, c_t, prev_in, reference)

% Inputs:
% inputs, sdpvar of inputs to the system
% states, sdpvar of the system states
% A_t, sdpvar of linearized state matrix at time t
% B_t, sdpvar of linearized input matrix at time t
% c_t, sdpvar of linearized output matrix at time t
% prev_in, sdpvar of the input at the previous timestep
% reference, sdpvar of the reference trajectory

% Outputs:
% controller, the YALMIP controller object for this optimization problem. Is created from all the other outputs, supplied here for convenience.
% constraints, objective, parameters_in, solutions_out: all the other YALMIP variables. These are only output from this script for debugging purposes,
% the controller variable is all that's needed.

% Recall, the states for this system are
% x
% z
% T (theta)
% dx
% dz
% dT

% Quick error check to make sure horizon is positive and at least 2
assert(N > 1, 'Horizon must be at least length 2')

%% Build constraints


%% Build objective


%% Build controller