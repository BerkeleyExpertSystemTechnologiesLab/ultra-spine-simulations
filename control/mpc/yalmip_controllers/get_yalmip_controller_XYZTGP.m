% get_yalmip_controller_XYZTGP.m
% Copyright 2015 Abishek Akella, Andrew P. Sabelhaus
% This function contains the objective and constraints that work for the ULTRA Spine MPC, with full weighting matrices
% passed in, instead of individual weights.
% Note that this is for a 4-vertebra (link == 3) spine system

function [controller, constraints, objective, parameters_in, solutions_out] = get_yalmip_controller_XYZTGP(N, inputs, states, ...
    A_t, B_t, c_t, prev_in, reference, Q_track, Q_smooth, r_smooth_mult, r_smooth_pow, opt_time_limit)

% Inputs:
% N, horizon length
% inputs, sdpvar of inputs to the system
% states, sdpvar of the system states
% A_t, sdpvar of (FILL IN)
% B_t, sdpvar of (FILL IN)
% c_t, sdpvar of (FILL IN)
% prev_in, sdpvar of the input at the previous timestep
% reference, sdpvar of the reference trajectory
% Q_track, the weighting matrix for reference tracking on the states
% Q_smooth, the weighting matrix for smoothness from one state to the next
% r_smooth_mult, scalar multiplicative weight for the infinity norm of successive control inputs. This is for control authority.
% r_smooth_pow, scalar power-law weight for the infinity norm of successive control inputs. This is for control authority.
% opt_time_limit, time limit for the optimization solver (in seconds).

% Outputs:
% controller, the YALMIP controller object for this optimization problem. Is created from all the other outputs, supplied here for convenience.
% constraints, objective, parameters_in, solutions_out: all the other YALMIP variables. These are only output from this script for debugging purposes,
% the controller variable is all that's needed.

% Recall: the states here are:
% 1: x
% 2: y
% 3: z
% 4: T
% 5: G
% 6: P
% 7: dx
% 8: dy
% 9: dz
% 10: dT
% 11: dG
% 12: dP
% repeated: x,y,z,T are then states 25:28 for the topmost tetrahedron

% Check if the passed-in trajectory is the correct length.
% Reference is a cell array of size 10, one for each state in the horizon. Each cell has 36 variables.
% So, check the size of the first cell.
assert(size(reference{1}, 1) == 36, 'Reference does not have 36 states.');

% The Q matrices passed in here are 36 x 36. Includes weights of both the longitudinal and angular coordinates.
% The R matrices pass in here are 24 x 24.

%% Build up the constraints
input_lim = .07*ones(24, 1); % Limit on length of cable allowed

constraints = [norm(inputs{1} - prev_in, inf) <= 0.02]; % Deviation from previous applied input to current input
for k = 1:(N-2)
    constraints = [constraints, states{k+1} == [A_t{:}]*states{k} + [B_t{:}]*inputs{k} + c_t, ...
        -input_lim <= inputs{k} <= input_lim, ...
        norm(inputs{k}(1:8) - inputs{1}(1:8), inf) <= 0.01, ... % Minimize deviation from first input (minimize linearization error)
        norm(inputs{k}(9:16) - inputs{1}(9:16), inf) <= 0.01, ...
        norm(inputs{k}(17:24) - inputs{1}(17:24), inf) <= 0.01];
end

constraints = [constraints, norm(inputs{N-1}(1:8) - inputs{1}(1:8), 2) <= 0.1, ...% Final input is given a wider tolerance
    norm(inputs{N-1}(9:16) - inputs{1}(9:16), 2) <= 0.1, ...
    norm(inputs{N-1}(17:24) - inputs{1}(17:24), 2) <= 0.1];

constraints = [constraints, states{N} == [A_t{:}]*states{N-1} + [B_t{:}]*inputs{N-1} + c_t, -input_lim <= inputs{N-1} <= input_lim];

for j = 1:(N-1)
    constraints = [constraints, norm(states{j}(1:6) - states{j+1}(1:6), inf) <= 0.02, ... % Minimize deviation from first state (minimize linearization error)
        norm(states{j}(13:18) - states{j+1}(13:18), inf) <= 0.03, ...
        norm(states{j}(25:30) - states{j+1}(25:30), inf) <= 0.04, ...
        states{j}(3) + .02 <= states{j}(15), ... % Maintain some distance between links to prevent collision
        states{j}(15) + .02 <= states{j}(27)];
end

constraints = [constraints, states{N}(3) + .02 <= states{N}(15), states{N}(15) + .02 <= states{N}(27)];

%% Build up the objective

% First, track the current state.
% NOTE that we do not do tracking on inputs here: instead, only do smoothness tracking on inputs.
% This is (x - r)' Q (x-r)
objective = (states{1} - reference{1})' * Q_track * (states{1} - reference{1});

% Then, track the states for the rest of the horizon.
% Weight successive states by 1/2 * weight ^2.
% Also weight the difference between successive states and successive inputs ("smoothness")
for k=2:(N-1)
    % State tracking
    objective = objective + (states{k} - reference{k})' * (1/2) * (Q_track.^k) * (states{k} - reference{k});
    % State smoothness
    objective = objective + (states{k} - states{k-1})' * (Q_smooth.^k) * (states{k} - states{k-1});
    % Control input smoothness. This is an infinity norm here, so only scalar weights are needed
    objective = objective + (r_smooth_mult)*(r_smooth_pow^k) * norm(inputs{k} - inputs{k-1}, inf);
end

% Add the objective for the last state in this horizon. There are no inputs for index k == N.
objective = objective + (states{N} - reference{N})' * (1/2) * (Q_track.^N) * (states{N} - reference{N});
objective = objective + (states{N} - states{N-1})' * (Q_smooth.^N) * (states{N} - states{N-1});

% For the controller, need to include input parameters variables and solutions output variables
parameters_in = {prev_in, states{1}, [A_t{:}], [B_t{:}], c_t, [reference{:}]};
solutions_out = {[inputs{:}], [states{:}]};

% Build controller object for faster computation during iteration
controller = optimizer(constraints, objective, sdpsettings('solver', 'gurobi', 'gurobi.TimeLimit', opt_time_limit, 'verbose', 1), parameters_in, solutions_out);



