% get_yalmip_controller_XYZT.m
% Copyright 2015 Abishek Akella, Andrew P. Sabelhaus
% This function contains the objective and constraints that work for the ULTRA Spine MPC when the reference trajectory has only nonzeros in X,Y,Z,T.
% Note that this is for a 4-vertebra (link == 3) spine system

function [controller, constraints, objective, parameters_in, solutions_out] = get_yalmip_controller_XYZT(N, inputs, states, ...
    A_t, B_t, c_t, prev_in, reference)

% Inputs:
% N, horizon length
% inputs, sdpvar of inputs to the system
% states, sdpvar of the system states
% A_t, sdpvar of (FILL IN)
% B_t, sdpvar of (FILL IN)
% c_t, sdpvar of (FILL IN)
% prev_in, sdpvar of the input at the previous timestep
% reference, sdpvar of the reference trajectory

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

% Build up the constraints
input_lim = .09*ones(24, 1); % Limit on length of cable allowed

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

% Build up the objective
% First, minimize deviations along trajectory for the topmost tetrahedron
% This controller is working on the first *four* states here, xyzT, which are reference{}(1:4) and states{}(25:28).
% For this current state:
objective = 25*norm(states{1}(25:28) - reference{1}(1:4), 2); 
% For the remaining states in this horizon:
for k = 2:(N-1)
    objective = objective + (1/2)*(25^k)*norm(states{k}(25:28) - reference{k}(1:4), 2) + (1/24)*(3^k)*norm(inputs{k} - inputs{k-1}, inf) ...
        + (3^k)*(norm(states{k}(25:28) - states{k-1}(25:28))); % Also minimize change in state/input for smooth motion
end
% Add the objective for the last state in this horizon. There are no inputs for index k ==  N.
objective = objective + (1/2)*(25^N)*norm(states{N}(25:28) - reference{N}(1:4), 2) + (3^N)*norm(states{N}(25:28) - states{N-1}(25:28));

parameters_in = {prev_in, states{1}, [A_t{:}], [B_t{:}], c_t, [reference{:}]};
solutions_out = {[inputs{:}], [states{:}]};

% Build controller object for faster computation during iteration
controller = optimizer(constraints, objective, sdpsettings('solver', 'gurobi', 'verbose', 1), parameters_in, solutions_out);