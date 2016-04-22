% get_yalmip_controller_XZG.m
% Copyright 2015 Abishek Akella, Andrew P. Sabelhaus
% This function contains the objective and constraints that work for the ULTRA Spine MPC when the reference trajectory has only nonzeros in X,Z,G.
% Note that this is for a 4-vertebra (link == 3) spine system

function [controller, constraints, objective, parameters_in, solutions_out] = get_yalmip_controller_XZG(N, inputs, states, ...
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

% Check if the passed-in trajectory is the correct length.
% Reference is a cell array of size 10, one for each state in the horizon. Each cell has 36 variables.
% So, check the size of the first cell.
assert(size(reference{1}, 1) == 36, 'Reference does not have 36 states.');

%% Define weights
% Power-function weight for the objectives, used on the reference-tracking terms, for the longitudinal coordinates x,y,z
obj_w_ref_xyz = 10; %used to be 5
% Power-function weight for the objectives, used on the reference-tracking terms, for the angle G
obj_w_ref_angle = 3; % used to be 25
% Multiplicative weight for the objectives, used on the successive-states terms (smooth motion)
obj_w_smooth = 3;
% Power-function weight for the objectives, used on the successive-input terms (control authority, how-strong-is-the-motor)
obj_w_input_pow = 1.5;
% Multiplicative weight for the objectives, used on the successive-input terms (control authority, how-strong-is-the-motor)
obj_w_input_mult = 1/24;

%% Build up the constraints
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

%% Build up the objective

% First, minimize deviations along trajectory for the topmost tetrahedron
% This controller is working on the X, Z, and G states here, for all three tetrahedra, which are states 1, 3, 5 for each tetrahedra.
% That's states{}(1), (3), (13), (15), (25), (27) for X and Z, and states{}(5), (17), (29) for gamma.

% For this current state (the first one in the set that are given to the optimizer):
objective = obj_w_ref_xyz * norm(states{1}(1) - reference{1}(1), 2) ...         %X1
            + obj_w_ref_xyz * norm(states{1}(3) - reference{1}(3), 2) ...       %Z1
            + obj_w_ref_xyz * norm(states{1}(13) - reference{1}(13), 2) ...     %X2
            + obj_w_ref_xyz * norm(states{1}(15) - reference{1}(15), 2) ...     %Z2
            + obj_w_ref_xyz * norm(states{1}(25) - reference{1}(25), 2) ...     %X3
            + obj_w_ref_xyz * norm(states{1}(27) - reference{1}(27), 2) ...     %Z3
            + obj_w_ref_angle * norm(states{1}(5) - reference{1}(5), 2) ...     %G1
            + obj_w_ref_angle * norm(states{1}(17) - reference{1}(17), 2) ...   %G2
            + obj_w_ref_angle * norm(states{1}(29) - reference{1}(29), 2);      %G3

% For the remaining states in this horizon:
for k = 2:(N-1)
    objective = objective ... % Add terms onto the objective, successively.
                + (1/2)*(obj_w_ref_xyz^k) * norm(states{k}(1) - reference{k}(1), 2) ... % Tracking X1
                + (1/2)*(obj_w_ref_xyz^k) * norm(states{k}(3) - reference{k}(3), 2) ... % Tracking Z1
                + (1/2)*(obj_w_ref_xyz^k) * norm(states{k}(13) - reference{k}(13), 2) ... % Tracking X2
                + (1/2)*(obj_w_ref_xyz^k) * norm(states{k}(15) - reference{k}(15), 2) ... % Tracking Z2
                + (1/2)*(obj_w_ref_xyz^k) * norm(states{k}(25) - reference{k}(25), 2) ... % Tracking X3
                + (1/2)*(obj_w_ref_xyz^k) * norm(states{k}(27) - reference{k}(27), 2) ... % Tracking Z3
                + (1/2)*(obj_w_ref_angle^k) * norm(states{k}(5) - reference{k}(5), 2) ...   % Tracking G1
                + (1/2)*(obj_w_ref_angle^k) * norm(states{k}(17) - reference{k}(17), 2) ... % Tracking G2
                + (1/2)*(obj_w_ref_angle^k) * norm(states{k}(29) - reference{k}(29), 2) ... % Tracking G3
                + (obj_w_input_mult)*(obj_w_input_pow^k) * norm(inputs{k} - inputs{k-1}, inf) ...      % Deviation of inputs, control authority
                + (obj_w_smooth^k) * norm(states{k}(1) - states{k-1}(1), 2) ...      % Smooth motion, X1
                + (obj_w_smooth^k) * norm(states{k}(3) - states{k-1}(3), 2) ...      % Smooth motion, Z1
                + (obj_w_smooth^k) * norm(states{k}(13) - states{k-1}(13), 2) ...      % Smooth motion, X2
                + (obj_w_smooth^k) * norm(states{k}(15) - states{k-1}(15), 2) ...      % Smooth motion, Z2
                + (obj_w_smooth^k) * norm(states{k}(25) - states{k-1}(25), 2) ...      % Smooth motion, X3
                + (obj_w_smooth^k) * norm(states{k}(27) - states{k-1}(27), 2) ...      % Smooth motion, Z3
                + (obj_w_smooth^k) * norm(states{k}(5) - states{k-1}(5), 2) ...        % Smooth motion, G1
                + (obj_w_smooth^k) * norm(states{k}(17) - states{k-1}(17), 2) ...      % Smooth motion, G2
                + (obj_w_smooth^k) * norm(states{k}(29) - states{k-1}(29), 2);         % Smooth motion, G3
end

% Add the objective for the last state in this horizon. There are no inputs for index k == N.
        
objective = objective ... % Add terms onto the objective, successively.
            + (1/2)*(obj_w_ref_xyz^N) * norm(states{N}(1) - reference{N}(1), 2) ... % Tracking X1
            + (1/2)*(obj_w_ref_xyz^N) * norm(states{N}(3) - reference{N}(3), 2) ... % Tracking Z1
            + (1/2)*(obj_w_ref_xyz^N) * norm(states{N}(13) - reference{N}(13), 2) ... % Tracking X2
            + (1/2)*(obj_w_ref_xyz^N) * norm(states{N}(15) - reference{N}(15), 2) ... % Tracking Z2
            + (1/2)*(obj_w_ref_xyz^N) * norm(states{N}(25) - reference{N}(25), 2) ... % Tracking X3
            + (1/2)*(obj_w_ref_xyz^N) * norm(states{N}(27) - reference{N}(27), 2) ... % Tracking Z3
            + (1/2)*(obj_w_ref_angle^N) * norm(states{N}(5) - reference{N}(5), 2) ...   % Tracking G1
            + (1/2)*(obj_w_ref_angle^N) * norm(states{N}(17) - reference{N}(17), 2) ... % Tracking G2
            + (1/2)*(obj_w_ref_angle^N) * norm(states{N}(29) - reference{N}(29), 2) ... % Tracking G3
            + (obj_w_smooth^N) * norm(states{N}(1) - states{N-1}(1), 2) ...      % Smooth motion, X1
            + (obj_w_smooth^N) * norm(states{N}(3) - states{N-1}(3), 2) ...      % Smooth motion, Z1
            + (obj_w_smooth^N) * norm(states{N}(13) - states{N-1}(13), 2) ...      % Smooth motion, X2
            + (obj_w_smooth^N) * norm(states{N}(15) - states{N-1}(15), 2) ...      % Smooth motion, Z2
            + (obj_w_smooth^N) * norm(states{N}(25) - states{N-1}(25), 2) ...      % Smooth motion, X3
            + (obj_w_smooth^N) * norm(states{N}(27) - states{N-1}(27), 2) ...      % Smooth motion, Z3
            + (obj_w_smooth^N) * norm(states{N}(5) - states{N-1}(5), 2) ...        % Smooth motion, G1
            + (obj_w_smooth^N) * norm(states{N}(17) - states{N-1}(17), 2) ...      % Smooth motion, G2
            + (obj_w_smooth^N) * norm(states{N}(29) - states{N-1}(29), 2);         % Smooth motion, G3

% For the controller, need to include input parameters variables and solutions output variables
parameters_in = {prev_in, states{1}, [A_t{:}], [B_t{:}], c_t, [reference{:}]};
solutions_out = {[inputs{:}], [states{:}]};

% Build controller object for faster computation during iteration
controller = optimizer(constraints, objective, sdpsettings('solver', 'gurobi', 'verbose', 1), parameters_in, solutions_out);