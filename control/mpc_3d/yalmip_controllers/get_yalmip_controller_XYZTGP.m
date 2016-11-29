% get_yalmip_controller_XYZTGP.m
% Copyright 2015 Abishek Akella, Andrew P. Sabelhaus
% This function contains the objective and constraints that work for the ULTRA Spine MPC, with full weighting matrices
% passed in, instead of individual weights.
% Note that this is for a 4-vertebra (link == 3) spine system

function [controller, constraints, objective, parameters_in, solutions_out] = get_yalmip_controller_XYZTGP(N, inputs, states, ...
    A_t, B_t, c_t, prev_in, reference, Q_track, Q_smooth, r_smooth_mult, r_smooth_pow, stab_const, opt_time_limit)

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
% stab_const, the distance that defines the region for the stability constraint. E.g, system error must be within stab_const. (at point N).
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
num_states = size(reference{1});

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

% Add the stability constraint.
% This is defined as: reference tracking error must be less than stab_const, for the final point in the horizon.
% Hard-code this for the moment to see if it works (2016-04-25):
% Result: Nope, didn't work at all, system stayed around equilibrium position.
%constraints = [constraints, norm(states{N}(25:27) - reference{N}(25:27), 2) <= stab_const];

%% Build up the objective

% Scalar Multiplication Version: 
% For some reason, YALMIP does not seem to do well with vector multiplications on arrays of sdpvars.
% Or maybe my math is just wrong again (Drew 2016-04-25).
% So, here's a version with all the objectives written out with 2-norms.
% Loop through each of the states, and if the diagonal on Q is nonzero, add it to the objective.
% NOTE THAT THIS ASSUMES THAT NONZEROS in Q_track ARE ALSO NONZEROS in Q_smooth
for s=1:num_states
    % For each state, check the diagonal on Q_track
    if ( Q_track(s,s)  ~= 0 )
        % Save the weight for this state
        w_track = Q_track(s,s);
        % Get the smoothing weight for this state
        w_smooth = Q_smooth(s,s);
        
        % Add this state for all horizon points:
        % tracking (for k=1:N)
        % smoothing (for k=2:N)

        % State s, horizon point 1: tracking only
        % Need to check if 'objective' exists and create it if not.
        % TO-DO: more efficient.
        if ~( exist('objective') )
            % create objective this step
            objective = w_track * norm(states{1}(s) - reference{1}(s), 2); % If using 2^k below, need 2 here.
        else
            % just add to objective.
            objective = objective + w_track * norm(states{1}(s) - reference{1}(s), 2);
        end
        
        % State s, horizon point k:
        for k= 2:(N)
            % tracking:
            % On 2016-06-03, changed power function to 2^k instead of w_track^k. Maybe this will help with computation?
            % Used to be: objective = objective + (1/2) * w_track^k * norm(states{k}(s) - reference{k}(s), 2);
            % On 2016-08-30, reverted. Tracking is now weighted much more heavily (since w_track >> 2).
            objective = objective + w_track^k * norm(states{k}(s) - reference{k}(s), 2);
            % smoothing:
            % On 2016-06-03, removed power function from w_smooth. Now, smoothness of all horizon lengths is weighted equally.
            % On 2016-08-30, reverted. The smoothness of future horizon points is now weighted more heavily.
            objective = objective + w_smooth^k * norm(states{k}(s) - states{k-1}(s), 2);
        end
        
        % Put an extra weight on the final horizon point.
        % TO-DO: parameterize this last term as a passed-in constant.
        % On 2016-08-30, removed this extra weight. This was a failed attempt at trying to force the controller to converge
        % instead of oscillating with noise from the system linearization.
        %objective = objective + 30 * w_track * 2^N * norm(states{N}(s) - states{N}(s), 2);
    end
end

% Also, add input penalties to the objective function, according to horizon point.
for k=2:(N-1)
    % penalize the largest of the inputs at this horizon point
    % TO-DO: maybe a 2-norm would work better here?
    objective = objective + r_smooth_mult * r_smooth_pow^k * norm(inputs{k} - inputs{k-1}, inf);
    %objective = objective + r_smooth_mult * r_smooth_pow^k * norm(inputs{k} - inputs{k-1}, 2);
end

%% Old backup versions of objective-build-up
% Vector Multiplication Version:
% On 2016-04-25, this does not seem to work.
% % First, track the current state.
% % NOTE that we do not do tracking on inputs here: instead, only do smoothness tracking on inputs.
% % This is (x - r)' Q (x-r), so a sqrt is needed to make it a 2-norm.
% objective = sqrt( (states{1} - reference{1})' * (Q_track.^1).^2 * (states{1} - reference{1}) ); % The 1 is included here as a reminder for k=1.
% 
% % Then, track the states for the rest of the horizon.
% % Weight successive states by 1/2 * weight ^2.
% % Also weight the difference between successive states and successive inputs ("smoothness")
% for k=2:(N-1)
%     % State tracking
%     objective = objective + sqrt((states{k} - reference{k})' * ((1/2) * Q_track.^k).^2 * (states{k} - reference{k}));
%     % State smoothness
%     objective = objective + sqrt((states{k} - states{k-1})' * (Q_smooth.^k).^2 * (states{k} - states{k-1}));
%     % Control input smoothness. This is an infinity norm here, so only scalar weights are needed
%     objective = objective + (r_smooth_mult)*(r_smooth_pow^k) * norm(inputs{k} - inputs{k-1}, inf);
% end
% 
% % Add the objective for the last state in this horizon. There are no inputs for index k == N.
% objective = objective + sqrt((states{N} - reference{N})' * ((1/2) * Q_track.^N).^2 * (states{N} - reference{N}));
% objective = objective + sqrt((states{N} - states{N-1})' * (Q_smooth.^N).^2 * (states{N} - states{N-1}));
% % result: system didn't even move.

% % Debugging 2016-04-25: copy in the older controller_XYZ code, modified to fit the 32-state reference:
% % For this current state:
% objective = 25*norm(states{1}(25:27) - reference{1}(25:27), 2);
% % For the remaining states in this horizon:
% for k = 2:(N-1)
%     objective = objective + (1/2)*(25^k)*norm(states{k}(25:27) - reference{k}(25:27), 2) + (1/24)*(3^k)*norm(inputs{k} - inputs{k-1}, inf) ...
%         + (3^k)*(norm(states{k}(25:27) - states{k-1}(25:27))); % Also minimize change in state/input for smooth motion
% end
% % Add the objective for the last state in this horizon. There are no inputs for index k ==  N.
% objective = objective + (1/2)*(25^N)*norm(states{N}(25:27) - reference{N}(25:27), 2) + (3^N)*norm(states{N}(25:27) - states{N-1}(25:27));
% % result: this one worked, for circular trajectory at least. Seems to stabilize in the XYZ directions, but everything else is unstable, as expected.
% 
% % Debugging 2016-06-04: copy in the older controller_XYZ code, modified to fit the 32-state reference, tracking angles too:
% % For this current state:
% objective = 25*norm(states{1}(25:30) - reference{1}(25:30), 2);
% % For the remaining states in this horizon:
% for k = 2:(N-1)
%     objective = objective + (1/2)*(25^k)*norm(states{k}(25:30) - reference{k}(25:30), 2) + (1/24)*(3^k)*norm(inputs{k} - inputs{k-1}, inf) ...
%         + (3^k)*(norm(states{k}(25:30) - states{k-1}(25:30))); % Also minimize change in state/input for smooth motion
% end
% % Add the objective for the last state in this horizon. There are no inputs for index k ==  N.
% objective = objective + (1/2)*(25^N)*norm(states{N}(25:30) - reference{N}(25:30), 2) + (3^N)*norm(states{N}(25:30) - states{N-1}(25:30));
% % result: this one worked, for circular trajectory at least. Still not necessarily stable.

% % Debugging 2016-06-04: copy in the older controller_XYZ code, modified to fit the 32-state reference, tracking angles too, all 3 vertebrae:
% % For this current state:
% objective = 25*norm(states{1}(1:6) - reference{1}(1:6), 2);
% objective = objective + 25*norm(states{1}(13:18) - reference{1}(13:18), 2);
% objective = objective + 25*norm(states{1}(25:30) - reference{1}(25:30), 2);
% % For the remaining states in this horizon:
% for k = 2:(N-1)
%     objective = objective + (1/2)*(25^k)*norm(states{k}(1:6) - reference{k}(1:6), 2) + ...
%         (1/2)*(25^k)*norm(states{k}(13:18) - reference{k}(13:18), 2) + ...
%         (1/2)*(25^k)*norm(states{k}(25:30) - reference{k}(25:30), 2) + ...
%         (1/24)*(3^k)*norm(inputs{k} - inputs{k-1}, inf) + ...
%         (3^k)*(norm(states{k}(1:6) - states{k-1}(1:6))) + ...
%         (3^k)*(norm(states{k}(13:18) - states{k-1}(13:18))) + ...
%         (3^k)*(norm(states{k}(25:30) - states{k-1}(25:30))); % Also minimize change in state/input for smooth motion
% end
% % Add the objective for the last state in this horizon. There are no inputs for index k ==  N.
% objective = objective + ...
%       (1/2)*(25^N)*norm(states{N}(1:6) - reference{N}(1:6), 2) + (3^N)*norm(states{N}(1:6) - states{N-1}(1:6)) + ...
%       (1/2)*(25^N)*norm(states{N}(13:18) - reference{N}(13:18), 2) + (3^N)*norm(states{N}(13:18) - states{N-1}(13:18)) + ...
%       (1/2)*(25^N)*norm(states{N}(25:30) - reference{N}(25:30), 2) + (3^N)*norm(states{N}(25:30) - states{N-1}(25:30));
% % result: ...


% % Debugging 2016-04-25: copy in the older controller_XYZ code, change it to matrix multiplications.
% % This controller is working on the first *three* states here, xyz, which are reference{}(1:3) and states{}(25:27).
% % The weighting matrix Q for debugging should be 36x36, and have '25' on its diagonal for states 25:27.
% % BUT since we need to take the square root to get a 2-norm, Q needs to be squared before use.
% % LEAVE THE INPUTS AS-IS FOR NOW.
% Q_debugging = zeros(36);
% Q_debugging(25:27, 25:27) = 25 * eye(3);
% objective = sqrt( (states{1} - reference{1})' * (Q_debugging.^1).^2 * (states{1} - reference{1}) ); % Include the raised-to-1 here just for a reminder
% % For the remaining states in this horizon:
% for k = 2:(N-1)
%     objective = objective + ...
%         (1/2) * sqrt( (states{k} - reference{k})' * (Q_debugging.^k).^2 * (states{k} - reference{k}) ) + ...
%         (1/24)*(3^k)*norm(inputs{k} - inputs{k-1}, inf) + ...
%         (3^k)*(norm(states{k}(25:27) - states{k-1}(25:27))); % Also minimize change in state/input for smooth motion
% end
% % Add the objective for the last state in this horizon. There are no inputs for index k ==  N.
% objective = objective + ...
%     (1/2) * sqrt( (states{N} - reference{N})' * (Q_debugging.^N).^2 * (states{N} - reference{N}) ) + ...
%     (3^N)*norm(states{N}(25:27) - states{N-1}(25:27));
% % result: tracking was very odd, didn't work.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used for the ACC 2017 Paper: 
% Currently does not track velocities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Copy in the older controller_XYZ code, modified to fit the 32-state reference, tracking angles too, all 3 vertebrae:
% % For this current state:
% objective = 25*norm(states{1}(1:6) - reference{1}(1:6), 2);
% %objective = objective + 25*norm(states{1}(13:18) - reference{1}(13:18), 2);
% %objective = objective + 25*norm(states{1}(25:30) - reference{1}(25:30), 2);
% % For the remaining states in this horizon:
% for k = 2:(N-1)
%     objective = objective + (1/2)*(25^k)*norm(states{k}(1:6) - reference{k}(1:6), 2) + ...
%         (1/2)*(25^k)*norm(states{k}(13:18) - reference{k}(13:18), 2) + ...
%         (1/2)*(25^k)*norm(states{k}(25:30) - reference{k}(25:30), 2) + ...
%         (1/24)*(3^k)*norm(inputs{k} - inputs{k-1}, inf) + ...
%         (3^k)*(norm(states{k}(1:6) - states{k-1}(1:6))) + ...
%         (3^k)*(norm(states{k}(13:18) - states{k-1}(13:18))) + ...
%         (3^k)*(norm(states{k}(25:30) - states{k-1}(25:30))); % Also minimize change in state/input for smooth motion
% end
% % Add the objective for the last state in this horizon. There are no inputs for index k ==  N.
% objective = objective + ...
%       (1/2)*(25^N)*norm(states{N}(1:6) - reference{N}(1:6), 2) + (3^N)*norm(states{N}(1:6) - states{N-1}(1:6)) + ...
%       (1/2)*(25^N)*norm(states{N}(13:18) - reference{N}(13:18), 2) + (3^N)*norm(states{N}(13:18) - states{N-1}(13:18)) + ...
%       (1/2)*(25^N)*norm(states{N}(25:30) - reference{N}(25:30), 2) + (3^N)*norm(states{N}(25:30) - states{N-1}(25:30));
% % result: ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Attach final parts of the objective, return.
% For the controller, need to include input parameters variables and solutions output variables
parameters_in = {prev_in, states{1}, [A_t{:}], [B_t{:}], c_t, [reference{:}]};
solutions_out = {[inputs{:}], [states{:}]};

% Build controller object for faster computation during iteration
controller = optimizer(constraints, objective, sdpsettings('solver', 'gurobi', 'gurobi.TimeLimit', opt_time_limit, 'verbose', 1), parameters_in, solutions_out);



