% simulate_mpc_traj.m
% Copyright 2016 Abishek Akella, Andrew P. Sabelhaus, Berkeley Emergent Space Technologies Lab
% This function simulates the MPC to build a reference trajectory for the LQR controller.
% This variant takes in a trajectory for the full 36-state system, not just for the top vertebra.

function [x_ref, u_ref, M] = simulate_mpc_traj(controller, systemStates, restLengths, links, ...
                                    dt, x, y, z, T, G, P, dx, dy, dz, dT, dG, dP, traj, N, noise)
% Simulate the MPC controller and find reference trajectory to follow
% provided waypoints. 

% Inputs:
% controller = yalmip controller object
% systemStates = a (links by 12, on 2016-04-19, 3x12) matrix of the current state of the system
% restLengths = a 1 x 8 vector of the rest lengths of vertical and saddle cables, per tetrahedron.
% links = scalar, the number of vertebra in the system
% dt = time between points in the trajectory
% x... dP = each, a 3x1 row vector of the states of the vertebra. TO-DO: this appears to be the same as systemStates. Remove.
% traj = the trajectory to track with the MPC controller. In this script, will be 36 x size(traj, 2).
% N = horizon length. As of 2016-04-19, = 10.
% noise = a flag, 0 or 1, which determines if a disturbance is added to the result of the dynamics calculations
%         (called "noise" here, but it's really a disturbance).

% Outputs:
% x_ref = the trajectory that the system followed under the influence of the MPC controller
% u_ref = the inputs used to get x_ref
% M = ?

% This function calls the following functions related to the spine
% dynamics:
%   linearize_dynamics
%   simulate_dynamics

disp('Simulating MPC trajectory...')

% Store the total length of the trajectory so it does not have to be re-calculated:
num_points_ref_traj = size(traj, 2);
% make it a string for displaying as output
num_points_ref_traj_str = num2str(num_points_ref_traj);

% Initializing states for each link based on default resting location of
% each individual link. Initial input to linearize around is u = {0}.
x_initial = [];
for k = 1:links
    x_initial = [x_initial; x(k); y(k); z(k); T(k); G(k); P(k); dx(k); dy(k); dz(k); dT(k); dG(k); dP(k)];
end
u_initial = zeros(8*links, 1);

% Linearize model around initial state/input.
[A, B, c] = linearize_dynamics(x_initial, u_initial, restLengths, links, dt);

M = 1;

prev_in = u_initial;

% Begin running MPC controller on the desired trajectory. The controller
% generates a dynamically feasible input and state trajectory which can
% later be followed using a faster control strategy.

for index = 1:(size(traj, 2) - N + 1)
    disp( strcat('Simulating MPC trajectory, timestep no.: ', num2str(index), ' out of ', num_points_ref_traj_str) );
    
    % Unroll the system states into a form that's used by the dynamics simulation
    for k = 1:links
        systemStates(k, 1) = x(k); systemStates(k, 2) = y(k); systemStates(k, 3) = z(k);
        systemStates(k, 4) = T(k); systemStates(k, 5) = G(k); systemStates(k, 6) = P(k);
        systemStates(k, 7) = dx(k); systemStates(k, 8) = dy(k); systemStates(k, 9) = dz(k);
        systemStates(k, 10) = dT(k); systemStates(k, 11) = dG(k); systemStates(k, 12) = dP(k);
    end
    
    tic;
    % Call YALMIP to get the controller inputs. Use N-step horizon controller
    outputs = controller{N}{{prev_in, reshape(systemStates', 36, 1), A, B, c, traj(:, index:(index+N-1))}}; 

    % Re-arrange the control inputs
    control = outputs{1}(:, 1);
    prev_in = control;

    % Simulate system using the results of the YALMIP optimization.
    % In other words, "apply the input from MPC"
    systemStates = simulate_dynamics(systemStates, restLengths, reshape(control, 8, 3)', dt, links, noise);
    
    % Reshape the output and store it.
    x_ref{M} = reshape(systemStates', 36, 1);
    u_ref{M} = control;
    
    % Linearize the dynamics around the new point (after having applied the control action)
    [A, B, c] = linearize_dynamics(reshape(systemStates', 36, 1), control, restLengths, links, dt);
    toc

    % Combine the system states back into their individual variables. This is for storing things later.
    for k = 1:links
        x(k) = systemStates(k, 1); y(k) = systemStates(k, 2); z(k) = systemStates(k, 3);
        T(k) = systemStates(k, 4); G(k) = systemStates(k, 5); P(k) = systemStates(k, 6);
        dx(k) = systemStates(k, 7); dy(k) = systemStates(k, 8); dz(k) = systemStates(k, 9);
        dT(k) = systemStates(k, 10); dG(k) = systemStates(k, 11); dP(k) = systemStates(k, 12);
    end

    % Increment M. This might have something to do with the movie?
    M = M + 1;
end

% Simulate states in the last N values of the desired trajectory (using the
% shrinking horizon controller)
% This is because the YALMIP optimizer must change for the end of the horizon.
% It's only tracking 9, then 8, then 7, then... 1 states ahead in the horizon.
% This is almost exactly the same as the above loop, except it indexes into the 
% cell array of YALMIP "controller"s according to horizon point.
for index = M:(size(traj, 2)-1)
    disp(strcat('Simulating MPC trajectory, timestep no.:',num2str(index)));
    for k = 1:links
        systemStates(k, 1) = x(k); systemStates(k, 2) = y(k); systemStates(k, 3) = z(k);
        systemStates(k, 4) = T(k); systemStates(k, 5) = G(k); systemStates(k, 6) = P(k);
        systemStates(k, 7) = dx(k); systemStates(k, 8) = dy(k); systemStates(k, 9) = dz(k);
        systemStates(k, 10) = dT(k); systemStates(k, 11) = dG(k); systemStates(k, 12) = dP(k);
    end

    tic;
    % Calculate number of states in traj from current state onwards
    shrinking_horizon = size(traj, 2) - M + 1;
    ref_traj = zeros(size(traj,1), N);    % THIS WAS CHANGED FROM 12 TO size(traj,1) ON 2016-04-22
    % Zero-pad reference trajectory passed into the controller (zero states
    % are ignored)
    ref_traj(:, 1:shrinking_horizon) = traj(:, index:end);
    % Use shrinking_horizon-step horizon controller
    outputs = controller{shrinking_horizon}{{prev_in, reshape(systemStates', 36, 1), A, B, c, ref_traj}}; 

    control = outputs{1}(:, 1);
    prev_in = control;

    % Simulate system dynamics based on MPC control input
    systemStates = simulate_dynamics(systemStates, restLengths, reshape(control, 8, 3)', dt, links, noise);
    x_ref{M} = reshape(systemStates', 36, 1);
    u_ref{M} = control;
    [A, B, c] = linearize_dynamics(reshape(systemStates', 36, 1), control, restLengths, links, dt);
    toc

    for k = 1:links
        x(k) = systemStates(k, 1); y(k) = systemStates(k, 2); z(k) = systemStates(k, 3);
        T(k) = systemStates(k, 4); G(k) = systemStates(k, 5); P(k) = systemStates(k, 6);
        dx(k) = systemStates(k, 7); dy(k) = systemStates(k, 8); dz(k) = systemStates(k, 9);
        dT(k) = systemStates(k, 10); dG(k) = systemStates(k, 11); dP(k) = systemStates(k, 12);
    end

    M = M + 1;
end
M = M - 1;

% end function.
end