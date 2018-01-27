% ultra_spine_mpc_2d.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This is the primary script file that runs MPC on a 2d spine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script outline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script initialization

clear all;
close all;
clc;

disp('ULTRA Spine MPC 2D')

% Add various paths

% System dynamics path, spine geometric parameters file is also in this
% directory
dynamics_path = '../../dynamics/2d-dynamics-symbolicsolver';
addpath(dynamics_path)
% Inverse kinematics solver
inv_kin_path = '../../kinematics/2D-inv-kinematics';
addpath(inv_kin_path)
% Reference trajectories path
ref_traj_path = './reference_trajectories';
addpath(ref_traj_path)
% YALMIP controllers path
yalmip_controllers_path = './yalmip_controllers';
addpath(yalmip_controllers_path)

% Since the individual simulations will save their own videos and data, this script
% passes in the paths to the data and videos folders to the mpc function.
% Call this struct 'paths'.

% The path to our videos repository now needs to be set, since we're no longer pushing videos to this simulations repository.
% Drew (on 2016-04-19) put the ultra-spine-videos repository in the same folder as ultra-spine-simulations, which means
% the path to the MPC videos folder will be:
% paths.path_to_videos_folder = '../../../ultra-spine-videos/simulations/control/mpc/';

% The path to the folder where we'll store the .mat data from this simulation also needs to be set:
paths.path_to_data_folder = '../../data/mpc_2d_data/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define optimization and spine parameters

% Load spine gemetric parameters. This loads a struct called
% spine_geometric_parameters into the workspace. Contains the fields 'g',
% 'N', 'l', 'h', and 'm'
load('two_d_geometry.mat')

% Create a struct of optimization parameters
opt_params.num_pts = 39;
opt_params.num_states = 6;
opt_params.num_inputs = 4;
opt_params.horizon_length = 4;
opt_params.opt_time_lim = 1.5;
opt_params.spine_params = two_d_geometry;
opt_params.dt = 1e-4;

% Define initial states
% xi_0 = [-0.05; 0.15; pi/4; 0; 0; 0];
% u_0 = zeros(opt_params.num_inputs,1);
% opt_params.xi = [-0.01; 0.1; 0; 0; 0; 0];
% opt_params.xi = [-0.001+0.01; 0.1237+0.01; 0.0657+0.01; 0; 0; 0];
opt_params.u = zeros(opt_params.num_inputs,1);
% opt_params.u = [0.12;0.12;0.06;0.06];


% Define output matrix
C = eye(opt_params.num_states);
opt_params.num_outputs = size(C,1);

prev_u = opt_params.u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get trajectory for top vertebra

% Generate state reference trajectory. Functions for doing this should be 
% placed in the reference_trajectories folder and called here.
% [traj, ~] = get_ref_traj_zero(opt_params.num_pts,opt_params.horizon_length,opt_params.num_states);
% [xi_traj, u_traj, ~] = get_ref_traj_eq(opt_params.num_pts,opt_params.horizon_length);
[xi_traj, ~] = get_ref_traj_invkin_XZG(0.1,opt_params.num_pts+opt_params.horizon_length+1,1,opt_params.dt);
u_traj = zeros(opt_params.num_inputs,opt_params.num_pts+opt_params.horizon_length+1);
for i = 1:opt_params.num_pts+opt_params.horizon_length+1
    [~, u_traj(:,i)] = getTensions(xi_traj(:,i),opt_params.spine_params,30);
%     disp(xi_traj(:,i))
end
opt_params.xi = xi_traj(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create controller

% xi = sdpvar(opt_params.num_states,opt_params.horizon_length+1,'full');
% u = sdpvar(opt_params.num_inputs,opt_params.horizon_length,'full');
% xi_ref = sdpvar(opt_params.num_states,opt_params.horizon_length+1,'full');
% u_ref = sdpvar(opt_params.num_inputs,opt_params.horizon_length+1,'full');
% A_k = sdpvar(opt_params.num_states,opt_params.num_states);
% B_k = sdpvar(opt_params.num_states,opt_params.num_inputs);
% c_k = sdpvar(opt_params.num_states,1);
% prev_u = sdpvar(opt_params.num_inputs,1);
% % xi_0 = sdpvar(opt_params.num_states,1);
% 
% % Create N-step controller
% [controller,~,~,~,~] = get_yalmip_controller_XZT(opt_params.horizon_length, ...
%     xi,u,xi_ref,u_ref,A_k,B_k,c_k,prev_u,opt_params.spine_params, ...
%     opt_params.opt_time_lim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MPC Setup

% Define matrices for saving all states and inputs
xi_step = zeros(opt_params.num_states,opt_params.horizon_length+1,opt_params.num_pts);
u_step = zeros(opt_params.num_inputs,opt_params.horizon_length,opt_params.num_pts);
xi_cl = zeros(opt_params.num_states,opt_params.num_pts+1);
u_cl = zeros(opt_params.num_inputs,opt_params.num_pts);
A_step = zeros(opt_params.num_states,opt_params.num_states,opt_params.num_pts);
B_step = zeros(opt_params.num_states,opt_params.num_inputs,opt_params.num_pts);
c_step = zeros(opt_params.num_states,opt_params.num_pts);

xi_cl(:,1) = opt_params.xi;
dyn_type = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MPC Loop

% Main loop iterating through the length of the trajectory
for i = 1:opt_params.num_pts
    
    fprintf('Iteration: %g\n',i)

    xi_ref = xi_traj(:,i:i+opt_params.horizon_length+1);
    u_ref = u_traj(:,i:i+opt_params.horizon_length);
    xi = sdpvar(opt_params.num_states,opt_params.horizon_length+1);
    u = sdpvar(opt_params.num_inputs,opt_params.horizon_length);
    
    % Linearize dynamics about current state and input
    [A_k, B_k, c_k] = linearize_dynamics_2d_new(opt_params.xi,opt_params.u,opt_params.dt,dyn_type);

    % Linearize dynamics about state and input references
%     [A_k, B_k, c_k] = linearize_dynamics_2d(opt_params.xi,u_traj(:,i),opt_params.dt,dyn_type);
    
    A_step(:,:,i) = A_k;
    B_step(:,:,i) = B_k;
    c_step(:,i) = c_k;
    
    % Solve and save results
    [outputs, opt_flag] = run_yalmip_opt_XZT_one_step(opt_params.horizon_length, ...
        xi, u, xi_ref, u_ref, opt_params.xi, A_k, B_k, c_k, ...
        prev_u, opt_params.spine_params, opt_params.opt_time_lim);
    
    u_step(:,:,i) = outputs{2};
    control = outputs{2}(:,1);
    u_cl(:,i) = control;
    opt_params.u = control;
    prev_u = control;
    
    xi_kp1 = simulate_2d_spine_dynamics(opt_params.xi,opt_params.u,opt_params.dt,1,dyn_type);
    xi_cl(:,i+1) = xi_kp1;
    opt_params.xi = xi_kp1;
        
    AM_step(:,:) = A_step(:,:,i);
    AM_eig(:,i) = eig(AM_step);
    
    [u,s,v] = svd(A_k);
    s_step(:,:,i) = s;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results

% Set up the window.
figure;
hold on;
axis([-0.2, 0.2, -0.1, 0.3]);
% Plot the first location of the spine:
%handles = plot_2d_spine(xi(:,1), two_d_geometry);
handles = plot_2d_tensegrity(xi_cl(:,1), opt_params.spine_params);
% Plot desired state trajectory
plot(xi_traj(1,:),xi_traj(2,:),'or','lineWidth',2)
drawnow;

for i=1:opt_params.num_pts
    fprintf('Plot iteration: %g\n',i)
    % Plot the vertebrae at this timestep:
    % clear the window
    for j = 1:length(handles)
        delete(handles{j});
    end
    % plot:
    %handles = plot_2d_spine(xi(:,i+1), two_d_geometry);
    handles = plot_2d_tensegrity(xi_cl(:,i+1), opt_params.spine_params);
    drawnow;
end

% Plot the state trajectory results of open loop control
figure;
subplot(6,1,1)
plot(1:opt_params.num_pts+1,xi_cl(1,:),'.')
hold on
plot(1:opt_params.num_pts+1,xi_traj(1,1:opt_params.num_pts+1))
ylabel('x')
legend('trajectory','reference','Location','Best')
title('State Trajectories')
subplot(6,1,2)
plot(1:opt_params.num_pts+1,xi_cl(2,:),'.')
hold on
plot(1:opt_params.num_pts+1,xi_traj(2,1:opt_params.num_pts+1))
ylabel('z')
subplot(6,1,3)
plot(1:opt_params.num_pts+1,xi_cl(3,:),'.')
hold on
plot(1:opt_params.num_pts+1,xi_traj(3,1:opt_params.num_pts+1))
ylabel('theta')
subplot(6,1,4)
plot(1:opt_params.num_pts+1,xi_cl(4,:),'.')
hold on
plot(1:opt_params.num_pts+1,xi_traj(4,1:opt_params.num_pts+1))
ylabel('v_x')
subplot(6,1,5)
plot(1:opt_params.num_pts+1,xi_cl(5,:),'.')
hold on
plot(1:opt_params.num_pts+1,xi_traj(5,1:opt_params.num_pts+1))
ylabel('v_z')
subplot(6,1,6)
plot(1:opt_params.num_pts+1,xi_cl(6,:),'.')
hold on
plot(1:opt_params.num_pts+1,xi_traj(6,1:opt_params.num_pts+1))
ylabel('omega')

% Plot the input trajectory results of open loop control
figure;
subplot(4,1,1)
plot(1:opt_params.num_pts,u_cl(1,:),'.')
hold on
plot(1:opt_params.num_pts,u_traj(1,1:opt_params.num_pts))
ylabel('u_1')
legend('trajectory','reference','best')
title('Input Trajectories')
subplot(4,1,2)
plot(1:opt_params.num_pts,u_cl(2,:),'.')
hold on
plot(1:opt_params.num_pts,u_traj(2,1:opt_params.num_pts))
ylabel('u_2')
subplot(4,1,3)
plot(1:opt_params.num_pts,u_cl(3,:),'.')
hold on
plot(1:opt_params.num_pts,u_traj(3,1:opt_params.num_pts))
ylabel('u_3')
subplot(4,1,4)
plot(1:opt_params.num_pts,u_cl(4,:),'.')
hold on
plot(1:opt_params.num_pts,u_traj(4,1:opt_params.num_pts))
ylabel('u_4')

%% Plot the tracking errors of the MPC results

% Plot the absolute errors, square errors and relative errors of each state
% along the trajectory tracking

% Define state-values as the relative-value of states from the start points
x_1 = abs(xi_cl(1,1)*ones(1,size(xi_cl,2))-xi_cl(1,:));
x_2 = abs(xi_cl(2,1)*ones(1,size(xi_cl,2))-xi_cl(2,:));
x_3 = abs(xi_cl(3,1)*ones(1,size(xi_cl,2))-xi_cl(3,:));
 
% Another way is to define state-value as the absolute value of each state
% x_1 = abs(xi_cl(1,:));
% x_2 = abs(xi_cl(2,:));
% x_3 = abs(xi_cl(3,:));

% Calculate absolute errors of states
% x_1 = xi_cl(1,:)-xi_traj(1,1:length(xi_cl(1,:)));
e_1 = xi_cl(1,:)-xi_traj(1,1:end-4);
e_2 = xi_cl(2,:)-xi_traj(2,1:end-4);
e_3 = xi_cl(3,:)-xi_traj(3,1:end-4);

% x_sq_1 = zeros(1,length(x_1));
% x_sq_2 = zeros(1,length(x_1));
% x_sq_3 = zeros(1,length(x_1));
% for i = 1:length(x_1)
%     x_sq_1(i) = x_1(i)^2;
%     x_sq_2(i) = x_2(i)^2;
%     x_sq_3(i) = x_3(i)^2;
% end

% Calculate square errors of states
e_sq_1 = e_1.*e_1;
e_sq_2 = e_2.*e_2;
e_sq_3 = e_3.*e_3;

% e_sq = e_1.*e_1+e_2.*e_2;

% Calculate absolute errors of states
e_abs_1 = abs(e_1);
e_abs_2 = abs(e_2);
e_abs_3 = abs(e_3);
 
% Calculate relative errors of states
% Using absolue state value: cannot reflect absolute value change of
% z_state, as it starts from 0.1 in case of 0 initial sweeping angle
e_rl_1 = e_abs_1./x_1;
e_rl_2 = e_abs_2./x_2;
e_rl_3 = e_abs_3./x_3;

% Plot the relevat errors of each state, namely the absolute errors over
% absolue velues
figure;
subplot(3,1,1)
plot(e_rl_1);
ylabel('RE(x)');
title('Relative Errors of State')
subplot(3,1,2)
plot(e_rl_2);
%ylim([0,0.1]);
ylabel('RE(z)');
subplot(3,1,3)
plot(e_rl_3);
ylabel('RE(theta)')

% % Plot the square errors of each state 
% figure;
% subplot(3,1,1)
% plot(e_sq_1);
% ylabel('SE(x)');
% title('Square Errors of States')
% subplot(3,1,2)
% plot(e_sq_2);
% ylabel('SE(z)');
% subplot(3,1,3)
% plot(e_sq_3);
% ylabel('SE(theta)')

% Plot the absolut error of 3 states, if necessary
figure;
subplot(3,1,1)
plot(e_abs_1);
ylabel('AE(x)');
title('Absolute Errors of States')
subplot(3,1,2)
plot(e_abs_2);
ylabel('AE(z)');
subplot(3,1,3)
plot(e_abs_3);
ylabel('AE(theta)')

% % Plot the eigenvalues of A matrix at each step, if necessary
% figure; 
% % Plot each entries specificly as well, if necessary
% subplot(7,1,1)
% plot(AM_eig)
% legend('lambda 1', 'lambda 2', 'lambda 3', 'lambda 4', 'lambda 5', 'lambda 6')
% title('Eigen Values of Matrix A')
% subplot(7,1,2)
% plot(abs(AM_eig(1,:)))
% subplot(7,1,3)
% plot(abs(AM_eig(2,:)))
% subplot(7,1,4)
% plot(abs(AM_eig(3,:)))
% subplot(7,1,5)
% plot(abs(AM_eig(4,:)))
% subplot(7,1,6)
% plot(abs(AM_eig(5,:)))
% subplot(7,1,7)
% plot(abs(AM_eig(6,:)))
