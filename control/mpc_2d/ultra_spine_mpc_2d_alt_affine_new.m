% ultra_spine_mpc_2d_alt_affine_new.m
%
% Copyright 2016-2018 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This is the primary script file that runs MPC on a 2d spine
% After much experimentation, this implementation of MPC is what's going 
% to be submitted for publication. -Drew 2018-01-27

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
% inv_kin_path = '../../kinematics/2D-inv-kinematics-new';
addpath(inv_kin_path)
% Reference trajectories path
ref_traj_path = './reference_trajectories';
addpath(ref_traj_path)
% YALMIP controllers path
yalmip_controllers_path = './yalmip_controllers';
addpath(yalmip_controllers_path)
% For plotting the vertebra:
plotting_path = './plotting';
addpath(plotting_path);

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
% To-do: try and keep the same size of an A matrix, but with the vertebrae
% going slower. maybe if we do num_pts = 800, then we could get something
% better out of dt=1e-3?
% Also, should try out the better inverse kinematics solution here to see
% if it makes a difference. Maybe if we keep the robot in 
% NOTE that this script really does num_pts+1 points, since includes
% initial position, so you'd set 399 for 400 points, for example.
%opt_params.num_pts = 399;
% testing the visualization:
opt_params.num_pts = 9;
opt_params.num_states = 6;
opt_params.num_inputs = 4;
opt_params.horizon_length = 4;
opt_params.opt_time_lim = 1.5;
opt_params.spine_params = two_d_geometry;
opt_params.dt = 0.01;
% worked decent with 1e-5.
opt_params.dt = 1e-5;
%opt_params.dt = 1e-4;
%opt_params.dt = 1e-3;
% opt_params.dt = 0.1/opt_params.num_pts;

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
[xi_traj, ~] = get_ref_traj_invkin_XZG_new(0.1,opt_params.num_pts+opt_params.horizon_length+1,-1,opt_params.dt);
u_traj = zeros(opt_params.num_inputs,opt_params.num_pts+opt_params.horizon_length+1);

% Use the inverse kinematics for the 2D spine to generate reference inputs.
% To-do: parameterize the pretension force? That's the third input here.
min_cable_tension = 5; % N, I think? Depends on units in inv kin.
% was 30 for the working version.

for i = 1:opt_params.num_pts+opt_params.horizon_length+1
    [~, u_traj(:,i)] = getTensions(xi_traj(:,i), opt_params.spine_params, ...
                            min_cable_tension);
%     disp(xi_traj(:,i))
end
opt_params.xi = xi_traj(:,1);
opt_params.xip1 = xi_traj(:,2);
% opt_params.xi(:,1:2) = xi_traj(:,1:2);



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
xi_clp1(:,1) = opt_params.xip1;

% For the forward simulation of the nonlinear dynamics, this script
% supports multiple types of simulations. They're types 1 through 4.
% See simulate_2d_spine_dynamics for more information. Here's a brief
% overview:
%   1 = no constraint on negative cable tensions, dynamics model allows
%       cables to exert negative tensions (not realistic.)
%   2 = rectified cable tensions, negative reset to zero. (realistic.)
%   3 = same as approach 1 but uses a different function. This option will
%       be used to compare the forward-simulation itself, and it not really
%       useful in the context of this MPC script.
%   4 = uses the logistically-smoothed cable tensions model, which
%       rectifies negative tensions back to 0, but is smooth!

% Using dynamics type 2 corresponds with the approach in the 2015 ACC
% paper.
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
    % IS THIS LINEARIZING AROUND THE ACTUAL TRAJECTORY OR THE REFERENCE???
    [A_k, B_k, c_k] = linearize_dynamics_2d_affine(opt_params.xi,opt_params.xip1, ...
        opt_params.u,opt_params.dt,dyn_type);
 %   [A_k, B_k, c_k] = linearize_dynamics_2d_affine(xi_traj(:,i),u_traj(:,i),opt_params.dt,dyn_type);

    xi_step(:,i) = opt_params.xi;
    
    % Linearize dynamics about state and input references
%     [A_k, B_k, c_k] = linearize_dynamics_2d(opt_params.xi,u_traj(:,i),opt_params.dt,dyn_type);
    
    A_step(:,:,i) = A_k;
%     B_step(:,:,i) = B_k;
%     c_step(:,i) = c_k;
    
    % Solve and save results
    [outputs, opt_flag] = run_yalmip_opt_XZT_one_step(opt_params.horizon_length, ...
        xi, u, xi_ref, u_ref, opt_params.xi, A_k, B_k, c_k, ...
        prev_u, opt_params.spine_params, opt_params.opt_time_lim);
    
    u_step(:,:,i) = outputs{2};
    control = outputs{2}(:,1);
    u_cl(:,i) = control;
    opt_params.u = control;
    prev_u = control;
    
%     xi_kp1 =    simulate_2d_spine_dynamics(opt_params.xi,opt_params.u,opt_params.dt,1,dyn_type);
    xi_kp1 = simulate_2d_spine_dynamics(opt_params.xi,opt_params.u,opt_params.dt,1,dyn_type);
    xi_kp2 = simulate_2d_spine_dynamics(opt_params.xip1,opt_params.u,opt_params.dt,1,dyn_type);
    
%     xi_kp1 = simulate_2d_spine_dynamics(opt_params.xi,opt_params.u,opt_params.dt,1,dyn_type)... 
%         -A_k*opt_params.xi - B_k*opt_params.u;
    xi_cl(:,i+1) = xi_kp1;
    opt_params.xi = xi_kp1;
    
    opt_params.xip1 = xi_kp2;
        
    AM_step(:,:) = A_step(:,:,i);
    AM_eig(:,i) = eig(AM_step);
    
    [u,s,v] = svd(A_k);
    s_step(:,:,i) = s;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results

% To-do: make the lines as the 3D MPC, and put the trajectory line "on top"
% of the cylinders so it's actually visible.
% Also to-do: make the surfaces a single color, since gray is NOT working
% well here.
% From the internet: something like
% h = surf(peaks(50));
%set(h,'edgecolor','none','facecolor',[.1 .9 .1])
% ...which sets to a green color. What's black, in a 3-vector?


% Specify the plot colors and line types
cable_color = 'r';
cable_thickness = 2;
trajectory_color = 'b';
% needs to be a bit larger here to be visible.
% Note: Seems to be weirdly inconsistent between dots and lines? The dots
% are big when zoomed out but small when zoomed in...
trajectory_thickness = 2;
mpc_result_color = 'c';
mpc_result_thickness = 2;

% Set up the window.
figure;
hold on;
axis equal;
% used to be: [-0.2, 0.2, -0.1, 0.3]
axis([-0.2, 0.2, -0.1, 0.2]);
% coloring:
cmaps = gray(512); 
colormap(cmaps(1:256,:));
% some shading settings:
shading interp
light
lighting phong

% Labels and text
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('Spine Model')

% Size everything properly
%xlim([-0.2 0.2])
%ylim([-0.2 0.2]) % Note, we use 'axis' above instead
%zlim([-0.1, 0.4])
fontsize = 20;
set(gca,'FontSize',fontsize)

% Plot the first location of the spine:
%handles = plot_2d_spine(xi(:,1), two_d_geometry);
%handles = plot_2d_tensegrity(xi_cl(:,1), opt_params.spine_params);
handles = plot_2d_tensegrity_surfaced(xi_cl(:,1), opt_params.spine_params, gca);

% For all the surfaces: change them to black only.
%black = [0,0,0];
%for i=1:length(handles)
%    % Failure here: should only set the surfaces to black! doesn't work for
%    % handles to lines.
%    % Also need to turn off all shading and lighting.
%    set(handles{i},'edgecolor','none','facecolor',black);
%end

% Plot desired state trajectory
% In order to make the line appear on top of the surf, offset it by some
% amount.
traj_line_offset = 0.1;
% Need to make it into a vector. Number of points is 2nd dimension
% (columns) of xi_cl.
% Number of actual points in closed loop:
pts_cl = size(xi_cl, 2);
line_z = traj_line_offset * ones(1, pts_cl);
%plot(xi_traj(1,:),xi_traj(2,:),'or','lineWidth',2)
% Pulling from the xi_traj, up until endpoint.
%plot3(xi_traj(1,1:pts_cl), xi_traj(2,1:pts_cl), line_z, '.', trajectory_color, 'LineWidth', trajectory_thickness);
% A very complicated plot command, which makes circles at each position of
% the center of mass, and draws a line between them:
plot3(xi_traj(1,1:pts_cl), xi_traj(2,1:pts_cl), line_z, '-o', ...
    'MarkerEdgeColor', trajectory_color, 'MarkerFaceColor', ...
    trajectory_color, 'MarkerSize', trajectory_thickness, ...
    'LineWidth', trajectory_thickness);
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
    %handles = plot_2d_tensegrity(xi_cl(:,i+1), opt_params.spine_params);
    handles = plot_2d_tensegrity_surfaced(xi_cl(:,i+1), opt_params.spine_params, gca);
    drawnow;
end

% Plot the trajectory of the centers of mass of the moving vertebra.
% Was, in 3D:
%plot3(refx(1, :), refx(2, :), refx(3, :), mpc_result_color, 'LineWidth', mpc_result_thickness);
%plot3(xi_cl(1,:), xi_cl(2,:), line_z, '.', mpc_result_color, 'LineWidth', mpc_result_thickness);
plot3(xi_cl(1,:), xi_cl(2,:), line_z, '-o', ...
    'MarkerEdgeColor', mpc_result_color, 'MarkerFaceColor', ...
    mpc_result_color, 'MarkerSize', mpc_result_thickness, ...
    'LineWidth', mpc_result_thickness);

% Try to turn off all lights?
% This seems to work!
delete(findall(gcf,'Type','light'));

%% Plot the state trajectory results of open loop control
figure;
subplot(3,1,1)
plot(1:opt_params.num_pts+1,100*xi_traj(1,1:opt_params.num_pts+1),'-b','LineWidth',1.5)
hold on
plot(1:opt_params.num_pts+1,100*xi_cl(1,:),'-xg','LineWidth',1.5)
ylabel('x /cm')
legend('reference','traj with u_{ref}','Location','Best')
title('State Trajectories')
grid
subplot(3,1,2)
plot(1:opt_params.num_pts+1,100*xi_traj(2,1:opt_params.num_pts+1),'-b','LineWidth',1.5)
hold on
plot(1:opt_params.num_pts+1,100*xi_cl(2,:),'-xg','LineWidth',1.5)
ylabel('z /cm')
grid
subplot(3,1,3)
plot(1:opt_params.num_pts+1,180/pi*xi_traj(3,1:opt_params.num_pts+1),'-b','LineWidth',1.5)
hold on
plot(1:opt_params.num_pts+1,180/pi*xi_cl(3,:),'-xg','LineWidth',1.5)
ylabel('\theta /��')
grid

%% Plot the input trajectory results of open loop control
figure;
subplot(4,1,1)
plot(1:opt_params.num_pts,u_traj(1,1:opt_params.num_pts),'b','LineWidth',1.5)
hold on
plot(1:opt_params.num_pts,u_cl(1,:),'-xg','LineWidth',1.5)
grid;
ylabel('u_1 /m')
legend('reference','trajectory','Location','Best')
title('Input Trajectories')
subplot(4,1,2)
plot(1:opt_params.num_pts,u_traj(2,1:opt_params.num_pts),'b','LineWidth',1.5)
hold on
plot(1:opt_params.num_pts,u_cl(2,:),'-xg','LineWidth',1.5)
ylabel('u_2 /m')
grid;
subplot(4,1,3)
plot(1:opt_params.num_pts,u_traj(3,1:opt_params.num_pts),'b','LineWidth',1.5)
hold on
plot(1:opt_params.num_pts,u_cl(3,:),'-xg','LineWidth',1.5)
grid;
ylabel('u_3 /m')
subplot(4,1,4)
plot(1:opt_params.num_pts,u_traj(4,1:opt_params.num_pts),'b','LineWidth',1.5)
hold on
plot(1:opt_params.num_pts,u_cl(4,:),'-xg','LineWidth',1.5)
grid;
ylabel('u_4 /m')

%% Analyse the errors of the MPC results
% Define the distance between the start point and end point of the
% reference as x_1, x_2, x_3:
x_1 = abs(xi_traj(1,1)*ones(1,size(xi_cl,2))-xi_traj(1,1:(end-opt_params.horizon_length)));
x_2 = abs(xi_traj(2,1)*ones(1,size(xi_cl,2))-xi_traj(2,1:(end-opt_params.horizon_length)));
x_3 = abs(xi_traj(3,1)*ones(1,size(xi_cl,2))-xi_traj(3,1:(end-opt_params.horizon_length)));

% % Prior Version: Use absolue state value as /x, which cannot reflect absolute value change of
% z_state, as it starts from 0.1 in case of 0 initial sweeping angle
% x_1 = abs(xi_cl(1,:));
% x_2 = abs(xi_cl(2,:));
% x_3 = abs(xi_cl(3,:));

% Calculate absolute errors of states
% x_1 = xi_cl(1,:)-xi_traj(1,1:length(xi_cl(1,:)));
e_1 = xi_cl(1,:)-xi_traj(1,1:end-opt_params.horizon_length);
e_2 = xi_cl(2,:)-xi_traj(2,1:end-opt_params.horizon_length);
e_3 = xi_cl(3,:)-xi_traj(3,1:end-opt_params.horizon_length);


% % TEST for trajectory lag:
% e_1p = xi_cl(1,:)/0.9 - xi_traj(1,1:end-opt_params.horizon_length);
% e_2p = (xi_cl(2,1)*ones(1,size(xi_cl,2)) - (xi_cl(2,1)*ones(1,size(xi_cl,2)) - xi_cl(2,:))/0.8) - xi_traj(2,1:end-opt_params.horizon_length);
% e_3p = xi_cl(3,:)/0.9 - xi_traj(3,1:end-opt_params.horizon_length);

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
e_rl_1 = e_abs_1./x_1;
e_rl_2 = e_abs_2./x_2;
e_rl_3 = e_abs_3./x_3;


% Plot the relevat errors of each state, namely the absolute errors over absolue velues
figure;
subplot(3,1,1)
plot(100*e_rl_1,'.-','LineWidth',1.25);
grid;
ylabel('Relative E(x)/ %');
title('Relative Tracking Errors of States')
subplot(3,1,2)
plot(100*e_rl_2,'.-','LineWidth',1.25);
ylabel('Relative E(z)/ %');
grid;
subplot(3,1,3)
plot(100*e_rl_3,'.-','LineWidth',1.25);
grid;
ylabel('Relative E(\theta)/ %')

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
plot(100*e_abs_1,'.-','LineWidth',1.25);
grid;
ylabel('AE(x)/cm');
title('Absolute Errors of States')
subplot(3,1,2)
plot(100*e_abs_2,'.-','LineWidth',1.25);
grid;
ylabel('AE(z)/cm');
subplot(3,1,3)
plot(180/pi*e_abs_3,'.-','LineWidth',1.25);
grid;
ylabel('AE(\theta)/��');
xlabel('steps');

%% PLOT THE X-Z POSITION
figure;
plot(100*xi_traj(1,1:opt_params.num_pts+1),100*xi_traj(2,1:opt_params.num_pts+1),'-ob','LineWidth',2.5);
% plot(100*xi_traj(1,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10),100*xi_traj(2,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10),'-o','LineWidth',2.5);
% plot(100*xi_traj(1,1:opt_params.num_pts+2-2*(opt_params.num_pts+1)/10),100*xi_traj(2,1:opt_params.num_pts+2-2*(opt_params.num_pts+1)/10),'LineWidth',2.5);
hold on;
plot(100*xi_cl(1,:), 100*xi_cl(2,:),'-xr','LineWidth',2.5);
grid on;
hold on;
plot(100*xi_cl(1,1), 100*xi_cl(2,1),'oy','LineWidth',3.5);
xlabel('X /cm');
ylabel('Z /cm');
% ylim([0,10]);
title('Plot of X-Z Position');
legend('reference','tractory','start point','location','best');
%% PLOT THE ANGULAR TRAJECTORY
figure;
plot(xi_traj(3,1:opt_params.num_pts+1)*180/pi,'-ob','LineWidth',2.5);
% plot(xi_traj(3,1:opt_params.num_pts+2-(opt_params.num_pts+1)/10)*180/pi,'-o','LineWidth',2.5);
% plot(xi_traj(3,1:opt_params.num_pts+1-1.5*(opt_params.num_pts+1)/10)*180/pi,'LineWidth',2.5);
% plot(xi_cl(3,:), '-x','LineWidth',2.5);
hold on;
plot(xi_cl(3,:)*180/pi, '-xr','LineWidth',2.5);
grid on;
hold on;
% plot(xi_cl(3,1),'o','LineWidth',3.5);
plot(xi_cl(3,1)*180/pi,'oy','LineWidth',3.5);
xlabel('steps');
% ylabel(' \theta /arc');
ylabel(' \theta /��');
title('Plot of \theta over steps');
legend('reference','tractory','start point','location','best');
%% PLOT X-Y-Theta IN 3D
figure;
plot3(100*xi_cl(1,:), 100*xi_cl(2,:),180/pi*xi_cl(3,:),'-bx','LineWidth',2.5);
hold on;
plot3(100*xi_traj(1,1:opt_params.num_pts+1), ...
    100*xi_traj(2,1:opt_params.num_pts+1), ...
    180/pi*xi_traj(3,1:opt_params.num_pts+1), ...
    '-ro','LineWidth',1.5);
% plot3(100*xi_traj(1,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10), ...
%     100*xi_traj(2,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10), ...
%     180/pi*xi_traj(3,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10), ...
%     '-o','LineWidth',2.5);
grid on;
hold on;
plot3(100*xi_cl(1,1), 100*xi_cl(2,1),180/pi*xi_cl(3,1),'yo','LineWidth',3.5);
xlabel('X /cm');
ylabel('Z /cm');
zlabel(' \theta /��');
title('Plot of states X - Z - \theta');
legend('tractory','reference','start point','location','best');
