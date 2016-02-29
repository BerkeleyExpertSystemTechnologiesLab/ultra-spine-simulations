% ultra_spine_mpc.m
% Abishek Akella, Andrew P. Sabelhaus
% This is the primary file for running the ULTRA Spine Model Predictive Control work. Run as a script.

%% Initialize system parameters
%clear variables;
clear all;
close all;
clc;

% The MATLAB functions that calculate the dynamics of the robot must be in
% the current path. If this file is being run outside the
% ultra-spine-simulations repository, change this path.
% Note that the files here are not directly called from this script, but
% instead are used in the following files (or those they call):
%   - simulate_reference_trajectory
%   - simulate_dynamics
%   - linearize_dynamics

path_to_dynamics = '../../dynamics/3d-dynamics-symbolicsolver';
addpath(path_to_dynamics);

% Time step for dynamics
dt = 0.001;

% Flag for adding noise to the forward simulation of the dynamics. noise = 1 turns it on.
noise = 0;

% Parameters for plotting:
% NOTE that the rod sizes here are only used for plotting, 
% since the dynamics used in this work is a point-mass model.

% Load parameters in from the saved file that accompanies the dynamics
spine_geometric_parameters_path = strcat(path_to_dynamics, '/spine_geometric_parameters.mat');
load(spine_geometric_parameters_path);
% Unroll this struct into individual variables. See the dynamics generation script for more information.
% Gravitational force
g = spine_geometric_parameters.g;
% Total number of spine tetrahedrons
N_tetras = spine_geometric_parameters.N;
% Length of one "leg" of the tetrahedron (straight-line distance from center to outer point)
l = spine_geometric_parameters.l;
% Height of one tetrahedtron
h = spine_geometric_parameters.h;
% total mass of one whole tetrahedron
m_t = spine_geometric_parameters.m_t;
% Factor-of-safety with respect to tetrahedron mass (note: this is unused in this script)
FoS = spine_geometric_parameters.FoS;
% mass of one node of the tetrahedron (there are five point masses per tetra)
%m = spine_geometric_parameters.m;

% Radius of a "leg" of the tetrahedron. NOTE that since this is a point-mass model, this parameter is only for visualization.
rad = 0.01;

% Number of links in addition to base link.
% NOTE that this must be consistent with the dynamics defined in
% duct_accel.m and associated files! Those dynamics are pre-calculated,
% and this parameter does NOT change them - it only affects the controller.
links = N_tetras-1; % Here, this is going to be 4-1 = 3.

% Tetrahedron vertical spacing. The initial z-distance between successive tetrahedra
tetra_vertical_spacing = 0.1; % meters

% Projection of leg length (e.g., the x or y coordinate of an endpoint of a node. See spineDynamics for more information.)
leg = (l^2 - (h/2)^2)^.5;

% Number of links in addition to base link.
% NOTE that this must be consistent with the dynamics defined in
% duct_accel.m and associated files! Those dynamics are pre-calculated,
% and this parameter does NOT change them - it only affects the controller.
links = 3;

 % Simulation time
time = 0:dt:500;

% Animation Frame Divisor
frame = 3;
 % Enable string plotting
stringEnable = 1;
% Position of the anchor location of the cables on a tetra node
anchor = [0 0 rad];

% To save a video, uncomment:
% - the following three initialization lines
% - the getframe call within the last loop
% - the open, save, and close lines at the end of this script
%videoObject = VideoWriter('../../videos/ultra-spine-mpc-topbending1.avi');
%videoObject.Quality = 100;
%videoObject.FrameRate = 5;

%% Initialize Plot
Figs = figure('Units','Normalized', 'outerposition', [0 0 1 1]);
M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

cmaps = summer(512);
ax = axes();

grid on;
axis equal;
xlim([-0.25 0.25])
ylim([-0.25 0.25])
zlim([-0.1, 0.5])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
colormap(cmaps(1:256,:))
shading interp
light
lighting phong
hold on
view(3)
title('ULTRA Spine Model')
m = 256;

%% Initialize the simulation
restLengths(1) = 0.1; % Vertical cable rest length
restLengths(2) = 0.1;
restLengths(3) = 0.1;
restLengths(4) = 0.1;
restLengths(5) = 0.187; % Saddle cable rest length
restLengths(6) = 0.187;
restLengths(7) = 0.187;
restLengths(8) = 0.187;

% There are 36 states in this simulation, as it stands: 3 bodies * 12 states each.
systemStates = zeros(links, 12);

% The initial state for all of these 36 variables.
x_initial = [];

% Initialize all the states of the system
% This script currently (2016-02-27) considers "k" to be a different index in different circumstances.
% Here, it's used for the system states as k=1 for the first moving tetrahedron, and k=3 for the topmost one.
% But, for plotting, each of these is shifted up: Tetra{1} and transform{1} are for the first (unmoving) tetra, and Tetra{4} is the topmost.

% Plot the first tetrahedron body (k=1 in the Tetra{} usage)
Tetra{1} = [(l^2 - (h/2)^2)^.5, 0, -h/2; ...
            -(l^2 - (h/2)^2)^.5, 0, -h/2; ...
            0, (l^2 - (h/2)^2)^.5, h/2; ...
            0, -(l^2 - (h/2)^2)^.5, h/2];

% Plot a visualization of this spine tetrahedron
[transform{1}, ~] = plotSpineLink(Tetra{1}, rad, ax);

% Perform 3 iterations: one for each moving tetrahedron.
for k = 1:links
    % Each tetrahedron starts completely still, centered at (x,y) = (0,0) with a z-offset
    x(k) = 0; 
    y(k) = 0.0; 
    z(k) = tetra_vertical_spacing * k; 
    T(k) = 0.0; 
    G(k) = 0.0; 
    P(k) = 0.0;
    dx(k) = 0; 
    dy(k) = 0; 
    dz(k) = 0; 
    dT(k) = 0; 
    dG(k) = 0; 
    dP(k) = 0;
    
    % Save the system states
    % The first state of this tetrahedron is at this starting location
    systemStates(k, :) = [x(k), y(k), z(k), T(k), G(k), P(k), dx(k), dy(k), dz(k), dT(k), dG(k), dP(k)];
    % Append this state to the vector of initial states.
    x_initial = [x_initial; x(k); y(k); z(k); T(k); G(k); P(k); dx(k); dy(k); dz(k); dT(k); dG(k); dP(k)];
    
    % Plot the tetrahedra.
    % Start the initial position of each tetrahedron at the bottom location: center at (0,0,0)
    Tetra{k+1} = [(l^2 - (h/2)^2)^.5, 0, -h/2; ...
                -(l^2 - (h/2)^2)^.5, 0, -h/2; ...
                0, (l^2 - (h/2)^2)^.5, h/2; ...
                0, -(l^2 - (h/2)^2)^.5, h/2];
            
    % Plot a visualization of this spine tetrahedron
    [transform{k+1}, ~] = plotSpineLink(Tetra{k+1}, rad, ax);
    
    % Then, move the tetrahedron body into place by updating the "transform" object.
    % The function below is generated by transforms.m
    % Recall, the system states are indexed starting from 1, but the "Tetra" points are indexed starting with 2 as the bottommost moving vertebra.
    RR{k} =  getHG_Tform(x(k),y(k),z(k),T(k),G(k),P(k));
    % Update the transform object
    set(transform{k+1},'Matrix',RR{k});

    % Finally, move the positions of the cables of this tetrahedra into place, by modifying the Tetra{} array.
    % We use this same transform matrix here.
    % First, append a column of "1"s to this set of coordinates, needed for applying the transform.
    % Note again that there are 4 node points per tetra.
    Tetra{k+1} = [Tetra{k+1}, ones(4,1)]; 
    
    % Move the cables of this tetra into position.
    Tetra{k+1} = RR{k}*Tetra{k+1}';
    % Needs a transpose to return it back to the original orientation
    Tetra{k+1} = Tetra{k+1}';
    % Remove the trailing "1" from the position vectors.
    Tetra{k+1} = Tetra{k+1}(:,1:3);

end

% Plot the cables for this spine position
if (stringEnable)    
    % Get the endpoints of the cables
    String_pts = get_spine_cable_points(Tetra, anchor);
    % Plot. Save the handle so we can delete this set of cables later
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','m');
end

%% Reference Trajectory
% Load in one of the trajectories

%[traj, L] = get_ref_traj_circletop();
%[traj, L] = get_ref_traj_quartercircletop();
%[traj, L] = get_ref_traj_topbending1(); %OPTIMIZATION FAILS AS OF 2016-02-27
%[traj, L] = get_ref_traj_topbending2();
%[traj, L] = get_ref_traj_toprotationtest(); %OPTIMIZATION FAILS AS OF 2016-02-27
[traj, L]  = get_ref_traj_zero();

% Plot this trajectory, for a visualization
plot3(traj(1, :), traj(2,:), traj(3, :), 'r', 'LineWidth', 2);

% Force the figure to draw. The figure at this point includes: tetra bodies, tetra cables, reference trajectory in (x,y,z).
drawnow;

%% Controller Initialization
disp('Controller Initialization')

N = 10; % Horizon length
inputs = sdpvar(repmat(8*links, 1, N-1), repmat(1, 1, N-1));
states = sdpvar(repmat(12*links, 1, N), repmat(1, 1, N));
A_t = sdpvar(repmat(12*links, 1, 12*links), repmat(1, 1, 12*links));
B_t = sdpvar(repmat(12*links, 1, 8*links), repmat(1, 1, 8*links));
c_t = sdpvar(36, 1);
prev_in = sdpvar(8*links, 1);
reference = sdpvar(repmat(12, 1, N), repmat(1, 1, N));

% Create the YALMIP controller for computation of the actual MPC optimizations.
% This function contains all the definitions of the constraints on the optimization, as well as the objective function.
[controller, ~, ~, ~, ~] = get_yalmip_controller_XYZ(N, inputs, states, A_t, B_t, c_t, prev_in, reference);
%[controller, ~, ~, ~, ~] = get_yalmip_controller_XYZT(N, inputs, states, A_t, B_t, c_t, prev_in, reference);

% Build controller object for faster computation during iteration
controller = optimizer(constraints, objective, sdpsettings('solver', 'gurobi', 'verbose', 1), parameters_in, solutions_out);

%% Forward simulate trajectory
disp('Forward simulate trajectory')

[x_ref, u_ref, M] = simulate_reference_traj(controller, systemStates, restLengths, links, dt, x, y, z, T, G, P, ...
    dx, dy, dz, dT, dG, dP, traj, N);

refx = [x_ref{:}];
plot3(refx(25, :), refx(26, :), refx(27, :), 'b-.', 'LineWidth', 2);

%% Build Iterative LQR Controller
disp('Build LQR Controllers for each timestep')

Q = zeros(12);
Q(1:6, 1:6) = eye(6);
Q_lqr = 5*kron(eye(3), Q);
R_lqr = 5*eye(8*links);

tic;
P0 = zeros(36);
[A, B, ~] = linearize_dynamics(x_ref{M}, u_ref{M}, restLengths, links, dt);
K{1} = -((R_lqr + B'*P0*B)^-1)*B'*P0*A;
P_lqr{1} = Q_lqr + K{1}'*R_lqr*K{1} + (A + B*K{1})'*P0*(A + B*K{1});
for k = (M-1):-1:1
    disp(strcat('Controller Build iteration:',num2str(k)))
    [A, B, ~] = linearize_dynamics(x_ref{k}, u_ref{k}, restLengths, links, dt);
    K{M-k+1} = -((R_lqr + B'*P_lqr{M-k}*B)^-1)*B'*P_lqr{M-k}*A;
    P_lqr{M-k+1} = Q_lqr + K{M-k+1}'*R_lqr*K{M-k+1} + (A + B*K{M-k+1})'*P_lqr{M-k}*(A + B*K{M-k+1});
end
toc;

disp('Starting simulation')
%% Perform forward simulation and plot results

% Plot the cables for this spine position
if (stringEnable)    
    % delete the previous plotted cables
    delete(string_handle);
    % Get the endpoints of the cables
    String_pts = get_spine_cable_points(Tetra, anchor);
    % Plot. Save the handle so we can delete these strings later.
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','m');
end

plot_dt = 0.01; % Slow down the animation slightly
offset = 30; % Used for animation so it stays still for a few seconds

% Reset all tetrahedra to their beginning states.
for k = 1:links
    % Each tetrahedron starts out where it was in the beginning. See code above in the initialization section.
    x(k) = 0; 
    y(k) = 0.0; 
    z(k) = tetra_vertical_spacing * k; 
    T(k) = 0.0; 
    G(k) = 0.0; 
    P(k) = 0.0;
    dx(k) = 0; 
    dy(k) = 0; 
    dz(k) = 0; 
    dT(k) = 0; 
    dG(k) = 0; 
    dP(k) = 0;
end

% Loop through each timestep, using the controller and plotting.
s = 1;
for t = 1:((M-1)+offset)
    if mod(t, frame) == 0
        tic
        
        while toc < (plot_dt*frame)
            % Wait until time has passed
        end
        
        % Update the visualization of the tetrahedra.
        % This is done by setting the matrix of "transform{k}" for each tetra.
        for k = 1:links
            % The function below is generated by transforms.m
            RR{k} =  getHG_Tform(x(k),y(k),z(k),T(k),G(k),P(k)); % Build graphical model of each link
            %set(transform{k},'Matrix',RR{k});
            set(transform{k+1},'Matrix',RR{k});
        end
        
        % Plot the cables
        if (stringEnable)
            % First, delete the old cables
            delete(string_handle);
            
            % Calculate the new position of the tetra's coordinates, for plotting, based on the transform from the current system states.
            % This section of code is the same as that in the initialization section, but instead, directly indexes the Tetra{} array.
            for k = 2:(links+1)
                % Reset this specific tetrahedron to the initial state of the bottom tetra.
                Tetra{k} = [(l^2 - (h/2)^2)^.5, 0, -h/2, 1; ...
                            -(l^2 - (h/2)^2)^.5, 0, -h/2, 1; ...
                            0, (l^2 - (h/2)^2)^.5, h/2, 1; ...
                            0, -(l^2 - (h/2)^2)^.5, h/2, 1];
                % Move the coordinates of the string points of this tetra into position. 
                % Note that the transforms are indexed as per the system states: set k=1 is for the first moving tetra,
                % AKA the second tetra graphically.
                Tetra{k} = RR{k-1}*Tetra{k}';
                % Needs a transpose!
                Tetra{k} = Tetra{k}';
                % Remove the trailing "1" from the position vectors.
                Tetra{k} = Tetra{k}(:,1:3);
            end
            
            % Get the coordinates of the spine cables
            String_pts = get_spine_cable_points(Tetra, anchor);
            % Plot the new strings
            string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','m');
        end
        drawnow;
    end
    
    % Update the systemStates matrix with the current states that are saved as individual variables
    for k = 1:links
        systemStates(k, 1) = x(k); 
        systemStates(k, 2) = y(k); 
        systemStates(k, 3) = z(k);
        systemStates(k, 4) = T(k); 
        systemStates(k, 5) = G(k); 
        systemStates(k, 6) = P(k);
        systemStates(k, 7) = dx(k); 
        systemStates(k, 8) = dy(k); 
        systemStates(k, 9) = dz(k);
        systemStates(k, 10) = dT(k); 
        systemStates(k, 11) = dG(k); 
        systemStates(k, 12) = dP(k);
    end
    
    % Forward simulate using the controller

    if ((t > offset) && (t < M + offset))
        % Calculate the input to the system dynamics
        control = K{M+offset-t}*(reshape(systemStates', 36, 1) - x_ref{t-offset}) + u_ref{t-offset};
        % Forward simulate using that control input
        systemStates = simulate_dynamics(systemStates, restLengths, reshape(control, 8, 3)', dt, links, noise);
        % Save the result as the next point in the performed trajectory
        actual_traj(:, s) = systemStates(links, :); s = s + 1;
    end
    
    % Unfoil the new system states back into the series of individual variables
    for k = 1:links
        x(k) = systemStates(k, 1); 
        y(k) = systemStates(k, 2); 
        z(k) = systemStates(k, 3);
        T(k) = systemStates(k, 4); 
        G(k) = systemStates(k, 5); 
        P(k) = systemStates(k, 6);
        dx(k) = systemStates(k, 7); 
        dy(k) = systemStates(k, 8); 
        dz(k) = systemStates(k, 9);
        dT(k) = systemStates(k, 10); 
        dG(k) = systemStates(k, 11); 
        dP(k) = systemStates(k, 12);
    end
    
    % Record this frame for a video
    videoFrames(t) = getframe(gcf);
end

% Plot the resultant trajectory, in full.
plot3(actual_traj(1, :), actual_traj(2, :), actual_traj(3, :), 'g', 'LineWidth', 2);

% Uncomment these lines to save the video.
%open(videoObject);
%writeVideo(videoObject, videoFrames);
%close(videoObject);
