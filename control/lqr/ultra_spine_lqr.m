% ultra_spine_lqr.m
% Andrew P. Sabelhaus, Abishek Akella
% This is the primary file for running time-varying LQR on the ULTRA Spine model. Run as a script.

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
% Note that this is the time between two successive points in the trajectory.
% So, the total time for this simulation is dt * length(trajectory).
dt = 1e-4;

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

% To save a video, set this flag to 1, and change the name of the output file.
save_video = 1;

if(save_video)
    %break;
    videoObject = VideoWriter( strcat('../../videos/ultra-spine-lqr_', datestr(datetime('now'))) );
    videoObject.Quality = 90;
    videoObject.FrameRate = 5;
end



%% Initialize Plot
% Note that this is performed at the beginning so the visualization of the terahedra bodies can be loaded properly.

% Create the figure window
%figure_handle = figure('position', [100, 100, 700 800],'Color','w');
figure_handle = figure('position', [0, 0, 600 700],'Color','w');

M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

% Set the color map
cmaps = summer(512);
colormap(cmaps(1:256,:))
%colormap(cool);
ax = axes();

grid on;
axis equal;

hold on;

% Rotate for a better visualization
view([-20, 14]);

% Labels
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('ULTRA Spine Model')

% Size everything properly
xlim([-0.2 0.2])
ylim([-0.2 0.2])
zlim([-0.1, 0.4])
set(gca,'FontSize',24)

shading interp
light
lighting phong

%m = 256;

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
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','r');
end

%% Reference Trajectory
% Load in one of the trajectories

% TO-DO 2016-04-18: change these to be trajectories in the full state space!
% For example, just something kinematic with all the tetrahedra.

%[traj, ~] = get_ref_traj_circletop();
%[traj, ~] = get_ref_traj_quartercircletop();
%[traj, ~] = get_ref_traj_topbending1(); % Has trajectories along angles. NOT WORKING WELL as of 2016-02-28...
%[traj, ~] = get_ref_traj_topbending2();
[traj, ~] = get_ref_traj_topbending_YZ();
%[traj, ~] = get_ref_traj_toprotationtest(); % Has trajectories along angles.
%[traj, ~]  = get_ref_traj_zero();

% Plot this trajectory, for a visualization
plot3(traj(1, :), traj(2,:), traj(3, :), 'b', 'LineWidth', 2);

% On 2016-04-18:
% Since we're only using a weighting matrix Q for the positions of the top tetrahedron, it *should* be OK to just pad the
% trajectory with zeroes on all the other states. After all, we're not tracking those errors!
% Pad with zeros for the middle two (links - 1) tetrahedra.
full_traj = [ zeros( (links - 1)*12, size(traj,2)) ; traj];

% Force the figure to draw. The figure at this point includes: tetra bodies, tetra cables, reference trajectory in (x,y,z).
drawnow;

%% Initialize the first controller, and get the first input u(0)

% The controller we're solving for optimizes the object J = (1/2) \sum_{k_0}^{\inf} [x(k)' Q x(k) + u(k)' R u(k)]
% Define our weighting matrices for the infinite-horizon LQR cost function.

% For Q:
% Version 1: weight on the position states for each tetrahedron equally (no weight on velocity)
Q_single_tetra = zeros(12);
Q_single_tetra(1:6, 1:6) = eye(6);
% Expand to a big enough Q for all tetras.
Q = 5*kron(eye(3), Q_single_tetra);

% Version 2: weight only on the position of the top tetrahedron
%Q = zeros(36);
%Q(25:30, 25:30) = 5 * eye(6); % I also tried to weight it as 50 * eye before...

% Version 3: weight on all states
% Q = eye(36); % Drew also tried 20 * I at some point before...

% For R:
R = eye(8*links); % used to be 5*eye

% reshape the states
% Actually just use x_initial here, it's the same as reshape(systemStates', 36, 1).

% The inputs here will be zero to start
% TO-DO: change this to have a better start. Use the cable lengths from the inverse kinematics script(s).
% There are 8*links inputs. 2016-04-18 that's 24 cables. In this simulation, the inputs are changes to rest length.
% NOTE that u_initial is ONLY used for the initial linearization: it's never actually applied to the system.
% It's really sort of u_(-1).
u_initial = zeros(8*links,1);

% All trajectories start at x_initial. So, there will be size(traj,2)-1 number of linearizations (we won't need a control input for the final state.)
%A_t = zeros(12*links, 12*links, size(traj,2)-1);
% Similarly, construct B, with 8 inputs :
%B_t = zeros(12*links, 8*

% For now, just use a variable that will change at each iteration, don't record.
[A_t, B_t, ~] = linearize_dynamics(x_initial, u_initial, restLengths, links, dt);

% Use MATLAB's built in infinite-horizon LQR solver to calculate our gain matrices.
%[K_t, P_t, ~] = dlqr(A_t, B_t, Q, R);

% Since MATLAB's dlqr throw errors when the system is unstable (we want to still try and run a control if that's the case, since stability
% may change with successive linearizations), let's try and use MATLAB's discrete algebrai riccati equation solver directly.
% The DARE is A'PA - P + Q = A'PB[R + B'PB]^(-1) * B'PA
[P_t, ~, G_t] = dare(A_t, B_t, Q, R);

% The control input can then be calculated from P_t. G_t is the proportional gain matrix, G = [R + B'PB]^(-1) * B'PA
% Since we're doing trajectory tracking, the first input is then
%u_t = K_t * (x_initial - full_traj(:,1));

% Let's just try regulation for now, trajectory tracking later.
u_t = -G_t * x_initial;

% well, u_1 should be zero. TO-DO: justify this?

% and t+1 is (i.e., timestep 1, the second in the trajectory, since our notation starts at zero but MATLAB doesn't)
%x_tp1 = simulate_dynamics(x_initial, restLengths, reshape(u_t, 8, 3)', dt, links, noise);
% Has to be in the form of:
x_tp1 = simulate_dynamics(reshape(x_initial, 12, links)', restLengths, reshape(u_t, 8, links)', dt, links, noise);

% Reshape x_tp1 to the state vector we require
x_tp1 = reshape(x_tp1',36,1);

% Save the result as the next point in the performed trajectory
% NOTE: This will be one state behind full_traj. E.g, the first column of actual_traj will be state 1, while the first colunm of traj is state zero.
% Actually, just augment it now, that's easier. Now it's on the same numbering.
actual_traj = x_initial;
actual_traj(:, 2) = x_tp1; 
% Save the control results for examination later
actual_control_inputs = u_t;

% There are num_timesteps points in our trajectory:
num_timesteps = size(full_traj,2);

%% Then, run LQR on the rest of the trajectory

for k = 2:num_timesteps-1 % Start from the second point in the trajectory, and go up until the end: then x_tp1 will be state num_timesteps.
    
    disp(k);
    
    % This is for timestep tp1. Calculate the input u_tp1.
    % Re-linearize. Use last timestep's input. TO-DO: formalize this use of input.
    [A_t, B_t, ~] = linearize_dynamics(x_tp1, u_t, restLengths, links, dt);
    % Calculate this step's LQR constants
    %[K_t, P_t, ~] = dlqr(A_t, B_t, Q, R);
    % Like above, try out matlab's DARE solver instead:
    [P_t, ~, G_t] = dare(A_t, B_t, Q, R);
    
    % Then, the u_t+1 input will be:
    %u_tp1 = K_t * (x_tp1 - full_traj(:,k));
    % Like above, let's just try regulation.
    u_tp1 = -G_t* x_tp1;
    
    % Update our variables, since now we have the tp1 timestep for both x and u.
    x_t = x_tp1;
    u_t = u_tp1;
    
    % Apply this control input and get a new x_tp1.
    x_tp1 = simulate_dynamics(reshape(x_t, 12, links)', restLengths, reshape(u_t, 8, links)', dt, links, noise);
    
    % (For plotting later): Unfoil the new system states back into the series of individual variables
    for i = 1:links
        x(i) = x_tp1(i, 1); 
        y(i) = x_tp1(i, 2); 
        z(i) = x_tp1(i, 3);
        T(i) = x_tp1(i, 4); 
        G(i) = x_tp1(i, 5); 
        P(i) = x_tp1(i, 6);
        dx(i) = x_tp1(i, 7); 
        dy(i) = x_tp1(i, 8); 
        dz(i) = x_tp1(i, 9);
        dT(i) = x_tp1(i, 10); 
        dG(i) = x_tp1(i, 11); 
        dP(i) = x_tp1(i, 12);
    end
    
    % Reshape x_tp1 to the state vector we require
    x_tp1 = reshape(x_tp1',36,1);
    
    % Save these values. Remember that at the end of this loop, k+1 will be the final timestep, and x_tp1 will be x_{num_timesteps}.
    actual_traj(:,k+1) = x_tp1;
    actual_control_inputs(:,k) = u_t;
    
    % Plot:
    
    % Update the visualization of the tetrahedra.
    % This is done by setting the matrix of "transform{k}" for each tetra.
    for i = 1:links
        % The function below is generated by transforms.m
        RR{i} =  getHG_Tform(x(i),y(i),z(i),T(i),G(i),P(i)); % Build graphical model of each link
        %set(transform{i},'Matrix',RR{i});
        set(transform{i+1},'Matrix',RR{i});
    end
    
    % First, delete the old cables
    delete(string_handle);

    % Calculate the new position of the tetra's coordinates, for plotting, based on the transform from the current system states.
    % This section of code is the same as that in the initialization section, but instead, directly indexes the Tetra{} array.
    for i = 2:(links+1)
        % Reset this specific tetrahedron to the initial state of the bottom tetra.
        Tetra{i} = [(l^2 - (h/2)^2)^.5, 0, -h/2, 1; ...
                    -(l^2 - (h/2)^2)^.5, 0, -h/2, 1; ...
                    0, (l^2 - (h/2)^2)^.5, h/2, 1; ...
                    0, -(l^2 - (h/2)^2)^.5, h/2, 1];
        % Move the coordinates of the string points of this tetra into position. 
        % Note that the transforms are indexed as per the system states: set k=1 is for the first moving tetra,
        % AKA the second tetra graphically.
        Tetra{i} = RR{i-1}*Tetra{i}';
        % Needs a transpose!
        Tetra{i} = Tetra{i}';
        % Remove the trailing "1" from the position vectors.
        Tetra{i} = Tetra{i}(:,1:3);
    end

    % Get the coordinates of the spine cables
    String_pts = get_spine_cable_points(Tetra, anchor);
    % Plot the new strings
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','r');
    
    drawnow;
    
    % Record this frame for a video
    videoFrames(k-1) = getframe(gcf);
    
end

% Plot the resultant trajectory, in full.
plot3(actual_traj(1, :), actual_traj(2, :), actual_traj(3, :), 'g', 'LineWidth', 2);

% Save the video, if the save_video flag is set
if(save_video)
    open(videoObject);
    writeVideo(videoObject, videoFrames);
    close(videoObject);
end

%% End script.

