% ultra_spine_mpc_single_simulation.m
% Copyright 2016 Andrew P. Sabelhaus, Abishek Akella, Berkeley Emergent Space Tensegrities Lab
% This function runs a single simulation of MPC.
% All parameters for MPC are passed in, as well as a reference trajectory and a (series of) controllers.

function [mpc_results] = ultra_spine_mpc_single_simulation(traj, controller, optimization_parameters, ...
                            flags, plotting_parameters, paths)                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
%% Function parameters
% Inputs:
%   traj = time series of points in the reference trajectory
%   controller = cell array of YALMIP controllers. There will be one controller for each horizon length,
%       since the horizon shrinks from horizon_length to length 1 during the last few iterations of MPC.
%   optimization_parameters = a struct with various parameters of the spine and of the optimization. See ultra_spine_mpc.
%   plotting_parameters = a struct with various parameters related to the figure and visualization. See ultra_spine_mpc.
%   paths = a struct with strings that are paths to various folder to save information in. As of 2016-04-24, videos and data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function initialization:
% Close all open figure windows, unroll the parameters and paths, intialize saving data and video if specified

% TO-DO: run some asserts here to do error checking.
close all;
function_start = clock;

% Record the current time in a string, used later in the path name for saving data and videos.
start_time_string = datestr(datetime('now'));
% Remove the colons and spaces from this string so that Windows doesn't complain when this repository is cloned
% colons become dashes, spaces become underscores. Regular expressions to the rescue!
start_time_string = regexprep(start_time_string, ':', '-');
start_time_string = regexprep(start_time_string, ' ', '_');

% See ultra_spine_mpc for descriptions of all these variables that are passed in as optimization_parameters.
% The following are unused in this script, except for possibly including for reference in a saved data file with the results: (2016-04-24)
%   optimization_weights (only used for the controller, which was created outside this function)
%   num_points_ref_traj (used to specify traj, which was created outside of this function)
%   direction (used to specify traj, which was created outside of this function)
%   g (already hard-coded into the dynamics script)
%   N_tetras (already hard-coded into the dynamics script)
%   m_t (already hard-coded into the dynamics script)
%   FoS (already hard-coded into the dynamics script)
%   opt_time_limit (only needed for controller)
% Unroll optimization_parameters:

dt = optimization_parameters.dt;

optimization_weights = optimization_parameters.optimization_weights;
links = optimization_parameters.links;
tetra_vertical_spacing = optimization_parameters.tetra_vertical_spacing;
frame = optimization_parameters.frame;
restLengths = optimization_parameters.restLengths;

num_points_ref_traj = optimization_parameters.num_points_ref_traj;
direction = optimization_parameters.direction;
horizon_length = optimization_parameters.horizon_length;
opt_time_limit = optimization_parameters.opt_time_limit;

% Unroll flags:
noise = flags.noise;
save_video = flags.save_video;
save_data = flags.save_data;
stringEnable = flags.stringEnable;
run_lqr = flags.run_lqr;
traj_is_full_system = flags.traj_is_full_system;

spine_geometric_parameters = optimization_parameters.spine_geometric_parameters;
% Unroll this struct into individual variables. See the dynamics generation script for more information.
% g = Gravitational force
% N_tetras = Total number of spine tetrahedrons
% l = Length of one "leg" of the tetrahedron (straight-line distance from center to outer point)
% h = Height of one tetrahedron
% m_t = total mass of one whole tetrahedron
% FoS = Factor-of-safety with respect to tetrahedron mass (note: this is unused in this script)
% m = mass of one node of the tetrahedron (there are five point masses per tetra)
g = spine_geometric_parameters.g;
N_tetras = spine_geometric_parameters.N; %unused as of 2016-04-24
l = spine_geometric_parameters.l;
h = spine_geometric_parameters.h;
m_t = spine_geometric_parameters.m_t;
FoS = spine_geometric_parameters.FoS; % unused as of 2016-04-24
%m = spine_geometric_parameters.m;

% Unroll paths:
path_to_videos_folder = paths.path_to_videos_folder;
path_to_data_folder = paths.path_to_data_folder;

% Unroll plotting_parameters:
time = plotting_parameters.time;
rad = plotting_parameters.rad;
figure_window_location = plotting_parameters.figure_window_location;
figure_window_color = plotting_parameters.figure_window_color;
cmaps = plotting_parameters.cmaps;
figure_rotation = plotting_parameters.figure_rotation;
fontsize = plotting_parameters.fontsize;
cable_color = plotting_parameters.cable_color;
cable_thickness = plotting_parameters.cable_thickness;
trajectory_color = plotting_parameters.trajectory_color;
trajectory_thickness = plotting_parameters.trajectory_thickness;
mpc_result_color = plotting_parameters.mpc_result_color;
mpc_result_thickness = plotting_parameters.mpc_result_thickness;
plot_dt = plotting_parameters.plot_dt;
plotting_offset = plotting_parameters.plotting_offset;
lqr_result_color = plotting_parameters.lqr_result_color;
lqr_result_thickness = plotting_parameters.lqr_result_thickness;
anchor = plotting_parameters.anchor;


% List the names of the variables that this script will save, if save_data is set.
% This is typed in manually.
% TO-DO: automatically generate this list from the 
% TO-DO: update this with new variables
% names of the variables passed in to the function (along with the others that are outputs.)
if(save_data)
    % A cell array of the names of all the variables to save.
    % There must be at least one variable to save.
    variables_to_save = { ...
        'direction', ...
        'dt', ...
        'elapsed_time', ...
        'flags', ...
        'links', ...
        'horizon_length', ...
        'noise', ...
        'num_points_ref_traj', ...
        'optimization_parameters', ...
        'paths', ...
        'plotting_parameters', ...
        'refx', ...
        'restLengths', ...
        'run_lqr', ...
        'spine_geometric_parameters', ...
        'start_time_string', ...
        'tetra_vertical_spacing', ...
        'traj', ...
        'traj_is_full_system', ...
        'u_ref', ...
        'x_initial', ...
        'x_ref'};
end

% Initialize the video, if flagged.
if(save_video)
    %break;
    % create the filename for this video by concatenating with the path to the video folders, defined above
    videoPath = strcat( path_to_videos_folder, 'ultra-spine-mpc_', start_time_string );
    %videoObject = VideoWriter( strcat('../../videos/ultra-spine-mpc_', datestr(datetime('now'))) );
    videoObject = VideoWriter( videoPath, 'Motion JPEG 2000' );
    %videoObject.Quality = 75;
    videoObject.FrameRate = 20;
    videoObject.CompressionRatio = 60;
    %videoObject.VideoCompressionMethod('Motion JPEG 2000');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize plot

% Create the figure window
figure_handle = figure('position', figure_window_location,'Color',figure_window_color);

% TO-DO: figure out what this M variable is and how it's used.
M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

% Configure the figure: colors, orientation, etc.
colormap(cmaps(1:256,:));
ax = axes();
grid on;
axis equal;
hold on;
view(figure_rotation);

% Labels and text
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('ULTRA Spine Model')
 
% Size everything properly
xlim([-0.2 0.2])
ylim([-0.2 0.2])
zlim([-0.1, 0.4])
set(gca,'FontSize',fontsize)

shading interp
light
lighting phong

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the spine at its rest configuration

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
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',cable_thickness,'Color',cable_color);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform some checks on the reference trajectory

% A bit of debugging
disp( 'Reference trajectory has: ');
disp( num2str(size(traj,2)));
disp('timesteps.');



% Plot this trajectory, for a visualization
% If traj is for top tetra only:
if ( ~traj_is_full_system)
    plot3(traj(1, :), traj(2,:), traj(3, :), trajectory_color, 'LineWidth', trajectory_thickness);
else
    % for the full-trajectory version:
    plot3(traj(1, :), traj(2,:), traj(3, :), trajectory_color, 'LineWidth', trajectory_thickness);
    plot3(traj(13, :), traj(14,:), traj(15, :), trajectory_color, 'LineWidth', trajectory_thickness);
    plot3(traj(25, :), traj(26,:), traj(27, :), trajectory_color, 'LineWidth', trajectory_thickness);
end

% Force the figure to draw. The figure at this point includes: tetra bodies, tetra cables, reference trajectory in (x,y,z).
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Model-Predictive Control on the ULTRA Spine 

disp('Starting MPC.')

% This function is now general enough to simulate both 12-state references as well as 36-state references.
[x_ref, u_ref, M] = simulate_mpc_traj(controller, systemStates, restLengths, links, dt, x, y, z, T, G, P, ...
    dx, dy, dz, dT, dG, dP, traj, horizon_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the results of the MPC

disp('Plotting resulting MPC trajectory for vertebrae centers.');

% Transform x_ref from a cell array into a matrix for easy indexing.
refx = [x_ref{:}];

% Check if this script should plot the paths of all 3 vertebrae, or just the top one.
if ( ~traj_is_full_system)
    % Plot top vertebra only
    plot3(refx(25, :), refx(26, :), refx(27, :), mpc_result_color, 'LineWidth', mpc_result_thickness);
else
    % for the full-trajectory version:
    plot3(refx(1, :), refx(2, :), refx(3, :), mpc_result_color, 'LineWidth', mpc_result_thickness);
    plot3(refx(13, :), refx(14, :), refx(15, :), mpc_result_color, 'LineWidth', mpc_result_thickness);
    plot3(refx(25, :), refx(26, :), refx(27, :), mpc_result_color, 'LineWidth', mpc_result_thickness);
end

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build Iterative LQR Controller
% This is now only done if flagged.

if (run_lqr)
    disp('Building LQR Controllers for each timestep.')

    % Create the weighting matrices Q and R
    % Version 1: weight on the position states for each tetrahedron equally (no weight on velocity)
    Q = zeros(12);
    Q(1:6, 1:6) = eye(6);
    Q_lqr = 5*kron(eye(3), Q);
    R_lqr = 5*eye(8*links);

    % Version 2: weight only on the position of the top tetrahedron
    % Q = zeros(36);
    % Q(25:30, 25:30) = 50 * eye(6);
    % Q_lqr = Q;
    % R_lqr = 2*eye(8*links);

    % Version 3: weight on all states
    % Q = 20 * eye(36);
    % Q_lqr = Q;
    % R_lqr = 2*eye(8*links);

    tic;
    P0 = zeros(36);
    % Linearize the dynamics around the final point in the trajectory (timestep M.)
    [A, B, ~] = linearize_dynamics(x_ref{M}, u_ref{M}, restLengths, links, dt);
    % Calculate the gain K for this timestep (at the final timestep, P == the 0 matrix.)
    K{1} = -((R_lqr + B'*P0*B)^-1)*B'*P0*A;
    % Calculate the first matrix P for use in the finite-horizon LQR below
    P_lqr{1} = Q_lqr + K{1}'*R_lqr*K{1} + (A + B*K{1})'*P0*(A + B*K{1});
    % Iterate in creating the finite-horizon LQRs, moving backwards from the end state
    for k = (M-1):-1:1
        disp(strcat('LQR Controller Build iteration:',num2str(k)))
        % Linearize the dynamics around this timestep
        [A, B, ~] = linearize_dynamics(x_ref{k}, u_ref{k}, restLengths, links, dt);
        % Calculate the gain K for this timestep, using the prior step's P
        K{M-k+1} = -((R_lqr + B'*P_lqr{M-k}*B)^-1)*B'*P_lqr{M-k}*A;
        % Calculate the P for this step using the gain K from this step.
        P_lqr{M-k+1} = Q_lqr + K{M-k+1}'*R_lqr*K{M-k+1} + (A + B*K{M-k+1})'*P_lqr{M-k}*(A + B*K{M-k+1});
    end
    toc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform forward simulation and plot results

disp('Plotting vertebrae in motion along trajectory.');
if run_lqr
    disp('run_lqr flag is set, LQR will be run at each step in this motion.');
end

% Plot the cables for this spine position
if (stringEnable)    
    % delete the previous plotted cables
    delete(string_handle);
    % Get the endpoints of the cables
    String_pts = get_spine_cable_points(Tetra, anchor);
    % Plot. Save the handle so we can delete these strings later.
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',cable_thickness,'Color',cable_color);
end

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

% Loop through each timestep, using the LQR controller if flag is set, and plotting the tetrahedra and their cables.
s = 1;
for t = 1:((M-1)+plotting_offset)
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
            string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',cable_thickness,'Color',cable_color);
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
    
    % If we are just doing straight MPC, replay the state trajectory directly.
    if (~run_lqr)
        if ((t > plotting_offset) && (t < M + plotting_offset))
            % just update systemStates to x_ref at t, since x_ref here is the output of MPC.
            systemStates = reshape(x_ref{t-plotting_offset}, 12, 3)';
            % Save the result as the next point in the performed trajectory.
            actual_traj(:, s) = systemStates(links, :);
            s = s + 1;
        end
    else
        % Run the LQR controllers.
        % Forward simulate using the controller
        if ((t > plotting_offset) && (t < M + plotting_offset))
            % Calculate the input to the system dynamics
            % A general representation of u(t) = K(t) * (x(t) - x_{ref}(t)) + u_{ref}(t)
            control = K{M+plotting_offset-t}*(reshape(systemStates', 36, 1) - x_ref{t-plotting_offset}) + u_ref{t-plotting_offset};
            % Forward simulate using that control input
            systemStates = simulate_dynamics(systemStates, restLengths, reshape(control, 8, 3)', dt, links, noise);
            % Save the result as the next point in the performed trajectory
            actual_traj(:, s) = systemStates(links, :); 
            % Save the control results for examination later
            actual_control_inputs(:,s) = control;
            s = s + 1;
        end
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
    % drawnow is needed here or else the MPC script with multiple iterations of MPC will destroy the video.
    drawnow;
    videoFrames(t) = getframe(gcf);
end

% Plot the resultant trajectory, in full.
if (run_lqr)
    plot3(actual_traj(1, :), actual_traj(2, :), actual_traj(3, :), lqr_result_color, 'LineWidth', lqr_result_thickness);
else
    % just running MPC
    % plot the x trajectory directly. Let's just do the top element:
    %plot3(x_ref{
    % this is going to require manipulating data. do it later.
    %plot3(actual_traj(1,:), actual_traj(2,:), actual_traj(3,:), 'g', 'LineWidth', 2);
end

videoFrames(t+1) = getframe(gcf);

% Save the last frame, now with full trajectory plotted.
% Save a few of them in a row for a better visualization.
% for frames = 1:5
%     videoFrames(t + frame) = getframe(gcf)
% end

% Record the time now, at the end of the script
function_end = clock;
% Record the difference, in minutes.
elapsed_time = etime(function_end, function_start) / 60;


% Save the data from this simulation, if the save_data flag is set.
if(save_data)
    % create the filename for this data by concatenating with the path to the data folder, defined above
    data_path = strcat( path_to_data_folder, 'ultra-spine-mpc_data_', start_time_string );
    % Since MATLAB doesn't seem to have a way to save a specific list of variables, this script
    % first creates a .mat file and then appends variables to it.
    save(data_path, variables_to_save{1});
    % Then, iterate over all the other variable names.
    for i=2:size(variables_to_save,2)
        % save each variable
        save(data_path, variables_to_save{i}, '-append');
    end
end

% Save the video, if the save_video flag is set
if(save_video)
    open(videoObject);
    writeVideo(videoObject, videoFrames);
    close(videoObject);
end

% Set the mpc_results cell array to have all the variables to save.
mpc_results = cell( size(variables_to_save,2), 1);
for i = 1:size(variables_to_save,2)
    % Store the designated variable into this cell array, with name-value pairs
    mpc_results{i} = { variables_to_save{i}, eval(variables_to_save{i}) };
end

% End function.







