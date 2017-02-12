% ultra_spine_mpc.m
% Copyright 2016 Andrew P. Sabelhaus, Abishek Akella, Berkeley Emergent Space Tensegrities Lab
% This is the primary file for running the ULTRA Spine Model Predictive Control work.
% This script calls one or more ultra_spine_mpc_single_simulation runs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script outline
% This script can run multiple iterations of MPC at a time. The general outline is:
% 1) Define the default set of parameters
% 2) Make a list of the trajectories, controllers, and parameters to run (based off default)
% 3) Loop through the following procedure: from 1 to end of list of runs,
%   3.a) Load in the specified trajectory
%   3.b) Create the specified controller
%   3.c) Run MPC 
% 4) to-do: automatic data analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script initialization
%clear variables
clear all;
close all;
clc;

disp('ULTRA Spine MPC');

% Record the start time of this script:
script_start = clock;

% Add various paths

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

% Path to the .mat file that holds the spine_geometric_parameters struct, as created by 
% the spineDynamics script in that folder.
spine_geometric_parameters_path = strcat(path_to_dynamics, '/spine_geometric_parameters.mat');

% Reference trajectories and controllers are now in subfolders:
path_to_reference_trajectories = './reference_trajectories';
addpath(path_to_reference_trajectories);
path_to_yalmip_controllers = './yalmip_controllers';
addpath(path_to_yalmip_controllers);

% Since the individual simulations will save their own videos and data, this script
% passes in the paths to the data and videos folders to the mpc function.
% Call this struct 'paths'.

% The path to our videos repository now needs to be set, since we're no longer pushing videos to this simulations repository.
% Drew (on 2016-04-19) put the ultra-spine-videos repository in the same folder as ultra-spine-simulations, which means
% the path to the MPC videos folder will be:
paths.path_to_videos_folder = '../../../ultra-spine-videos/simulations/control/mpc/';

% The path to the folder where we'll store the .mat data from this simulation also needs to be set:
paths.path_to_data_folder = '../../data/mpc_data/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define default optimization parameters
% One large struct with multiple parameters. This includes parameters for the reference trajectory, only used if parameterizable.

% dt = time step for dynamics
% optimization_weights = a struct, various weights for the MPC objective function
% flags = a struct, various flags for the optimization procedure
% spine_geometric_parameters =  a struct, the size and shape of the spine. Loaded in from the dynamics generation script.
% links = number of moving vertebrae (tetrahedrons.) Equals 1 - N_tetras from the spine geometry struct.
%   NOTE that this script is not yet generalized past links == 3 (four tetra total.)
% tetra_vertical_spacing = Tetrahedron vertical spacing. The initial z-distance between successive tetrahedra. in meters.
% frame = animation frame divisor
% rest_lengths = the rest lengths of the cables, for the dynamics simulation, for a single vertebra. 
%   The first four rest lengths are for vertical cables, the second four are for saddle.
% num_points_ref_traj_tracking = the length of the first part of the reference trajectory, during the trajectory tracking part (e.g. movement) 
%   In terms of number of timesteps. ONLY USED FOR CERTAIN TRAJs
% num_points_ref_traj_regulation = the length of the second part of the reference trajectory, during the regulation part (e.g. stabilizing 
%   after movement.) In terms of number of timesteps. Augmenting the trajectory with a regulation portion only happens after the initial
%   reference trajectory is returned. See add_regulation_to_traj for more details.
%   NOTE that the total trajectory length is num_points_ref_traj = num_points_ref_traj_tracking + num_points_ref_traj_regulation.
%   This calculation is only performed in mpc_single_simulation - this outer script does not have the combined num_points_ref_traj.
% direction = the bending direction of a reference trajectory, either clockwise or counterclockwise. ONLY USED FOR CERTAIN TRAJs
%   direction is 1 for cw, -1 for ccw.
% horizon_length = the length of the horizon for MPC. This is the number of steps forward that the controller optimizes over. Was = 10.
% opt_time_limit = time limit enforced on the optimization solver. See the controller generation function's sdpsettings.

% Optimization weights:
% NOTE that this is ONLY used for a few yalmip controllers.
% In particular, on 2016-04-23 controllers that use this struct are get_yalmip_controller_XZG.
% obj_w_ref_xyz = Power-function weight for the objectives, used on the reference-tracking terms, for the longitudinal coordinates x,y,z
% obj_w_ref_angle = Power-function weight for the objectives, used on the reference-tracking terms, for the angle 
% obj_w_smooth = Power-function weight for the objectives, used on the successive-states terms (smooth motion)
% obj_w_input_pow = Power-function weight for the objectives, used on the successive-input terms (control authority, how-strong-is-the-motor)
% obj_w_input_mult = Multiplicative weight for the objectives, used on the successive-input terms (control authority, how-strong-is-the-motor)
% weighting_ratio = ratio for re-scaling the weights of each rigid body.
%   The lowest rigid body (vertebra 1) is not re-weighted, and the others are re-weighted linearly up to highest = weighting_ratio.
%   This value is multiplicative on top of the weights already proscribed. So, weights are (for ex., top tetra) obj_w_ref_xyz * weighting_ratio.
% vertebrae_do_not_track = cell array of numerical values of the vertebrae to NOT track. For ex., {} tracks all, and {1,2} is top only.
%   See generate_Q_rigidbody for more information.
% stab_const = distance that defines the region of the stability constraint on the optimization. 
%   The squared error at the last point in the horizon is enforced to be less than this. TO-DO: feasibility?

% The weights for Abishek's original MPC run were, respectively, 25, 0, 3, 3, 1/24. These were for only the top tetrahedron.

% Other notes about parameters:
% 'links' must be consistent with the dynamics defined in duct_accel.m and associated files! Those dynamics are pre-calculated,
%   and this parameter does NOT change them - it only affects the controller. If the controller and dynamics are inconsistent,
%   bad things will happen.
%   In the future, links will be defined as links = N_tetras-1.

optimization_parameters.dt = 0.001;

optimization_weights.obj_w_ref_xyz = 25;
optimization_weights.obj_w_ref_angle = 20;
optimization_weights.obj_w_smooth = 3;
optimization_weights.obj_w_input_pow = 1; % working with = 1 to 3, at least. 1 means no power-based weighting
optimization_weights.obj_w_input_mult = 1; % working with = 1/24 to 1, at least.
optimization_weights.weighting_ratio = 1;
optimization_weights.vertebrae_do_not_track = {};
optimization_weights.stab_const = 0;    % UNUSED AS OF 2016-04-29

optimization_parameters.optimization_weights = optimization_weights;

load(spine_geometric_parameters_path);
optimization_parameters.spine_geometric_parameters = spine_geometric_parameters;

optimization_parameters.links = 3;
optimization_parameters.tetra_vertical_spacing = 0.1;
optimization_parameters.frame = 3;

restLengths(1:4) = 0.1;
restLengths(5:8) = 0.187;
optimization_parameters.restLengths = restLengths;

% As of 2016-05-02, this script no longer defines the full length of a traj, but instead the length of each of the two parts.
%optimization_parameters.num_points_ref_traj = 80; 
optimization_parameters.num_points_ref_traj_tracking = 80; % for invkin, 80 gives a 0.0015 straight-line dist b/w points.
optimization_parameters.num_points_ref_traj_regulation = 0;
optimization_parameters.direction = -1; % 1 for cw, -1 for ccw.
optimization_parameters.horizon_length = 10;
optimization_parameters.opt_time_limit = 8; % seconds

% Flags:
% noise = adds noise to the forward simulation of the dynamics. noise = 1 turns it on.
% save_video = saves a video file, see function for more details.
% save_data = saves a .mat file with the simulation results, see function for more details.
% stringEnable = set to 1 plots the cables of the robot in the visualization
% run_lqr = flag to run the finite-time receding-horizon LQR after MPC. 0 = MPC only, 1 = run lqr also
% traj_is_full_system = flag representing whether traj is size 12 * whatever or 36 * whatever.
%   NOTE that this full_system flag is calculated automatically below after loading in the reference traj.
%   It is then inserted into the flags struct, before passing in to ultra_spine_mpc_single_simulation.

flags.noise = 1;
flags.save_video = 1;
flags.save_data = 1;
flags.stringEnable = 1;
flags.run_lqr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define default plotting parameters

% figure_window_location = the pixels of the two corners of the figure window
% figure_window_color = the color of the figure window. Usually white, 'w'.
% time = vector of time points for the visualization. This is NOT for the mpc or dynamics simulation.
%   time is NOT the time for the mpc, it's only for the frames for the visualization/video. 
%   the dynamics simulation length is determined by the reference trajectory number of timesteps * dt.
% rad = radius of a "leg" of the tetrahedron, for visualization. 
%   NOTE that this is a point mass model, so this variable is only graphical, not related to the mpc.
% cmaps = the color map used for plotting. Mostly for the tetrahedra colors.
% figure_rotation = the azimuthal rotations of the figure. Changes viewpoint, for better visualization.
% fontsize = the size of the text used in the figure
% cable_color = the color used in plotting the robot's cables
% cable_thickness = the thickness of the lines used to plot the robot's cables
% trajectory_color = the color of the line drawn of the centers of the tetrahedra from the reference traj
% trajectory_thickness = the thickness of the line drawn of the centers of the tetrahedra from the ref traj
% mpc_result_color = the color of the resulting trajectory output from MPC. Compare this one to traj.
%   NOTE that this isn't really 'color', it's matlab's designation for plotting, so can include symbols like . and -
% mpc_result_thickness = thickness of the resulting trajectory output from MPC. Compare to traj.
% plot_dt = timestep for animation. Used to slow the simulation down.
% plotting_offset = frame offset for animation. See ultra_spine_mpc for its use.
% lqr_result_color = color of the lines used to plot the centers of the vertebrae after running LQR.
% lqr_result_thickness = thickness of the lines used to plot the centers of the vertebrae after running LQR.
% anchor = position of the anchor location of the cables on one of a vertebra's nodes
% video_quality = option that controls what type of video to output. 'low' is Motion JPEG 2000, 'high' is 

plotting_parameters.figure_window_location = [0, 0, 600 700];
plotting_parameters.figure_window_color = 'w';
plotting_parameters.time = 0:optimization_parameters.dt:500;
plotting_parameters.rad = 0.01;
plotting_parameters.cmaps = gray(512); % summer(512);
plotting_parameters.figure_rotation = [-20, 14];
plotting_parameters.fontsize = 24;
plotting_parameters.cable_color = 'r';
plotting_parameters.cable_thickness = 2;
plotting_parameters.trajectory_color = 'b';
plotting_parameters.trajectory_thickness = 2;
plotting_parameters.mpc_result_color = 'c-';
plotting_parameters.mpc_result_thickness = 2;
plotting_parameters.plot_dt = 0.01;
plotting_parameters.plotting_offset = 30;
plotting_parameters.lqr_result_color = 'g';
plotting_parameters.lqr_result_thickness = 2;
plotting_parameters.anchor = [0 0 plotting_parameters.rad];
plotting_parameters.video_quality = 'low';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create list of variations from default parameters that will be run

% Save the default parameters:
default_optimization_parameters = optimization_parameters;
default_flags = flags;
default_plotting_parameters = plotting_parameters;

% Number of iterations of MPC. 
% NOTE that this must be carried throughout all the list-making: if there aren't the correct number of
% parameter sets, reference trajectories, and controllers, errors will be thrown.

num_mpc_runs = 1;

optimization_parameters_by_iteration = cell(num_mpc_runs, 1);
flags_by_iteration = cell(num_mpc_runs, 1);
plotting_parameters_by_iteration = cell(num_mpc_runs, 1);

% Populate these cell arrays with the default sets of parameters
for i=1:num_mpc_runs
    optimization_parameters_by_iteration{i} = default_optimization_parameters;
    flags_by_iteration{i} = default_flags;
    plotting_parameters_by_iteration{i} = default_plotting_parameters;
end

% Manually change these according to the runs of MPC that this script should make.

% Run 1:
plotting_parameters_by_iteration{1}.video_quality = 'high';
%optimization_parameters_by_iteration{1}.vertebrae_do_not_track = {1,2};
% for the circletop traj:
%optimization_parameters_by_iteration{1}.num_points_ref_traj_tracking = 180;
%optimization_parameters_by_iteration{1}.optimization_weights.vertebrae_do_not_track = {1,2};
%optimization_parameters_by_iteration{1}.num_points_ref_traj_tracking = 240;

% Run 2:
%optimization_parameters_by_iteration{2}.num_points_ref_traj_tracking = 80;
%flags_by_iteration{2}.noise = 1;
%optimization_parameters_by_iteration{2}.optimization_weights.vertebrae_do_not_track = {1,2};
%optimization_parameters_by_iteration{2}.optimization_weights.obj_w_input_mult = 1/4;
%optimization_parameters_by_iteration{2}.optimization_weights.obj_w_ref_xyz = 100;

% Run 3+
%optimization_parameters_by_iteration{3}.num_points_ref_traj_tracking = 160;
%flags_by_iteration{3}.noise = 1;
%optimization_parameters_by_iteration{3}.optimization_weights.vertebrae_do_not_track = {1,2};
% For MPC with larger numbers of timesteps, put a higher penalty on inputs
% that will hopefully make the motion more stable...
% optimization_parameters_by_iteration{3}.num_points_ref_traj_tracking = 150;
% optimization_parameters_by_iteration{4}.num_points_ref_traj_tracking = 200;
% optimization_parameters_by_iteration{5}.num_points_ref_traj_tracking = 300;
% optimization_parameters_by_iteration{6}.num_points_ref_traj_tracking = 500;
% optimization_parameters_by_iteration{7}.num_points_ref_traj_tracking = 700;
% optimization_parameters_by_iteration{8}.num_points_ref_traj_tracking = 1000;
% optimization_parameters_by_iteration{9}.num_points_ref_traj_tracking = 1500;

% a bit of error checking to make sure this simulation runs correctly.
assert( num_mpc_runs == size(optimization_parameters_by_iteration, 1), ...
            'Error! num_mpc_runs not equal to size of optimization_parameters_by_iteration. MPC simulation will fail.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create list of trajectories that should be run

% SOME NOTES ON ROTATION DIRECTION: these were found by changing ONLY these variables. Will not combine independently (recall, Euler angles.)
% These are all from the perspective of a ray starting at zero looking down the positive axis. (I've denoted it X+, Y+, etc.)
% T+ clockwise around the X+ axis (rotates in the Y,Z plane)
% T- counterclockwise around the X+ axis (rotates in the Y,Z plane)
% G+ clockwise around Y+ axis (rotates in the X,Z plane)
% G- counterclockwise around Y+ axis (rotates in the X,Z plane)
% P+ clockwise around the Z+ axis (rotates in the X,Y plane)
% P- counterclockwise around the Z+ axis (rotates in the X,Y plane)

% TO-DO: better format for loading in the trajectories, error checking.
% NOTE THAT THIS IS CURRENTLY BROKEN: Since the different 'get_ref_traj' commands have different inputs, we can't
% generalize a call to all of them. So, as of 2016-04-24, only allow the default inv_kin_XZG to run.

% Trajectories for *top tetrahedron only*
%[traj, ~] = get_ref_traj_circletop();
%[traj, ~] = get_ref_traj_quartercircletop();
%[traj, ~] = get_ref_traj_topbending1(); % Has trajectories along angles. NOT WORKING WELL as of 2016-02-28...
%[traj, ~] = get_ref_traj_topbending2();
%[traj, ~] = get_ref_traj_topbending_YZ();
%[traj, ~] = get_ref_traj_toprotationtest(); % Has trajectories along angles.
%[traj, ~]  = get_ref_traj_zero();

% Trajectories for *all tetrahedra*
% These have more input variables, define them here.
%[traj, ~] = get_ref_traj_allbending_ccw_XZG(tetra_vertical_spacing)
% [traj, ~] = get_ref_traj_invkin_XZG(optimization_parameters.tetra_vertical_spacing, ...
%                 optimization_parameters.num_points_ref_traj_tracking, ...
%                 optimization_parameters.direction);

% Create a cell array of strings that represent the trajectories to run.
% These will be eval'd to get the names of the functions to call.
% Each function starts with 'get_ref_traj_', so exclude that part.
trajectories_list = cell(num_mpc_runs, 1);

% Name a default trajectory. The trajectory from the inverse kinematics script is a good default:
default_traj = 'invkin_XZG';
% The "partial" trajectory was my attempt to have a very very small amount of movement.
% It might look like nothing is happening with that trajectory, since the movement is so small.
%default_traj = 'invkin_XZG_partial';

% Populate the list with this default
for i=1:num_mpc_runs
    trajectories_list{i} = default_traj;
    %trajectories_list{i} = 'invkin_XZG_partial';
end

% Manually change these according to the runs of MPC that this script should make.

% Run 1:

% try the circletop trajectory. This should match Abishek's original results.
%trajectories_list{1} = 'circletop_allvertebrae';

% Upright equilibrium trajectory.
%trajectories_list{1} = 'upright';

% Partial inverse kinematics trajectory
%trajectories_list{1} = 'invkin_XZG_partial';

% Run 2:
%trajectories_list{i} = ;

% Run 3+

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create list of controllers to use


% AS OF 2016-04-24, THIS IS NOT IMPLEMENTED, since the XYZTGP controller is quite general already.
% TO-DO: find a way to have the XYZTGP controller have different weights per dimension, not just one for all XYZ etc.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iterate over MPC runs. Begin the primary loop of this script

% a cell array to store the results of each MPC run
mpc_results = cell(num_mpc_runs, 1);

% a bit of debugging
disp('There will be:');
disp(num_mpc_runs);
disp('iterations of MPC.');

for mpc_iteration = 1:num_mpc_runs
    
    disp('MPC iteration:');
    disp(mpc_iteration);
    
    %% Reference Trajectory: Load in one of the trajectories

    disp('Loading reference trajectory:');
    disp(trajectories_list{mpc_iteration});
    
    % In order to call the desired trajectory-returning function, this script
    % creates a string of the entire command to call, then evaluates that string.
    % So, create strings representing all the outputs and arguments to pass in to this function.
    % Recall that each of these arguments are also a function of the iteration of MPC.
    % TO-DO: can we implement MPC in python to make this less painful? Object-oriented programming would be wonderful here.
    outputs = '[traj, ~] = ';
    arguments = strcat('(', ...
        'optimization_parameters_by_iteration{mpc_iteration}.tetra_vertical_spacing', ...
        ',', ...
        'optimization_parameters_by_iteration{mpc_iteration}.num_points_ref_traj_tracking', ...
        ',', ...
        'optimization_parameters_by_iteration{mpc_iteration}.direction', ...
        ')');

    % Create the full string of the command to execute
    reference_loading_command = strcat( outputs, 'get_ref_traj_', trajectories_list{mpc_iteration}, arguments, ';');
    % Execute the command. This should load in the traj matrix.
    eval(reference_loading_command);
    
    % Trajectories for *all tetrahedra*
    % These have more input variables, define them here.
    %[traj, ~] = get_ref_traj_allbending_ccw_XZG(tetra_vertical_spacing)
    %[traj, ~] = get_ref_traj_invkin_XZG(optimization_parameters.tetra_vertical_spacing, ...
    %                optimization_parameters.num_points_ref_traj_tracking, ...
    %                optimization_parameters.direction);

    % Automatically check if the trajectory that was loaded is for the full spine (3 vertebrae) or just the top one.
    % Declare a flag variable:
    traj_is_full_system = 0;
    if ( size(traj,1) == 12)
        % Small system, top vertebra only
        traj_is_full_system = 0;
        disp('Reference trajectory is for the top vertebra only (12 states.)');
    elseif ( size(traj, 1) == 36)
        % large system, 3 vertebrae
        traj_is_full_system = 1;
        disp('Reference trajectory is for all three vertebra (36 states.)');
    else
        error('Script currently only configured for trajectories of size 12 and 36! Loaded trajectory is not.');
    end

    % Append this flag to the 'flags' struct that will be passed in to ultra_spine_mpc_single_simulation
    %flags.traj_is_full_system = traj_is_full_system;
    % The flags can be varied for each run, also:
    flags_by_iteration{mpc_iteration}.traj_is_full_system = traj_is_full_system;
    
    % Augment the trajectory with a series of copies of the final state. This creates a period of 'regulation' around the final state,
    % and is useful for empirical tests of controller convergence / 'overshoot' etc.
    [traj, ~] = add_regulation_to_traj(traj, optimization_parameters_by_iteration{mpc_iteration}.num_points_ref_traj_regulation);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    %% Controller Initialization
    disp('Initializing YALMIP Controller...')

    % TO-DO: have these functions create their own YALMIP variables. This copying the yalmip things around might be causing slow-downs.
    % TO-DO: remove the hard-coded '36' for number of states and replace with 12 * links.

    % Some notes:
    %   The number of cables per tetrahedron (8) and state variables per tetrahedron (12, rigid body states) are hard-coded here.
    %   Copy optimization_parameters.links here into a variable named just 'links' for ease of writing out these variable names.

    links = optimization_parameters_by_iteration{mpc_iteration}.links;

    % Initialize yalmip variables for changing controller parameters
    % Horizon length:
    %horizon_length = 10;
    horizon_length = optimization_parameters_by_iteration{mpc_iteration}.horizon_length;
    
    disp('This iteration has horizon length:');
    disp(horizon_length);
    
    % YALMIP variables:
    inputs = sdpvar(repmat(8*links, 1, horizon_length-1), repmat(1, 1, horizon_length-1));
    states = sdpvar(repmat(12*links, 1, horizon_length), repmat(1, 1, horizon_length));
    A_t = sdpvar(repmat(12*links, 1, 12*links), repmat(1, 1, 12*links));
    B_t = sdpvar(repmat(12*links, 1, 8*links), repmat(1, 1, 8*links));
    c_t = sdpvar(36, 1);

    prev_in = sdpvar(8*links, 1);

    % Two different behaviors here: reference trajectory variables will be smaller if the trajectory is only the top tetra.
    if( ~traj_is_full_system)
        reference = sdpvar(repmat(12, 1, horizon_length), repmat(1, 1, horizon_length));
    else
        % all 36 states (3 vertebra) are present in the trajectory to track
        reference = sdpvar(repmat(36, 1, horizon_length), repmat(1, 1, horizon_length));
    end

    % Create the YALMIP controller for computation of the actual MPC optimizations.
    % This function contains all the definitions of the constraints on the optimization, as well as the objective function.

    % *TO-DO* have the controllers throw an error if the reference passed in is the wrong size?

    % Unroll the weights here for ease of referring to them.
    optimization_weights = optimization_parameters_by_iteration{mpc_iteration}.optimization_weights;
    obj_w_ref_xyz = optimization_weights.obj_w_ref_xyz; 
    obj_w_ref_angle = optimization_weights.obj_w_ref_angle;
    obj_w_smooth = optimization_weights.obj_w_smooth;
    obj_w_input_pow = optimization_weights.obj_w_input_pow;
    obj_w_input_mult = optimization_weights.obj_w_input_mult;
    weighting_ratio = optimization_weights.weighting_ratio;
    vertebrae_do_not_track = optimization_weights.vertebrae_do_not_track;
    stab_const = optimization_weights.stab_const;
    % Unroll the time limit too
    opt_time_limit = optimization_parameters_by_iteration{mpc_iteration}.opt_time_limit;
    
    % A bit of debugging
    disp('This iteration will not track the following vertebrae:');
    if( size(vertebrae_do_not_track,1) == 0 )
        disp('(all vertebrae will be tracked.)');
    else
        disp(vertebrae_do_not_track);
    end

    % Create weighting matrices for get_yalmip_controller_XYZTGP
    % Even though these are only used for one controller out of the list below, there is no harm in making a few more matrices.
    % TO-DO: record these weighting matrices with each iteration of MPC, so it's easier to see exactly what happened.
    Q_track = generate_Q_rigidbody( [obj_w_ref_xyz; obj_w_ref_angle; 0; 0], weighting_ratio, vertebrae_do_not_track, links);
    Q_smooth = generate_Q_rigidbody( [obj_w_smooth; obj_w_smooth; 0; 0], 1, vertebrae_do_not_track, links);

    % Cell array of controllers are generated from 2-step to N-step horizon to
    % deal with situations at the end of the reference trajectory
    for k = 2:horizon_length

        disp('Generating YALMIP controller for horizon length:');
        disp(k);
        % Controllers for 12 states
        %[controller{k}, ~, ~, ~, ~] = get_yalmip_controller_XYZ(k, inputs, states, A_t, B_t, c_t, prev_in, reference);
        [controller{k}, ~, ~, ~, ~] = get_yalmip_controller_XYZT(k, inputs, states, A_t, B_t, c_t, prev_in, reference);
        %[controller{k}, ~, ~, ~, ~] = get_yalmip_controller_XYZG(k, inputs, states, A_t, B_t, c_t, prev_in, reference);
        %[controller{k}, ~, ~, ~, ~] = get_yalmip_controller_G(k, inputs, states, A_t, B_t, c_t, prev_in, reference);
        %[controller{k}, ~, ~, ~, ~] = get_yalmip_controller_T(k, inputs, states, A_t, B_t, c_t, prev_in, reference);
        %[controller{k}, ~, ~, ~, ~] = get_yalmip_controller_P(k, inputs, states, A_t, B_t, c_t, prev_in, reference);

        % Controllers for 36 states
        % Note that these take in the parameters for the weights of the optimization. See above definition.
        %[controller{k}, ~, ~, ~, ~] = get_yalmip_controller_XZG(k, inputs, states, A_t, B_t, c_t, prev_in, reference, optimization_weights);
        [controller{k}, ~, ~, ~, ~] = get_yalmip_controller_XYZTGP(k, inputs, states, A_t, B_t, c_t, prev_in, reference, ...
                                            Q_track, Q_smooth, obj_w_input_mult, obj_w_input_pow, stab_const, opt_time_limit);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    %% Run this iteration of ultra_spine_mpc_single_simulation.
    
    % Note: the weighting matrices Q are passed in to the single simulation only so that they can be saved in the .mat data file
    % for that specific run. They are NOT used in that routine.
    
    other_data_to_save = {};
    other_data_to_save{1} = {'Q_track', Q_track};
    other_data_to_save{2} = {'Q_smooth', Q_smooth};

    %mpc_results{mpc_iteration} = ultra_spine_mpc_single_simulation(traj, controller, optimization_parameters_by_iteration{mpc_iteration}, ...
    %                flags, plotting_parameters_by_iteration{mpc_iteration}, other_data_to_save, paths);
    
    % Pass in the flags per iteration also:
    mpc_results{mpc_iteration} = ultra_spine_mpc_single_simulation(traj, controller, optimization_parameters_by_iteration{mpc_iteration}, ...
                    flags_by_iteration{mpc_iteration}, plotting_parameters_by_iteration{mpc_iteration}, other_data_to_save, paths);
    
    % End of this iteration of MPC.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    %% Delete some variables so excess memory is not taken up between iterations.
    % In particular, clear all the YALMIP stuff.
    clear inputs states A_t B_t c_t prev_in reference controller
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%% Analyze the results.

script_end = clock;

% End of script.






