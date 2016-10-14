This directory contains the files required to run model-predictive control 
for the Underactuated Lightweight Tensegrity Robotic Assistive Spine (ULTRA Spine.)

The structure of these files are the following.
Each sub-section means "this function is called by the one above it."

ultra_spine_mpc.m: the primary script. This sets up multiple simulations of MPC. 
                   It creates the reference trajectories and the YALMIP optimizer.
    get_ref_traj_inkvin_XZG.m: under reference_trajectories. Creates a sequence of points corresponding to a trajectory we want to track,
                               adapted from an earlier inverse kinematics simulation of the spine.
    get_yalmip_controller_XYZTGP.m: under yalmip_controllers. Creates all the optimizer objects for the YALMIP optimizations.
                                    This is where the optimization solver is specified, currently as "gurobi".
    ultra_spine_mpc_single_simulation.m: This function runs one of the individual simulations of MPC. 
                                         It mostly does the plotting/visualization.
	plotSpineLink.m: creates the visualization of a spine vertebra (also called link in various places in these scripts),
			 based on the "Tetra" matrix of the vectors to each of the point masses of a vertebra.
	getHG_Tform.m: creates the rotation matrix for a specific vertebra based on the system states.
        simulate_mpc_traj.m: this function actually executes the YALMIP controller. It applies a YALMIP "optimizer" object, as output
                             by get_yalmip_controller_XYZTGP, which are stored in the controllers{k} struct, to a set of system states.
                             It also calls the dynamics function to forward-simulate the states after each step of the controller is called.
            linearize_dynamics.m: this function creates the specific A,B, and c matrices for a particular timestep in the simulation.
                                  These constant matrices are then plugged-in to the YALMIP optimizer. These are for the At, Bt, and ct in
                                  the "controller formulation" section of the ACC 2017 paper, section IV.
            simulate_dynamics.m: this function forward-simulates the system for one timestep. It includes hard-coded constants for
                                 the added disturbance. It uses the Runge-Kutta method for forward simulation.
                                 This is the function that calls the underlying system dynamics in 
                                 the ../../dynamics/3d-dynamics-symbolicsolver folder, which are:
                getTensions.m: function that calculates the tensions in the cables, given the positions of the vertebrae
                duct_accel.m: calculates x_dot for each vertebra, given the forces acting on it (the tensions in the cables, plus gravity)
                              Note that this function is named after one of Jeffrey Friesen's older robots, the "DuCTT", 
                              which was originally simulated using these scripts. We haven't renamed the file yet.

As of this commit, the ultra_spine_mpc.m script should create the results from the ACC 2017 paper.
However, users will have to 
(1) change the optimization solver, if gurobi is not installed, and
(2) change the path to the videos folder to wherever else you'd like to store the video results from the simulations.
(3) make sure YALMIP is installed properly and added to the MATLAB path.