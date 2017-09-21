This folder contains code that generates a set of dynamics equations for the ULTRA Spine robot.
It is derived from Jeffrey Friesen's work on the DuCTT robot.

The dynamics are present in function spine_accel.m, and the functions lengths.m and dlengths_dt.m provide support.
The script spineDynamics_underactuated.m generates these files.

Note that spineDynamics_underactuated.m requires the 'fulldiff.m' script.
Credit goes to Tim Jorris for that function.

To use this code:
run the spineDynamics_underactuated script to generate (or update) the three functions
Test out the dynamics. You might need to adapt some code from the MPC directory:
control/mpc_3d/simulate_mpc_traj
or
control/mpc_3d/simulate_spine_dynamics (...?)
These should show how to use the three functions in this folder to create a visualization of the spine
moving around. It might be a good idea to make a script or function to test out 
a bunch of initial conditions, which lets us see the changes in the spine behavior when the 
dynamics symbolic solver changes.