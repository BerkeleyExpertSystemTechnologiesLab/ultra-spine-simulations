# ultra-spine
All software for the Underactuated Lightweight Tensegrity Robotic Assistive Spine project. This repository contains both published code (public) and unpublished code (private), and all code for the project should be here in a folder, possibly as a submodule.

Copyright 2015 Berkeley Emergent Space Technologies Lab.

Previous ULTRA Spine software and code has been organized into the following folders:

1) 2DSpineSim (Zeerek Ahmad) -> simulations/dynamics/ultra-spine-2D-analytical-dynamics
2) ULTRA-Spine (Abishek Akella) -> multiple locations.
	    *Dynamics work, with Jeff Friesen: simulations/dynamics/ultra-spine-symbolicsolver-dynamics
	    *Model Predictive Controls work -> simulations/control/ultra-spine-mpc
	    Sequential Quadratic Programming trajectory generation -> simulations/control/ultra-spine-sqp-traj-gen
3) ultraSpineDynamics (ChanWoo Yang) -> simulations/dynamics/ultra-spine-3D-analytical-dynamics
4) *TensegritySpineInverseKinematics (Drew Sabelhaus, Jeff Friesen) -> simulations/kinematics/invkinematics-flemons
5) Julia code for EE249A project (Abishek Akella, Aldrich Ong) -> hardware/1segment-servo

*This work is a submodule.
Other work will move to a public repository and reconfigured as a submodule once it's time to publish it.
