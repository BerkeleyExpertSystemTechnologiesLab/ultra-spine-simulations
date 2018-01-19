# ultra-spine-simulations\kinematics

Copyright 2015-2018 Berkeley Emergent Space Technologies Lab.

The general-solver folder is the most up-to-date version of our inverse kinematics code.
Here is a description of the files/folders as of 2017-01-18:

- Brian Comparison: data that compares the 'equality-constraint' method with Brian Cera's results. Brian's code not present here, only results. Deprecated, unused.

- Old Solver versions: deprecated.

NOW, for the code that's used regularly:

- BodyForceReader: a function, used to check Inv Kin results versus static equilibrium.
  Use this function for checking results. Should show no accelerations if Inv Kin is applied correctly.

- generalInverseKinematicsSolverv2 and v3: scripts.
  First attempts at 3D solver. Hard-coded, but uses a connectivity matrix.
  Equality-constraint method.
  To-do: check using BodyForceReader.
  Difference between v2 and v3: v2 does not do moments, but v3 does include moments.

- InvKin.m: function, first attempt at turning general solver into a function instead of a script.
  May not work correctly, or at least might be limited by equality-constraint method.

- InvKin_Horizontal_Spine.m: script, attempt to adapt general solver to a 3D spine that is horizontal instead of vertical.
  Produces errors, since the constraints are not treated properly right now.

- InvKin_pseudoinv.m: function, CURRENT WORK.
  Adapting InvKin.m to use inequality constraints, via Jeff Friesens' 2014 paper.
  Hopefully, this function will allow for feasible solutions to be found for a wider variety of configurations/positions of the spine.

- InvKin_Test_Run / 6bar: scripts that call the InvKin function.
  To-do: are these set up correctly??

- InvKin_Test_Vertical_Spine: script, CURRENT WORK.
  This script has Schek's connectivity matrix for a 3D vertical spine.
  It calls InvKin properly (we will be testing with InvKin_pseudoinv.)
  At the moment (2018-01-18), does not find feasible solution, which is not correct.

- simpleBarTestCase: script, uses a very simple tensegrity for InvKin.
  The script finds a feasible solution here, seemingly.


TO-DO: move these parts to different branches or folders...