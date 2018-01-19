# ultra-spine-simulations\kinematics

Copyright 2015-2018 Berkeley Emergent Space Technologies Lab.

There is a variety of inverse kinematics code in these folders.
As of 2018-01-18, here's a description of the files, what they do, and if they're deprecated.

- inv-kinematics: Software used in the 2015 IDETC paper.
  This is *deprecated* and not maintained.
  It uses Jeff Friesen's approach, which is adapted for more versatile situations later.

- 2D-inv-kinematics: software developed in 2016-2017 from scratch, for a single-vertebra 2D spine.
  This uses equality constraints, and is hard-coded with the dimensions for one topology.
  It does not use Schek's connectivity matrix.
  Also *deprecated*

- general-solver: Most recent work on the inverse kinematics of a tensegrity spine.
  This is an adaptation of 2D-inv-kinematics to be more general, and is what's currently maintained.
  Files described in the README in that folder.


