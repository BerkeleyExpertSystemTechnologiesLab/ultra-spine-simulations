This folder contains code that generates a set of dynamics equations for the ULTRA Spine robot.
It is derived from Jeffrey Friesen's work on the DuCTT robot.

The dynamics are present in function duct_accel.m, and the functions lengths.m and dlengths_dt.m provide support.
The script spineDynamics.m generates these files.

Note that spineDynamics.m requires the 'fulldiff.m' script.
Credit goes to Tim Jorris for that function.
