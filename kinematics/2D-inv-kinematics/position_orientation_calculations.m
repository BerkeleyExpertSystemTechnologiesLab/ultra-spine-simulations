% Calculation of constants for inverse kinematics position and orientation tracking
% Andrew P. Sabelhaus 2016-04-22
% This is just a copy of a bunch of commands I ran from the MATLAB prompt. They may not work.

% Run SpineInverseKinematics first.

% According to the notes I took, c1 here is the constant for the theta-timestep curve,
% and c2 is the constant for the gamma-timestep curve. Theta is the angle drawn out by
% the arc following the successive positions of a tetrahedron center, with respect to the origin
% (think: two lines drawn from origin to tetra center, angle between these two lines.)
% Gamma is the rotation of a single tetrahedron around its own axis.
% This script is NOT equiped to run anything other than the version of SpineInverseKinematics
% with hard-coded constants as of 2016-04-22.

x3 = squeeze(centersHistory(3,2,:));
z3 = squeeze(centersHistory(3,3,:));
x4 = squeeze(centersHistory(4,2,:));
z4 = squeeze(centersHistory(4,3,:));
x5 = squeeze(centersHistory(5,2,:));
z5 = squeeze(centersHistory(5,3,:));
x2 = squeeze(centersHistory(2,2,:));
z2 = squeeze(centersHistory(2,3,:));
theta2 = atan(x2./z2);
theta4 = atan(x4./z4);
theta3 = atan(x3./z3);
theta5 = atan(x5./z5);
c2_2_history = zeros(19,1);
c2_3_history = zeros(19,1);
c2_4_history = zeros(19,1);
c2_5_history = zeros(19,1);
c1_2_history = zeros(19,1);
c1_3_history = zeros(19,1);
c1_4_history = zeros(19,1);
c1_5_history = zeros(19,1);
for i=2:20
c1_2_history(i-1) = theta2(i)-theta2(i-1);
c1_3_history(i-1) = theta3(i)-theta3(i-1);
c1_4_history(i-1) = theta4(i)-theta4(i-1);
c1_5_history(i-1) = theta5(i)-theta5(i-1);
end
for i=2:20
c2_2_history(i-1) = rotationHistory(2,i)-rotationHistory(2,i-1);
c2_3_history(i-1) = rotationHistory(3,i)-rotationHistory(3,i-1);
c2_4_history(i-1) = rotationHistory(4,i)-rotationHistory(4,i-1);
c2_5_history(i-1) = rotationHistory(5,i)-rotationHistory(5,i-1);
end
c1_2 = c1_2_history(1);
c1_3 = c1_3_history(1);
c1_4 = c1_4_history(1);
c1_5 = c1_5_history(1);
disp(c2_2_history);
disp(c2_3_history);
disp(c2_4_history);
disp(c2_5_history);
c2_2 = -0.01;
c2_3 = -0.02;
c2_4 = -0.03;
c2_5 = -0.04;
% Then, the curve relating theta and gamma is a straight line 
% with slope c2/c1.
c_rot_2 = c2_2 / c1_2;
c_rot_3 = c2_3 / c1_3;
c_rot_4 = c2_4 / c1_4;
x_rot_5 = c2_5 / c1_5;

% calculating the angle changes between tetras, for displacement.
% This is just searching for whatever pattern might be present/
% I noticed that there is a constant offset at timestep 20 (the end of the trajectory)
% from the inverse kinematics: e.g., theta(20) for tetra N is 0.005 * (20) * N + theta_offset, 
% where theta_offset is the SAME for all tetras at time=20. It's a weird value though, something like pi/35.
theta_offset2 = zeros(20,1);
theta_offset3 = zeros(20,1);
theta_offset4 = zeros(20,1);
theta_offset5 = zeros(20,1);
for i=1:20
    theta_offset2(i) = theta2(i) - 0.005 * (1) * i;
    theta_offset3(i) = theta3(i) - 0.005 * (2) * i;
    theta_offset4(i) = theta4(i) - 0.005 * (3) * i;
    theta_offset5(i) = theta5(i) - 0.005 * (4) * i;
end

% let's see if these offsets are nice, clean angles:
theta_offset2_coeff = pi ./theta_offset2;
theta_offset3_coeff = pi ./theta_offset3;
theta_offset4_coeff = pi ./theta_offset4;
theta_offset5_coeff = pi ./theta_offset5;
