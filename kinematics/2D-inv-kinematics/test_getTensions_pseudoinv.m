% test_getTensions_pseudoinv.m
% A short script that loads some parameters and tries out the relaxed
% inverse kinematics formulation, as will be used for the 2D MPC.
% Andrew P. Sabelhaus 2018-03-29

% prepare the workspace
clear all;
close all;
clc;

% load in the parameters for the tetrahedral 2D spine
dynamics_path = '../../dynamics/2d-dynamics-symbolicsolver';
addpath(dynamics_path);
load('two_d_geometry.mat'); % loads a struct named "two_d_geometry"

% define an example system state. Make the vertebra 10 cm separated,
% as per the 2D MPC script.
% This is x, z, g, dx, dz, dg.
xi = [ 0; 0.1; 0; 0; 0; 0];
min_cable_tension = 5;

% Call the fcn and let's see what happens:
[tensions, restLengths, A, p, A_skelton_c, qOpt, qOpt_lax] = getTensions_pseudoinv(xi, two_d_geometry, min_cable_tension);