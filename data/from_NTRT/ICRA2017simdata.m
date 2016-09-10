% ICRA2017simdata.m
% A helper script file that takes in data from 
% the NASA Tensegrity Robotics Toolkit simulations of the ULTRA Spine
% and its force plates.
%   Andrew P. Sabelhaus
%   Berkeley Emergent Space Tensegrities Lab

% Clean up the workspace
clear all;
close all;
clc;

% Add the path to the hline and vline functions.
% @TODO make this more robust!
addpath('./hline_vline');

% The log file base path:
% (hard-coded to Drew's computer)
logfile_base = '~/repositories/NTRTsim/resources/src/forcePlate/forcePlateDemo/logs/';
% The timestamp for the file to read in
% Copied from the name of the log file itself
logfile_timestamp = '09102016_113901';
% The calibration factor for the force plate readings.
% Currently 5, since gravity = 98.1 and scale = 0.5.
calib_factor = 5;
% A flag to control making plots or not
make_plots = 1;

% Call the parser function
fpdata = parseNTRTForcePlateData(logfile_base, logfile_timestamp, calib_factor, make_plots);
