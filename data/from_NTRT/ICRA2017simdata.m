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
%logfile_base = '~/repositories/NTRTsim/resources/src/forcePlate/forcePlateDemo/logs/';
logfile_base = '~/repositories/NTRTsim/resources/src/forcePlate/AppHorizontalSpine/logs/';
%logfile_base = '~/repositories/NTRTsim/resources/src/forcePlate/AppRotatingVertebraSpine/logs/';
% The timestamp for the file to read in
% Copied from the name of the log file itself
% Top Left:
%logfile_timestamp = '09102016_113901';
% Top Right:
%logfile_timestamp = '09102016_150150';
% Corrected top right, with symmetric spine:
%logfile_timestamp = '09122016_214927';
% Corrected top left, with symmetric spine:
%logfile_timestamp = '09122016_215445';
% Corrected top right, with symmetric spine, slightly larger:
%logfile_timestamp = '09122016_221708'
% Corrected top left, with symmetric spine, slightly larger:
%logfile_timestamp = '09122016_221235'
% Testing:
% Top right:
%logfile_timestamp = '09132016_184234';
% Top left:
%logfile_timestamp = '09132016_103617';
% The calibration factor for the force plate readings.
% Currently 5, since gravity = 98.1 and scale = 0.5.
%calib_factor = 5;
%calib_factor = 1;
% A flag to control making plots or not
make_plots = 1;

% To get the calibration factor:
% Manually add up the total force from all the robot feet
% at the start of the simulation.
% This is fpdata{i}.y_foranalaysis(1).
%robot_totalF_NTRT = 35.79 + 38.6715 + 67.33 + 55.81 % for rotating vertebra
robot_totalF_NTRT = 72.3028 % for bending
% The total weight of the robot in hardware is (kg from kitchen scales * g)
robot_totalF_hardware = 1.365 * 9.81
% Since we're diving by the calibration factor, the converstion
% is then (NTRT weight / hardware weight), or
calib_factor = robot_totalF_NTRT / robot_totalF_hardware;
%calib_factor = 1;

% For the rotating vertebra:
logfile_timestamp = '02282017_194142';

% For single-bending:
logfile_timestamp = '03012017_140509';

% Call the parser function
fpdata = parseNTRTForcePlateData(logfile_base, logfile_timestamp, calib_factor, make_plots);
