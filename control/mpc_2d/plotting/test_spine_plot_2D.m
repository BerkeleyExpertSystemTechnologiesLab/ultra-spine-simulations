% test_spine_plot_2D.m
% script that tests out the updated plot_2d_tensegrity_surfaced function.
% Andrew P. Sabelhaus, Feb 2018

clear all;
close all;
clc;

% call the plot_2d_tensegrity_surface function.
% Need to load a geometry struct first...
% System dynamics path, spine geometric parameters file is also in this
% directory
dynamics_path = '../../../dynamics/2d-dynamics-symbolicsolver';
addpath(dynamics_path)
% creates the two_d_geometry struct.
load('two_d_geometry.mat')

% Set up the window.
figure;
hold on;
axis([-0.2, 0.2, -0.1, 0.3]);

% for surf:
ax = axes();

% pick out a simple system state: vertebra at some vertical displacement.
xi = [0; 0.1; 0; 0; 0; 0];
% to-do: seems not plotting the right vertebra???

%handles = plot_2d_tensegrity_surfaced(xi, two_d_geometry, ax);
handles = plot_2d_tensegrity(xi, two_d_geometry);
drawnow;