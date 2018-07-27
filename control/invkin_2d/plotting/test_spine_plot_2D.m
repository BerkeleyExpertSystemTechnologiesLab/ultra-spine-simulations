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
load('two_d_geometry.mat');

% From the 3D version: we set a bunch of plotting options.
% Copy these here, and see what's actually used. Coloring?
figure_window_location = [0, 0, 600 700];
figure_window_color = 'w';
cmaps = gray(512); % summer(512);
figure_rotation = [-20, 14];
fontsize = 24;
cable_color = 'r';
cable_thickness = 2;
trajectory_color = 'b';
trajectory_thickness = 2;
mpc_result_color = 'c-';
mpc_result_thickness = 2;
plotting_offset = 30;
lqr_result_color = 'g';
lqr_result_thickness = 2;
video_quality = 'low';

% Here's their setup:
% Create the figure window
%figure_handle = figure('position', figure_window_location,'Color',figure_window_color);

% TO-DO: figure out what this M variable is and how it's used.
%M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

% Configure the figure: colors, orientation, etc.
colormap(cmaps(1:256,:));
%ax = axes();
grid on;
axis equal;
hold on;
%view(figure_rotation);

% Labels and text
%xlabel('X (m)')
%ylabel('Y (m)')
%zlabel('Z (m)')
%title('ULTRA Spine Model')
 
% Size everything properly
%xlim([-0.2 0.2])
%ylim([-0.2 0.2])
%zlim([-0.1, 0.4])
%set(gca,'FontSize',fontsize)

shading interp
light
lighting phong

% Set up the window.
%figure;
%hold on;
%axis([-0.2, 0.2, -0.1, 0.3]);

% for surf:
%ax = axes();
% so we don't create a new axis, just get the old one:
ax = gca;

% pick out a simple system state: vertebra at some vertical displacement.
xi = [0; 0.1; 0; 0; 0; 0];
% to-do: seems not plotting the right vertebra???

handles = plot_2d_tensegrity_surfaced(xi, two_d_geometry, ax);
%handles = plot_2d_tensegrity(xi, two_d_geometry);
drawnow;