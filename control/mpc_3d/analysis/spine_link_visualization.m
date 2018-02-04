% spine_link_visualization.m
% Abishek Akella, Andrew P. Sabelhaus
% This is a file for plotting the spine link with corresponding labels

%% Initialize system parameters
%clear variables;
clear all;
close all;
clc;


path_to_dynamics = '../../../dynamics/3d-dynamics-symbolicsolver';
addpath(path_to_dynamics);
% Parameters for plotting:
% NOTE that the rod sizes here are only used for plotting, 
% since the dynamics used in this work is a point-mass model.

% parent folder, needs plotSpineLink
addpath('..');

% Load parameters in from the saved file that accompanies the dynamics
spine_geometric_parameters_path = strcat(path_to_dynamics, '/spine_geometric_parameters.mat');
load(spine_geometric_parameters_path);
% Unroll this struct into individual variables. See the dynamics generation script for more information.
% Gravitational force
g = spine_geometric_parameters.g;
% Total number of spine tetrahedrons
N_tetras = spine_geometric_parameters.N;
% Length of one "leg" of the tetrahedron (straight-line distance from center to outer point)
l = spine_geometric_parameters.l;
% Height of one tetrahedtron
h = spine_geometric_parameters.h;
% total mass of one whole tetrahedron
m_t = spine_geometric_parameters.m_t;
% Factor-of-safety with respect to tetrahedron mass (note: this is unused in this script)
FoS = spine_geometric_parameters.FoS;
% mass of one node of the tetrahedron (there are five point masses per tetra)
%m = spine_geometric_parameters.m;

% Radius of a "leg" of the tetrahedron. NOTE that since this is a point-mass model, this parameter is only for visualization.
rad = 0.01;


% Tetrahedron vertical spacing. The initial z-distance between successive tetrahedra
tetra_vertical_spacing = 0.1; % meters

%% Initialize Plot
% Note that this is performed at the beginning so the visualization of the terahedra bodies can be loaded properly.

% Create the figure window
%figure_handle = figure('position', [100, 100, 700 800],'Color','w');
figure_handle = figure('position', [50, 50, 600, 600],'Color','w');

time = 0:0.01:500;
M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

% Set the color map
cmaps = summer(512);
colormap(cmaps(1:256,:))
%colormap(cool);
ax = axes();

grid on;
axis equal;

hold on;

% Rotate for a better visualization
view([-20, 14]);

% Labels

% Size everything properly
xlim([-0.3 0.3])
ylim([-0.3 0.3])
zlim([-0.1, 0.3])
set(gca,'FontSize',24)

shading interp
light
lighting phong

%m = 256;

%% Initialize the simulation
restLengths(1) = 0.1; % Vertical cable rest length
restLengths(2) = 0.1;
restLengths(3) = 0.1;
restLengths(4) = 0.1;
restLengths(5) = 0.187; % Saddle cable rest length
restLengths(6) = 0.187;
restLengths(7) = 0.187;
restLengths(8) = 0.187;

% The initial state for all of these 36 variables.
x_initial = [];

% Initialize all the states of the system
% This script currently (2016-02-27) considers "k" to be a different index in different circumstances.
% Here, it's used for the system states as k=1 for the first moving tetrahedron, and k=3 for the topmost one.
% But, for plotting, each of these is shifted up: Tetra{1} and transform{1} are for the first (unmoving) tetra, and Tetra{4} is the topmost.

% Plot the first tetrahedron body (k=1 in the Tetra{} usage)
Tetra{1} = [(l^2 - (h/2)^2)^.5, 0, -h/2; ...
            -(l^2 - (h/2)^2)^.5, 0, -h/2; ...
            0, (l^2 - (h/2)^2)^.5, h/2; ...
            0, -(l^2 - (h/2)^2)^.5, h/2];

% Plot a visualization of this spine tetrahedron
[transform{1}, ~] = plotSpineLink(Tetra{1}, rad, ax);
plotCurvedArrows();

%% Graphics Options
ax1 = gca;
yruler = ax1.YRuler;
yruler.Axle.Visible = 'off';
xruler = ax1.XRuler;
xruler.Axle.Visible = 'off';
zruler = ax1.ZRuler;
zruler.Visible = 'off';
set(gca, 'YTick', []);
set(gca, 'XTick', []);
set(gca, 'ZTick', []);

view([1, 1, 0.75]);
zlabel = text(0, 0.05, 0.3, 'Z'); zlabel.FontWeight = 'bold';
xlabel = text(0.3, 0, 0.05, 'X'); xlabel.FontWeight = 'bold';
ylabel = text(0, 0.3, 0.05, 'Y'); ylabel.FontWeight = 'bold';

plabel = text(0, -0.09, 0.18, '\phi'); plabel.FontWeight = 'bold'; plabel.FontSize = 12;
tlabel = text(0.2, 0, 0.065, '\tau'); tlabel.FontWeight = 'bold'; tlabel.FontSize = 12;
glabel = text(0, 0.2, 0.065, '\gamma'); glabel.FontWeight = 'bold'; glabel.FontSize = 12;

annotation(gcf,'arrow',[0.430625 0.430625],...
    [0.364261813537677 0.349936143039593],'Color',[0 0 1]);
annotation(gcf,'arrow',[0.604375 0.604375],...
    [0.362984674329503 0.348659003831419],'Color',[0 0 1]);
annotation(gcf,'arrow',[0.532500000000001 0.525625],...
    [0.619689655172416 0.628352490421457],'Color',[0 0 1]);