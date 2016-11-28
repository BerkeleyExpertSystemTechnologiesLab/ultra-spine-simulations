%% two_d_spine_simulation_example.m
% Andrew P. Sabelhaus and the Berkeley Emergent Space Tensegrities Lab 
% Copyright 2016

% This script simulates the two-dimensional "inverted-Y" tensegrity spine,
% as an example for others to use.

%% 1) Prepare the workspace

clear all;
close all;
clc;

%% 2) Declare constants, read in the geometric parameters of the spine

% Simulation timestep
dt = 0.001;
% Geometry of the spine, for plotting purposes
spine_geometric_parameters_path = 'spine_geometric_parameters_2D.mat';
load(spine_geometric_parameters_path);

%% 3) Set up the simulation

% The initial conditions of the spine will be: not moving,
% slightly upwards in z, zero in x, with some rotation.
% These are just guesses for now.
xi_0 = [0; 0.1; 0.1; 0; 0; 0];
% The inputs to the system will be constant here:
% let's not change the rest lengths, just have them all be tight for the moment.
% The first two rest lengths are for the vertical cables, the second two are for the saddle.
u = [0.01; 0.02; 0.01; 0.01];

% We'll simulate for the following amount of time, in seconds:
t = 1;
% That means, with our dt, we'll have the following iterations of 
% forward-Euler simulation:
steps = t / dt;

% Let's make matrices to store the system states over time
xi = zeros(6, steps+1);
% Put in xi_0 into this matrix
xi(:,1) = xi_0;

%% 4) Perform the simulation

for i=1:steps
    % Forward simulate this step
    % Note that we're not using any of the multiple-step functionality that's 
    % built-in to simulate_2d_spine_dynamics, e.g., only one forward-Euler integration
    % per call to simulate_2_spine_dynamics.
    xi(:,i+1) = simulate_2d_spine_dynamics(xi(:,i), u, dt, 1);
end

%% 5) Plot the results

% For ease of indexing, make a big matrix of dt times.
dt_list = 0:dt:t;

% Plot each system state
figure;
subplot(6, 1, 1);
plot(dt_list, xi(1,:));
ylabel('x');
subplot(6, 1, 2);
plot(dt_list, xi(2,:));
ylabel('z');
subplot(6, 1, 3);
plot(dt_list, xi(3,:));
ylabel('theta');
subplot(6, 1, 4);
plot(dt_list, xi(4,:));
ylabel('dx');
subplot(6, 1, 5);
plot(dt_list, xi(5,:));
ylabel('dz');
subplot(6, 1, 6);
plot(dt_list, xi(6,:));
ylabel('dtheta');

% end.


