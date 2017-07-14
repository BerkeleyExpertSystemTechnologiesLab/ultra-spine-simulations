% mpc_error_analysis_2D
% Copyright 2017 Andrew P. Sabelhaus and Berkeley Emergent Space Tensegrities Lab
% This script does the error analysis of the 2D Model-Predictive Controller for the spine.
% As of 2017-04-16, just makes the plots of Shirley's data for Drew's qualifying exam.

% This script ASSUMES that the following variables are in the workspace:
% xi_cl, xi_traj, u_cl, u_traj, opt_params

% Modeling things after the mpc_error_analysis 3D script,

fontsize = 12;
% Calculate the time in seconds at each point.
% We'll use the following variables a few times here, so pick it out for cleanliness:
num_pts = opt_params.num_pts;
dt = opt_params.dt;
t = dt : dt : num_pts*dt;
% Adjust time so that everything is in milliseconds. This makes things clearer.
t = t * 1000;

%% States plot:

% A good figure window setup is 'Position',[100,100,500,300].

% Plot the XZG errors on one plot.
xzg_handle = figure;
hold on;
set(gca,'FontSize',fontsize);
set(xzg_handle,'Position',[100,100,550,400]);

% For the first subplot (x), ...
subplot_handle = subplot(3, 1, 1);
hold on;
set(gca,'FontSize',fontsize);
% Plot x, adjusting the distances to cm
% Reference is
plot(t, xi_traj(1, 1:num_pts)*100, 'b.-', 'LineWidth', 1.5 , 'markersize', 15);
% Closed-loop behavior is
plot(t, xi_cl(1, 1:num_pts)*100, 'g.-', 'LineWidth', 1.5 , 'markersize', 10);
% Labels
title('State Tracking with Inv. Kin. Inputs');
ylabel('X (cm)');
%xlabel('time (msec)');
% Make the legend
legend('Reference', 'Result', 'Location', 'Southwest');

% For the second subplot (z), ...
subplot_handle = subplot(3, 1, 2);
hold on;
set(gca,'FontSize',fontsize);
% Plot z, adjusting the distances to cm
% Reference is
plot(t, xi_traj(2, 1:num_pts)*100, 'b.-', 'LineWidth', 1.5 , 'markersize', 15);
% Closed-loop behavior is
plot(t, xi_cl(2, 1:num_pts)*100, 'g.-', 'LineWidth', 1.5 , 'markersize', 10);
% Labels
ylabel('Z (cm)');
%xlabel('time (msec)');

% For the third subplot (theta), ...
subplot_handle = subplot(3, 1, 3);
hold on;
set(gca,'FontSize',fontsize);
% Plot theta, adjusting radians to degrees
% Reference is
plot(t, xi_traj(3, 1:num_pts)*180/pi, 'b.-', 'LineWidth', 1.5 , 'markersize', 15);
% Closed-loop behavior is
plot(t, xi_cl(3, 1:num_pts)*180/pi, 'g.-', 'LineWidth', 1.5 , 'markersize', 10);
% Labels
ylabel('\theta (deg)');
% Put the time label below the bottom graph.
xlabel('time (msec)');

%% Inputs plot.

% Plot the XZG errors on one plot.
u_handle = figure;
hold on;
set(gca,'FontSize',fontsize);
set(u_handle,'Position',[100,100,500,500]);

% For the first subplot (u1), ...
subplot_handle = subplot(4, 1, 1);
hold on;
set(gca,'FontSize',fontsize);
% Plot input 1 length, adjusting the distances to cm
% Reference is
plot(t, u_traj(1, 1:num_pts)*100, 'b.-', 'LineWidth', 1.5 , 'markersize', 15);
% Closed-loop behavior is
plot(t, u_cl(1, 1:num_pts)*100, 'g.-', 'LineWidth', 1.5 , 'markersize', 10);
% Labels
title('Input Tracking with Inv. Kin. Inputs');
ylabel('U1 (cm)');
%xlabel('time (msec)');
% Make the legend
legend('Reference', 'Result', 'Location', 'Best');

% For the second subplot (u2), ...
subplot_handle = subplot(4, 1, 2);
hold on;
set(gca,'FontSize',fontsize);
% Plot u2, adjusting the distances to cm
% Reference is
plot(t, u_traj(2, 1:num_pts)*100, 'b.-', 'LineWidth', 1.5 , 'markersize', 15);
% Closed-loop behavior is
plot(t, u_cl(2, 1:num_pts)*100, 'g.-', 'LineWidth', 1.5 , 'markersize', 10);
% Labels
ylabel('U2 (cm)');
%xlabel('time (msec)');

% For the third subplot (u3), ...
subplot_handle = subplot(4, 1, 3);
hold on;
set(gca,'FontSize',fontsize);
% Plot u3, adjusting the distances to cm
% Reference is
plot(t, u_traj(3, 1:num_pts)*100, 'b.-', 'LineWidth', 1.5 , 'markersize', 15);
% Closed-loop behavior is
plot(t, u_cl(3, 1:num_pts)*100, 'g.-', 'LineWidth', 1.5 , 'markersize', 10);
% Labels
ylabel('U3 (cm)');
%xlabel('time (msec)');

% For the fourth subplot (u4), ...
subplot_handle = subplot(4, 1, 4);
hold on;
set(gca,'FontSize',fontsize);
% Plot u4
% Reference is
plot(t, u_traj(4, 1:num_pts)*100, 'b.-', 'LineWidth', 1.5 , 'markersize', 15);
% Closed-loop behavior is
plot(t, u_cl(4, 1:num_pts)*100, 'g.-', 'LineWidth', 1.5 , 'markersize', 10);
% Labels
ylabel('U4 (cm)');
% Put the time label below the bottom graph.
xlabel('time (msec)');


