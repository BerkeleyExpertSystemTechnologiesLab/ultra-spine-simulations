% mpc_XZ_plot_2d.m
% A quick script to make an x-z plot from Shirley's data.
% This was used by Drew for AHS 2017 to try to re-create Fig.9 from the IEEE T-CST paper, with larger text.

% prepare the workspace
clear all
close all
clc

% constants from the 3d error analysis script
fontsize = 12;

% load in the data we'll analyze: 200 pts, 0 to pi/8, 1e-6. 
% Should be 1e-5 but I can't find that one anywhere???
load affine_200_2e-5_eps_1e-6_0-pi-8;

% Plot the reference. Remove the last four states. N=4 here, right? That should be saved in the .mat file but I don't see it.
xz_handle = figure;
hold on;
set(xz_handle,'Position',[100,100,500,300]);
set(gca,'FontSize',fontsize);
plot( xi_traj(1,:)*100, xi_traj(2, :)*100, 'b','LineWidth',3);
plot( xi_cl(1,:)*100, xi_cl(2,:)*100, 'g','LineWidth',3);
legend('Reference', 'Result','Location','Best');
xlabel('Position in X (cm)');
ylabel('Position in Z (cm)');
title('Vertebra X-Z Position: Input Tracking');
set(gca,'FontSize',fontsize);