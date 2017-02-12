% gen_report_plots.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This script generates the figures for the report

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script initialization

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data

data_path = './data';
addpath(data_path)

with_ref = load('MPC_results_with_invkin_ref.mat');
without_ref = load('MPC_results_without_invkin_ref.mat');

% with_ref = load('MPC_results_morecomplex_with_invkin_ref.mat');
% without_ref = load('MPC_results_morecomplex_without_invkin_ref.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate plots

num_pts = with_ref.opt_params.num_pts;

%% 1. Error plot with and without input reference tracking

error_w_ref = abs(with_ref.xi_cl(1:3,:)-with_ref.xi_traj(1:3,1:num_pts+1));
error_wo_ref = abs(without_ref.xi_cl(1:3,:)-without_ref.xi_traj(1:3,1:num_pts+1));

figure;
subplot(3,1,1)
plot(1:num_pts+1,error_w_ref(1,:),'LineWidth',1)
hold on
plot(1:num_pts+1,error_wo_ref(1,:),'LineWidth',1)
xlim([0 num_pts+1])
grid on
ylabel('x')
title('Absolute Tracking Error')
subplot(3,1,2)
plot(1:num_pts+1,error_w_ref(2,:),'LineWidth',1)
hold on
plot(1:num_pts+1,error_wo_ref(2,:),'LineWidth',1)
grid on
xlim([0 num_pts+1])
ylabel('z')
subplot(3,1,3)
plot(1:num_pts+1,error_w_ref(3,:),'LineWidth',1)
hold on
plot(1:num_pts+1,error_wo_ref(3,:),'LineWidth',1)
grid on
xlim([0 num_pts+1])
ylabel('\theta')
legend('With u_{ref}','Without u_{ref}','Location','best')

%% 2. [x z] State trajectory plot with and without input reference tracking

state_ref = with_ref.xi_traj;
state_w_ref = with_ref.xi_cl;
state_wo_ref = without_ref.xi_cl;

figure;
plot(state_ref(1,:),state_ref(2,:),'LineWidth',1)
hold on
plot(state_w_ref(1,:),state_w_ref(2,:),'LineWidth',1)
plot(state_wo_ref(1,:),state_wo_ref(2,:),'LineWidth',1)
plot(state_ref(1,1),state_ref(2,1),'o','LineWidth',1)
plot(state_ref(1,end),state_ref(2,end),'p','LineWidth',1)
plot(state_w_ref(1,end),state_w_ref(2,end),'p','LineWidth',1)
plot(state_wo_ref(1,end),state_wo_ref(2,end),'p','LineWidth',1)
grid on
title('x-z Position Trajectory')
legend('\xi_{ref}','\xi_{traj,w}','\xi_{traj,wo}','Location','best')
xlabel('x');
ylabel('z');

%% 3. All state trajectories with and without input reference tracking

figure;
subplot(3,1,1)
plot(1:num_pts+1,state_ref(1,1:num_pts+1),'LineWidth',1)
hold on
plot(1:num_pts+1,state_w_ref(1,:),'LineWidth',1)
plot(1:num_pts+1,state_wo_ref(1,:),'LineWidth',1)
xlim([0 num_pts+1])
ylabel('x')
grid on
title('State Trajectories')
subplot(3,1,2)
plot(1:num_pts+1,state_ref(2,1:num_pts+1),'LineWidth',1)
hold on
plot(1:num_pts+1,state_w_ref(2,:),'LineWidth',1)
plot(1:num_pts+1,state_wo_ref(2,:),'LineWidth',1)
xlim([0 num_pts+1])
ylabel('z')
grid on
subplot(3,1,3)
plot(1:num_pts+1,state_ref(3,1:num_pts+1),'LineWidth',1)
hold on
plot(1:num_pts+1,state_w_ref(3,:),'LineWidth',1)
plot(1:num_pts+1,state_wo_ref(3,:),'LineWidth',1)
xlim([0 num_pts+1])
ylabel('\theta')
grid on
legend('State reference','With u_{ref}','Without u_{ref}','Location','best')

%% 4. All input trajectories

input_ref = with_ref.u_traj;
input_w_ref = with_ref.u_cl;

figure;
subplot(4,1,1)
plot(1:num_pts,input_ref(1,1:num_pts),'LineWidth',1)
hold on
plot(1:num_pts,input_w_ref(1,:),'LineWidth',1)
xlim([0 num_pts])
ylabel('u_1')
grid on
title('Input trajectories')
subplot(4,1,2)
plot(1:num_pts,input_ref(2,1:num_pts),'LineWidth',1)
hold on
plot(1:num_pts,input_w_ref(2,:),'LineWidth',1)
xlim([0 num_pts])
ylabel('u_2')
grid on
subplot(4,1,3)
plot(1:num_pts,input_ref(3,1:num_pts),'LineWidth',1)
hold on
plot(1:num_pts,input_w_ref(3,:),'LineWidth',1)
xlim([0 num_pts])
ylabel('u_3')
grid on
subplot(4,1,4)
plot(1:num_pts,input_ref(4,1:num_pts),'LineWidth',1)
hold on
plot(1:num_pts,input_w_ref(4,:),'LineWidth',1)
xlim([0 num_pts])
ylabel('u_4')
grid on
legend('Input reference','Input trajectory','Location','best')