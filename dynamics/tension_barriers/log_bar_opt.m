% log_bar_opt.m
% Logistic barrier optimization
% Copyright 2016 Andrew P. Sabelhaus
% Calculates optimum values of the parameters to the logistic function
% so that the potential energy of a logistic * cable tension = 
% the original potential energy as close as possible.

clear all;
close all;
clc;

% This uses YALMIP.
T0 = sdpvar(1);
% Choose some large k1 (the logistic slope.)
% to-do: how to optimize over this when its optimal
% value probably equals infinity?
k1 = 500;
% For the spring constant: choose some arbitrary k2.
% We use 2000 N/m in lots of tensegrity work.
k2 = 2000;

% For now, try this minimization at various values of x.
% Later: should this be a minimzation over the integral 
% of potential energy of all potential spring lengths?
% Or something else?
% 10 cm:
x = 0.01;

% The objective function has two parts:
% the potential energy for the ideal spring,
% subtracted from the potential energy for our new spring.
PE_ideal = 0.5*k2*x^2;
PE_approx = -x + (1/(k1*k2*x - k1*T0))*log(1+ exp(k1*k2*x - k1*T0));
% Will this solve?
optimize([], PE_ideal - PE_approx);
% Solves, but totally not correct. T0 = 20 is silly. 20 newtons of offset??

% To-do:
% Maybe this potential energy isn't correct, and we need to do it
% as a function of both x and \dot x, like the equations say.
% Or, maybe the solution to that integral of the exponential function
% isn't correct.
% Or, I've got some mixed-up negative signs somewhere.
% MOST IMPORTANTLY: need to think about what the objective function
% should look like. How should we incorporate x and \dot x?
% Should we just look for some ranges that are appropriate for the situation,
% then integrate over the error?
% Should we use some machine learning techniques, or Bayesian statistics, to 
% minimze some expected value of the error?


