% linearize_dynamics.m
% Abishek Akella

function [A, B, c] = linearize_dynamics(xbar, ubar, restLengths, links, dt)
% Calculates linearized approximation to the dynamics given the current
% state and input to linearize around. Linearization is calculated
% numerically by determining model matrices A and B via the Jacobians of
% the nonlinear dynamics function x' = f(x, u).
%
% Inputs:
%   xbar = (12*links)x1 vector of states, organized from bottom link to top link
%   ubar = (8*links)x1 vector of inputs, organized from bottom link to top link
%   restLengths = given lengths of each of the cables in the resting
%       (equilibrium position)
%   links = number of links in spine
%   dt = timestep
%
% This function calls the following functions related to the spine
% dynamics:
%   simulate_dynamics

eps = .001;
A = zeros(12*links);

for k = 1:(12*links)
    xlinU = xbar;
    xlinL = xbar;
    xlinU(k) = xbar(k) + eps;
    xlinL(k) = xbar(k) - eps;
    A(:, k) = (reshape(simulate_dynamics(reshape(xlinU, 12, links)', restLengths, reshape(ubar, 8, links)', dt, links, 0)', 12*links, 1) ...
        - reshape(simulate_dynamics(reshape(xlinL, 12, links)', restLengths, reshape(ubar, 8, links)', dt, links, 0)', 12*links, 1))/(2*eps);
end

B = zeros(12*links, 8*links);
for k = 1:(8*links)
    ulinU = ubar;
    ulinL = ubar;
    ulinU(k) = ubar(k) + eps;
    ulinL(k) = ubar(k) - eps;
    B(:, k) = (reshape(simulate_dynamics(reshape(xbar, 12, links)', restLengths, reshape(ulinU, 8, links)', dt, links, 0)', 12*links, 1) ...
        - reshape(simulate_dynamics(reshape(xbar, 12, links)', restLengths, reshape(ulinL, 8, links)', dt, links, 0)', 12*links, 1))/(2*eps);
end

c = reshape(simulate_dynamics(reshape(xbar, 12, links)', restLengths, reshape(ubar, 8, links)', dt, links, 0)', 12*links, 1) - A*xbar - B*ubar;
end