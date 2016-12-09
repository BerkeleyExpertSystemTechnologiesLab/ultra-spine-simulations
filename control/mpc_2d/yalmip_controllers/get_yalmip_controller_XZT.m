% get_yalmip_controller_XZT.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This function contains the formulation of the 2D spine MPC optimization
% problem (This is for a 2 vertebra system)

function [controller, constraints, objective, parameters_in, solution_out] = get_yalmip_controller_XZT ...
    (N, xi, u, xi_ref, u_ref, A_t, B_t, c_t, prev_in, spine_geo)

% Inputs:
% inputs, sdpvar of inputs to the system
% states, sdpvar of the system states
% A_t, sdpvar of linearized state matrix at time t
% B_t, sdpvar of linearized input matrix at time t
% c_t, sdpvar of linearized output matrix at time t
% prev_in, sdpvar of the input at the previous timestep
% reference, sdpvar of the reference trajectory

% Outputs:
% controller, the YALMIP controller object for this optimization problem. Is created from all the other outputs, supplied here for convenience.
% constraints, objective, parameters_in, solutions_out: all the other YALMIP variables. These are only output from this script for debugging purposes,
% the controller variable is all that's needed.

% Recall, the states for this system are
% x
% z
% T (theta)
% dx
% dz
% dT

% Quick error check to make sure horizon is positive and at least 2
assert(N > 1, 'Horizon must be at least length 2')

nx = size(xi,1);
nu = size(u,1);

u_lim_L = -ones(nu,1);
u_lim_U = ones(nu,1);
xi_lim_L = -ones(nx,1);
xi_lim_U = ones(nx,1);

%% Build constraints

constraints = [];
for i = 1:N
%     constraints = [constraints u_lim_L<=u(:,i)<=u_lim_U ...
%         xi_lim_L<=xi(:,i)<=xi_lim_U ...
%         xi(:,i+1)==A_t*xi(:,i)+B_t*u(:,i)+c_t ...
%         xi(2,i)>=spine_geo.h/2];
    
    constraints = [constraints xi(:,i+1)==A_t*xi(:,i)+B_t*u(:,i)+c_t ...
        xi(2,i)>=spine_geo.h/2];
end
constraints = [constraints xi(2,N+1)>=spine_geo.h/2];

%% Build objective

% Q = blkdiag(eye(nx/2),zeros(nx/2));
Q = eye(nx/2);
R = eye(nu);
% P = blkdiag(eye(nx/2),zeros(nx/2));
P = eye(nx/2);

objective = 0;
for i = 1:N
   objective = objective + (xi(1:3,i)-xi_ref(1:3,i))'*Q*(xi(1:3,i)-xi_ref(1:3,i)) ...
       + (u(:,i)-u_ref(:,i))'*R*(u(:,i)-u_ref(:,i));
end
objective = objective + (xi(1:3,N+1)-xi_ref(1:3,N+1))'*P*(xi(1:3,N+1)-xi_ref(1:3,N+1));

%% Build controller

parameters_in = {A_t B_t c_t xi_ref u_ref xi(:,1)};
solution_out = {xi u};
options = sdpsettings('verbose',1,'solver','gurobi');

controller = optimizer(constraints,objective,options,parameters_in,solution_out);

