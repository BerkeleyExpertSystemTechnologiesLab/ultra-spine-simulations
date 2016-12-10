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

u_lim_L = zeros(nu,1);
u_lim_U = 10*ones(nu,1);
xi_lim_L = -2*ones(nx,1);
xi_lim_U = 2*ones(nx,1);

u_norm_lim = 1;
pos_norm_lim = 1;
vel_norm_lim = 1;

%% Build constraints

constraints = [];
for i = 1:N
    %     constraints = [constraints u_lim_L<=u(:,i)<=u_lim_U ...
    %         xi_lim_L<=xi(:,i)<=xi_lim_U ...
    %         xi(:,i+1)==A_t*xi(:,i)+B_t*u(:,i)+c_t ...
    %         xi(2,i)>=spine_geo.h/2];
    
    constraints = [constraints xi(:,i+1)==A_t*xi(:,i)+B_t*u(:,i)+c_t ...
        xi(2,i)>=spine_geo.h/2 ...
        u_lim_L<=u(:,i)];
end
% 
% for i = 1:N-1
%     constraints = [constraints norm(u(:,i+1)-u(:,i))<=u_norm_lim ...
%         norm(xi(1:3,i+1)-xi(1:3,i),inf)<=pos_norm_lim ...
%         norm(xi(4:6,i+1)-xi(4:6,i),inf)<=vel_norm_lim];
% end

constraints = [constraints xi(2,N+1)>=spine_geo.h/2];

% constraints = [constraints ...
%     norm(xi(1:3,N+1)-xi(1:3,N),inf)<=pos_norm_lim ...
%     norm(xi(4:6,N+1)-xi(4:6,N),inf)<=vel_norm_lim];

%% Build objective

% Q = blkdiag(eye(nx/2),zeros(nx/2));
Q = eye(nx/2);
R = eye(nu);
% P = blkdiag(eye(nx/2),zeros(nx/2));
P = eye(nx/2);

objective = (xi(1:3,1)-xi_ref(1:3,1))'*Q*(xi(1:3,1)-xi_ref(1:3,1)) ...
        + (u(:,1)-u_ref(:,1))'*R*(u(:,1)-u_ref(:,1));
for i = 2:N
        objective = objective + (xi(1:3,i)-xi_ref(1:3,i))'*Q*(xi(1:3,i)-xi_ref(1:3,i)) ...
        + (u(:,i)-u_ref(:,i))'*R*(u(:,i)-u_ref(:,i));
%     objective = objective + (xi(1:3,i)-xi_ref(1:3,i))'*Q*(xi(1:3,i)-xi_ref(1:3,i)) ...
%         + (u(:,i)-u_ref(:,i))'*R*(u(:,i)-u_ref(:,i)) ...
%         + (1/2)*norm(xi(1:3,i)-xi(1:3,i-1)) ...
%         + (1/24)*(3^i)*norm(u(:,i)-u(:,i-1),inf);
%     objective = objective + (xi(:,i)-xi_ref(:,i))'*Q*(xi(:,i)-xi_ref(:,i)) ...
%         + (u(:,i)-u_ref(:,i))'*R*(u(:,i)-u_ref(:,i));
end
objective = objective + (xi(1:3,N+1)-xi_ref(1:3,N+1))'*P*(xi(1:3,N+1)-xi_ref(1:3,N+1));
% objective = objective + (xi(1:3,N+1)-xi_ref(1:3,N+1))'*P*(xi(1:3,N+1)-xi_ref(1:3,N+1)) ...
%     + (3^i)*norm(xi(1:3,N+1)-xi(1:3,N));
% objective = objective + (xi(:,N+1)-xi_ref(:,N+1))'*P*(xi(:,N+1)-xi_ref(:,N+1));

%% Build controller

parameters_in = {A_t B_t c_t xi_ref u_ref xi(:,1)};
% parameters_in = {A_t B_t c_t xi_ref xi(:,1)};
solution_out = {xi u};
options = sdpsettings('verbose',1,'solver','gurobi');

controller = optimizer(constraints,objective,options,parameters_in,solution_out);

