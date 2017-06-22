% get_yalmip_controller_XZT.m
%
% Copyright 2016 Mallory Daly, Andrew P. Sabelhaus, Ellande Tang, Shirley
% Zhao, Edward Zhu, Berkeley Emergent Space Tensegrities Lab
%
% This function contains the formulation of the 2D spine MPC optimization
% problem (This is for a 2 vertebra system)

function [solution_out, opt_flag] = run_yalmip_opt_XZT_one_step_affine ...
    (N, xi, u, xi_ref, u_ref, xi_0, A_t, B_t, prev_in, spine_geo, opt_time_lim)

% Inputs:
% A_t, sdpvar of linearized state matrix at time t
% B_t, sdpvar of linearized input matrix at time t
% c_t, sdpvar of linearized output matrix at time t
% prev_in, sdpvar of the input at the previous timestep
% reference, sdpvar of the reference trajectory

% Outputs:

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
nu = size(u+1,1);

u_lim_L = zeros(nu,1);
u_lim_U = 0.3*ones(nu,1);
% xi_lim_L = -2*ones(nx,1);
% xi_lim_U = 2*ones(nx,1);

w1 = 0.02;
w2 = 0.01;
w3 = 0.03;
w4 = 0.04;
w5 = 3;
w6 = 1;

%% Build constraints

% Build system dynamics constraints, input constraints, and collision
% contstraints
constraints = [xi(:,1)==xi_0];
for i = 1:N
    %     constraints = [constraints u_lim_L<=u(:,i)<=u_lim_U ...
    %         xi_lim_L<=xi(:,i)<=xi_lim_U ...
    %         xi(:,i+1)==A_t*xi(:,i)+B_t*u(:,i)+c_t ...
    %         xi(2,i)>=spine_geo.h/2];
    
    constraints = [constraints xi(:,i+1)==A_t*xi(:,i)+B_t*u(:,i) ...
        xi(2,i)>=spine_geo.h/2 ...
        u_lim_L<=u(:,i)];
end
constraints = [constraints xi(2,N+1)>=spine_geo.h/2];

% Build constraints from ACC 2017 paper
% constraints = [constraints norm(u(:,1)-prev_in,inf)<=w1];
% for i = 1:N-1
%     constraints = [constraints norm(u(:,i+1)-u(:,i),inf)<=w2 ...
%         norm(xi(1:3,i+1)-xi(1:3,i),inf)<=w3 ...
%         norm(xi(4:6,i+1)-xi(4:6,i),inf)<=w4];
% end
% constraints = [constraints ...
%     norm(xi(1:3,N+1)-xi(1:3,N),inf)<=w3 ...
%     norm(xi(4:6,N+1)-xi(4:6,N),inf)<=w4];

%% Build objective

% Q = blkdiag(eye(nx/2),zeros(nx/2));
Q = eye(nx/2);
R = 2*eye(nu);
% P = blkdiag(eye(nx/2),zeros(nx/2));
P = eye(nx/2);

% Build state and input reference tracking objective
objective = (xi(1:3,1)-xi_ref(1:3,1))'*Q*(xi(1:3,1)-xi_ref(1:3,1)) ...
        + (u(:,1)-u_ref(:,1))'*R*(u(:,1)-u_ref(:,1));
for i = 2:N+1
    objective = objective + (xi(1:3,i)-xi_ref(1:3,i))'*Q*(xi(1:3,i)-xi_ref(1:3,i)) ...
        + (u(:,i)-u_ref(:,i))'*R*(u(:,i)-u_ref(:,i));
%     objective = objective + (xi(1:3,i)-xi_ref(1:3,i))'*Q*(xi(1:3,i)-xi_ref(1:3,i)) ...
%         + (u(:,i)-u_ref(:,i))'*R*(u(:,i)-u_ref(:,i)) ...
%         + w5*norm(xi(1:3,i)-xi(1:3,i-1)) ...
%         + w6*norm(u(:,i)-u(:,i-1),inf);
%     objective = objective + (xi(1:3,i)-xi_ref(1:3,i))'*Q*(xi(1:3,i)-xi_ref(1:3,i)) ...
%         + u(:,i)'*R*u(:,i) ...
%         + (1/2)*norm(xi(1:3,i)-xi(1:3,i-1)) ...
%         + (1/24)*(3^i)*norm(u(:,i)-u(:,i-1));
end
objective = objective + (xi(1:3,N+1)-xi_ref(1:3,N+1))'*P*(xi(1:3,N+1)-xi_ref(1:3,N+1));
% objective = objective + (xi(1:3,N+1)-xi_ref(1:3,N+1))'*P*(xi(1:3,N+1)-xi_ref(1:3,N+1)) ...
%     + (3^i)*norm(xi(1:3,N+1)-xi(1:3,N));

%% Optimize

options = sdpsettings('verbose',1,'gurobi.TimeLimit',opt_time_lim,'solver','gurobi');

opt_flag = optimize(constraints,objective,options);

solution_out = {value(xi), value(u)};
