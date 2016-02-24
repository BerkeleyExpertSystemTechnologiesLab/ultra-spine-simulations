% ultra_spine_mpc.m
% Abishek Akella, Andrew P. Sabelhaus
% This is the primary file for running the ULTRA Spine Model Predictive
% Control work.

%% Initialize system parameters
%clear variables;
clear all;
close all;
clc;

% The MATLAB functions that calculate the dynamics of the robot must be in
% the current path. If this file is being run outside the
% ultra-spine-simulations repository, change this path.
% Note that the files here are not directly called from this script, but
% instead are used in the following files (or those they call):
%   - simulate_reference_trajectory
%   - simulate_dynamics
%   - linearize_dynamics

path_to_dynamics = '../../dynamics/3d-dynamics-symbolicsolver';
addpath(path_to_dynamics);

% Time step for dynamics
dt = 0.001;

% Parameters for plotting:
% NOTE that the lengths and sizes here are only used for plotting, 
% since the dynamics used in this work is a point-mass model.
% Animation Frame Divisor
frame = 3;
 % Enable string plotting
stringEnable = 1;
% Radius of Link.
rad = 0.01;
% Gravitational force
g = 9.81;
% Height of Link
h = 0.15;
% Length of leg
l = 0.15;
% Projection of leg length
leg = (l^2 - (h/2)^2)^.5;

% Number of links in addition to base link.
% NOTE that this must be consistent with the dynamics defined in
% duct_accel.m and associated files! Those dynamics are 
links = 3;

 % Simulation time
time = 0:dt:500;

% To save a video, uncomment:
% - the following three initialization lines
% - the getframe call within the last loop
% - the open, save, and close lines at the end of this script
%videoObject = VideoWriter('../videos/UTLTRASpineMPC4TetraCircular.avi');
%videoObject.Quality = 100;
%videoObject.FrameRate = 5;

%% Initialize Plot
Figs = figure('Units','Normalized', 'outerposition', [0 0 1 1]);
M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

cmaps = summer(512);
ax = axes();

grid on;
axis equal;
xlim([-0.25 0.25])
ylim([-0.25 0.25])
zlim([-0.1, 0.5])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
colormap(cmaps(1:256,:))
shading interp
light
lighting phong
hold on
view(3)
title('ULTRA Spine Model')
m = 256;

%% Simulation
restLengths(1) = 0.1; % Vertical cable rest length
restLengths(2) = 0.1;
restLengths(3) = 0.1;
restLengths(4) = 0.1;
restLengths(5) = 0.187; % Saddle cable rest length
restLengths(6) = 0.187;
restLengths(7) = 0.187;
restLengths(8) = 0.187;

systemStates = zeros(links, 12);

x_initial = [];
for k = 1:links
    x(k) = 0; y(k) = 0.0; z(k) = 0.1*k; T(k) = 0.0; G(k) = 0.0; P(k) = 0.0;
    dx(k) = 0; dy(k) = 0; dz(k) = 0; dT(k) = 0; dG(k) = 0; dP(k) = 0;
    
    Tetra{k} = [(leg + rad/2), 0, -h/2; -(leg + rad/2), 0, -h/2; 0, (leg + rad/2), h/2; 0, -(leg + rad/2), h/2];
    [transform{k}, handle{k}] = plotSpineLink(Tetra{k}, rad, ax);
    
    systemStates(k, :) = [x(k), y(k), z(k), T(k), G(k), P(k), dx(k), dy(k), dz(k), dT(k), dG(k), dP(k)];
    x_initial = [x_initial; x(k); y(k); z(k); T(k); G(k); P(k); dx(k); dy(k); dz(k); dT(k); dG(k); dP(k)];
end
Tetra{links+1} = Tetra{links};
[transform{links+1}, handle{links+1}] = plotSpineLink(Tetra{links+1}, rad, ax);

%% Reference Trajectory
L = 180;
theta = linspace(-pi, pi, L);
r = .04;
traj = [r*cos(theta) + r; ...
        zeros(1, L); ...
        r*sin(theta) + 0.3; ...
        zeros(1, L); ...
        zeros(1, L); ... 
        zeros(1, L); ...
        zeros(1, L); ...
        zeros(1, L); ...
        zeros(1, L); ...
        zeros(1, L); ...
        zeros(1, L); ...
        zeros(1, L)];
plot3(traj(1, :), zeros(1, L), traj(3, :), 'r', 'LineWidth', 2);

%% Controller Initialization
disp('Controller Initialization')
pause(0.5)
disp('unpaused')

N = 10; % Horizon length
inputs = sdpvar(repmat(8*links, 1, N-1), repmat(1, 1, N-1));
states = sdpvar(repmat(12*links, 1, N), repmat(1, 1, N));
A_t = sdpvar(repmat(12*links, 1, 12*links), repmat(1, 1, 12*links));
B_t = sdpvar(repmat(12*links, 1, 8*links), repmat(1, 1, 8*links));
c_t = sdpvar(36, 1);
prev_in = sdpvar(8*links, 1);
reference = sdpvar(repmat(12, 1, N), repmat(1, 1, N));

input_lim = .07*ones(24, 1); % Limit on length of cable allowed

constraints = [norm(inputs{1} - prev_in, inf) <= 0.02]; % Deviation from previous applied input to current input
for k = 1:(N-2)
    constraints = [constraints, states{k+1} == [A_t{:}]*states{k} + [B_t{:}]*inputs{k} + c_t, ...
        -input_lim <= inputs{k} <= input_lim, ...
        norm(inputs{k}(1:8) - inputs{1}(1:8), inf) <= 0.01, ... % Minimize deviation from first input (minimize linearization error)
        norm(inputs{k}(9:16) - inputs{1}(9:16), inf) <= 0.01, ...
        norm(inputs{k}(17:24) - inputs{1}(17:24), inf) <= 0.01];
end
constraints = [constraints, norm(inputs{N-1}(1:8) - inputs{1}(1:8), 2) <= 0.1, ...% Final input is given a wider tolerance
    norm(inputs{N-1}(9:16) - inputs{1}(9:16), 2) <= 0.1, ...
    norm(inputs{N-1}(17:24) - inputs{1}(17:24), 2) <= 0.1];
constraints = [constraints, states{N} == [A_t{:}]*states{N-1} + [B_t{:}]*inputs{N-1} + c_t, -input_lim <= inputs{N-1} <= input_lim];

for j = 1:(N-1)
    constraints = [constraints, norm(states{j}(1:6) - states{j+1}(1:6), inf) <= 0.02, ... % Minimize deviation from first state (minimize linearization error)
        norm(states{j}(13:18) - states{j+1}(13:18), inf) <= 0.03, ...
        norm(states{j}(25:30) - states{j+1}(25:30), inf) <= 0.04, ...
        states{j}(3) + .02 <= states{j}(15), ... % Maintain some distance between links to prevent collision
        states{j}(15) + .02 <= states{j}(27)];
end
constraints = [constraints, states{N}(3) + .02 <= states{N}(15), states{N}(15) + .02 <= states{N}(27)];

objective = 25*norm(states{1}(25:27) - reference{1}(1:3), 2); % Minimize deviations along trajectory
for k = 2:(N-1)
    objective = objective + (1/2)*(25^k)*norm(states{k}(25:27) - reference{k}(1:3), 2) + (1/24)*(3^k)*norm(inputs{k} - inputs{k-1}, inf) ...
        + (3^k)*(norm(states{k}(25:27) - states{k-1}(25:27))); % Also minimize change in state/input for smooth motion
end
objective = objective + (1/2)*(25^N)*norm(states{N}(25:27) - reference{N}(1:3), 2) + (3^N)*norm(states{N}(25:27) - states{N-1}(25:27));

parameters_in = {prev_in, states{1}, [A_t{:}], [B_t{:}], c_t, [reference{:}]};
solutions_out = {[inputs{:}], [states{:}]};

% Build controller object for faster computation during iteration
controller = optimizer(constraints, objective, sdpsettings('solver', 'gurobi', 'verbose', 1), parameters_in, solutions_out);

%% Forward simulate trajectory
disp('Forward simulate trajectory')
pause(0.5)
disp('unpaused')

[x_ref, u_ref, M] = simulate_reference_traj(controller, systemStates, restLengths, links, dt, x, y, z, T, G, P, ...
    dx, dy, dz, dT, dG, dP, traj, N);

refx = [x_ref{:}];
plot3(refx(25, :), refx(26, :), refx(27, :), 'b-.', 'LineWidth', 2);

%% Build Iterative LQR Controller
disp('Build Iterative LQR Controller')
pause(0.5)
disp('unpaused')

Q = zeros(12);
Q(1:6, 1:6) = eye(6);
Q_lqr = 5*kron(eye(3), Q);
R_lqr = 5*eye(8*links);

tic;
P0 = zeros(36);
[A, B, ~] = linearize_dynamics(x_ref{M}, u_ref{M}, restLengths, links, dt);
K{1} = -((R_lqr + B'*P0*B)^-1)*B'*P0*A;
P_lqr{1} = Q_lqr + K{1}'*R_lqr*K{1} + (A + B*K{1})'*P0*(A + B*K{1});
for k = (M-1):-1:1
    disp(strcat('Controller Build iteration:',num2str(k)))
    [A, B, ~] = linearize_dynamics(x_ref{k}, u_ref{k}, restLengths, links, dt);
    K{M-k+1} = -((R_lqr + B'*P_lqr{M-k}*B)^-1)*B'*P_lqr{M-k}*A;
    P_lqr{M-k+1} = Q_lqr + K{M-k+1}'*R_lqr*K{M-k+1} + (A + B*K{M-k+1})'*P_lqr{M-k}*(A + B*K{M-k+1});
end
toc;

disp('Starting simulation')
%% Plotting
if (stringEnable)    
    anchor1=[0 0 rad];
    anchor2=[0 0 rad];
    String_pts = [];
    for k = 1:links
       String_pts = [String_pts; (Tetra{k}(1,:)+anchor1); (Tetra{k+1}(1,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(3,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); ...
                    (Tetra{k}(2,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(4,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(1,:)-anchor2)];
    end
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','m');
    linkdata on
end

plot_dt = 0.01; % Slow down the animation slightly
offset = 30; % Used for animation so it stays still for a few seconds

for k = 1:links
    x(k) = 0; y(k) = 0.0; z(k) = 0.1*k; T(k) = 0.0; G(k) = 0.0; P(k) = 0.0;
    dx(k) = 0; dy(k) = 0; dz(k) = 0; dT(k) = 0; dG(k) = 0; dP(k) = 0;
end

s = 1;
for t = 1:((M-1)+offset)
    String_pts = [];
    if mod(t, frame) == 0
        tic
        
        while toc < (plot_dt*frame)
            % Wait until time has passed
        end
        
        for k = 1:links
            RR{k} =  getHG_Tform(x(k),y(k),z(k),T(k),G(k),P(k)); % Build graphical model of each link
            set(transform{k},'Matrix',RR{k});
        end
        
        if (stringEnable)
            for k = 2:(links+1)
                Tetra{k} = [(l^2 - (h/2)^2)^.5, 0, -h/2, 1; -(l^2 - (h/2)^2)^.5, 0, -h/2, 1; 0, (l^2 - (h/2)^2)^.5, h/2, 1; 0, -(l^2 - (h/2)^2)^.5, h/2, 1];
                Tetra{k} = RR{k-1}*Tetra{k}';
                Tetra{k} = Tetra{k}';
                Tetra{k} = Tetra{k}(:,1:3);
            end
            anchor2 = [0 0 rad];
            for k = 1:links
                String_pts = [String_pts; (Tetra{k}(1,:)+anchor1); (Tetra{k+1}(1,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(3,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); ...
                    (Tetra{k}(2,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(4,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(1,:)-anchor2)];
            end
            refreshdata(string_handle);
        end
        drawnow;
    end
    
    for k = 1:links
        systemStates(k, 1) = x(k); systemStates(k, 2) = y(k); systemStates(k, 3) = z(k);
        systemStates(k, 4) = T(k); systemStates(k, 5) = G(k); systemStates(k, 6) = P(k);
        systemStates(k, 7) = dx(k); systemStates(k, 8) = dy(k); systemStates(k, 9) = dz(k);
        systemStates(k, 10) = dT(k); systemStates(k, 11) = dG(k); systemStates(k, 12) = dP(k);
    end
    
    %% Controller 
    noise = 0;
    if ((t > offset) && (t < M + offset))
        control = K{M+offset-t}*(reshape(systemStates', 36, 1) - x_ref{t-offset}) + u_ref{t-offset};
        systemStates = simulate_dynamics(systemStates, restLengths, reshape(control, 8, 3)', dt, links, noise);
        actual_traj(:, s) = systemStates(links, :); s = s + 1;
    end
    
    %% End Controller Design
    for k = 1:links
        x(k) = systemStates(k, 1); y(k) = systemStates(k, 2); z(k) = systemStates(k, 3);
        T(k) = systemStates(k, 4); G(k) = systemStates(k, 5); P(k) = systemStates(k, 6);
        dx(k) = systemStates(k, 7); dy(k) = systemStates(k, 8); dz(k) = systemStates(k, 9);
        dT(k) = systemStates(k, 10); dG(k) = systemStates(k, 11); dP(k) = systemStates(k, 12);
    end
    
    % Record this frame for a video
    videoFrames(t) = getframe(gcf);
end
plot3(actual_traj(1, :), actual_traj(2, :), actual_traj(3, :), 'g', 'LineWidth', 2);

% Uncomment these lines to save the video.
%open(videoObject);
%writeVideo(videoObject, videoFrames);
%close(videoObject);
