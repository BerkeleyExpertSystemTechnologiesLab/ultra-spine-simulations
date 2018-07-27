%% Analyse the errors of the MPC results

x_1 = abs(xi_traj(1,1)*ones(1,size(xi_cl,2))-xi_traj(1,1:end-opt_params.horizon_length));
x_2 = abs(xi_traj(2,1)*ones(1,size(xi_cl,2))-xi_traj(2,1:end-opt_params.horizon_length));
x_3 = abs(xi_traj(3,1)*ones(1,size(xi_cl,2))-xi_traj(3,1:end-opt_params.horizon_length));

% % Prior: Use absolue state value as x, which cannot reflect absolute value change of
% z_state, as it starts from 0.1 in case of 0 initial sweeping angle
% x_1 = abs(xi_cl(1,:));
% x_2 = abs(xi_cl(2,:));
% x_3 = abs(xi_cl(3,:));

% Calculate absolute errors of states
% x_1 = xi_cl(1,:)-xi_traj(1,1:length(xi_cl(1,:)));
e_1 = xi_cl(1,:)-xi_traj(1,1:end-opt_params.horizon_length);
e_2 = xi_cl(2,:)-xi_traj(2,1:end-opt_params.horizon_length);
e_3 = xi_cl(3,:)-xi_traj(3,1:end-opt_params.horizon_length);


% % TEST for trajectory lag:
% e_1p = xi_cl(1,:)/0.9 - xi_traj(1,1:end-opt_params.horizon_length);
% e_2p = (xi_cl(2,1)*ones(1,size(xi_cl,2)) - (xi_cl(2,1)*ones(1,size(xi_cl,2)) - xi_cl(2,:))/0.8) - xi_traj(2,1:end-opt_params.horizon_length);
% e_3p = xi_cl(3,:)/0.9 - xi_traj(3,1:end-opt_params.horizon_length);

% x_sq_1 = zeros(1,length(x_1));
% x_sq_2 = zeros(1,length(x_1));
% x_sq_3 = zeros(1,length(x_1));
% for i = 1:length(x_1)
%     x_sq_1(i) = x_1(i)^2;
%     x_sq_2(i) = x_2(i)^2;
%     x_sq_3(i) = x_3(i)^2;
% end

% Calculate square errors of states
e_sq_1 = e_1.*e_1;
e_sq_2 = e_2.*e_2;
e_sq_3 = e_3.*e_3;

% e_sq = e_1.*e_1+e_2.*e_2;

% Calculate absolute errors of states
e_abs_1 = abs(e_1);
e_abs_2 = abs(e_2);
e_abs_3 = abs(e_3);

e_abs_1p = abs(e_1p);
e_abs_2p = abs(e_2p);
e_abs_3p = abs(e_3p);
 
% Calculate relative errors of states
e_rl_1 = e_abs_1./x_1;
e_rl_2 = e_abs_2./x_2;
e_rl_3 = e_abs_3./x_3;

e_rl_1p = e_abs_1p./x_1;
e_rl_2p = e_abs_2p./x_2;
e_rl_3p = e_abs_3p./x_3;

% Plot the relevat errors of each state, namely the absolute errors over absolue velues
figure;
subplot(3,1,1)
plot(100*e_rl_1,'.-','LineWidth',1.25);
grid;
ylabel('Relative E(x)/ %');
title('Relative Tracking Errors of States')
subplot(3,1,2)
plot(100*e_rl_2,'.-','LineWidth',1.25);
ylabel('Relative E(z)/ %');
grid;
subplot(3,1,3)
plot(100*e_rl_3,'.-','LineWidth',1.25);
grid;
ylabel('Relative E(theta)/ %')

% % Plot the square errors of each state 
% figure;
% subplot(3,1,1)
% plot(e_sq_1);
% ylabel('SE(x)');
% title('Square Errors of States')
% subplot(3,1,2)
% plot(e_sq_2);
% ylabel('SE(z)');
% subplot(3,1,3)
% plot(e_sq_3);
% ylabel('SE(theta)')

% Plot the absolut error of 3 states, if necessary
figure;
subplot(3,1,1)
plot(e_abs_1);
ylabel('AE(x)');
title('Absolute Errors of States')
subplot(3,1,2)
plot(e_abs_2);
ylabel('AE(z)');
subplot(3,1,3)
plot(e_abs_3);
ylabel('AE(theta)')

%% 
figure;
subplot(3,1,1)
plot(100*e_rl_1p,'.-','LineWidth',1.25);
grid;
ylabel('Relative E(x)/ %');
title('Relative Tracking Errors of States')
subplot(3,1,2)
plot(100*e_rl_2p,'.-','LineWidth',1.25);
ylabel('Relative E(z)/ %');
grid;
subplot(3,1,3)
plot(100*e_rl_3p,'.-','LineWidth',1.25);
grid;
ylabel('Relative E(theta)/ %')
%% PLOT THE X-Z POSITION
figure;
plot(100*xi_traj(1,1:opt_params.num_pts),100*xi_traj(2,1:opt_params.num_pts),'-o','LineWidth',2.5);
% plot(100*xi_traj(1,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10),100*xi_traj(2,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10),'-o','LineWidth',2.5);
% plot(100*xi_traj(1,1:opt_params.num_pts+2-2*(opt_params.num_pts+1)/10),100*xi_traj(2,1:opt_params.num_pts+2-2*(opt_params.num_pts+1)/10),'LineWidth',2.5);
hold on;
plot(100*xi_cl(1,:), 100*xi_cl(2,:),'-x','LineWidth',2.5);
grid on;
hold on;
plot(100*xi_cl(1,1), 100*xi_cl(2,1),'o','LineWidth',3.5);
xlabel('X /cm');
ylabel('Z /cm');
% ylim([0,10]);
title('Plot of X-Z Position');
legend('reference','tractory','start point','location','best');
%% PLOT THE ANGULAR TRAJECTORY
figure;
plot(xi_traj(3,1:opt_params.num_pts+1)*180/pi,'-o','LineWidth',2.5);
% plot(xi_traj(3,1:opt_params.num_pts+2-(opt_params.num_pts+1)/10)*180/pi,'-o','LineWidth',2.5);
% plot(xi_traj(3,1:opt_params.num_pts+1-1.5*(opt_params.num_pts+1)/10)*180/pi,'LineWidth',2.5);
% plot(xi_cl(3,:), '-x','LineWidth',2.5);
hold on;
plot(xi_cl(3,:)*180/pi, '-x','LineWidth',2.5);
grid on;
hold on;
% plot(xi_cl(3,1),'o','LineWidth',3.5);
plot(xi_cl(3,1)*180/pi,'o','LineWidth',3.5);
xlabel('steps');
% ylabel(' \theta /arc');
ylabel(' \theta /°');
title('Plot of \theta over steps');
legend('reference','tractory','start point','location','best');
%% PLOT X-Y-Theta IN 3D
figure;
plot3(100*xi_cl(1,:), 100*xi_cl(2,:),180/pi*xi_cl(3,:),'-x','LineWidth',2.5);
hold on;
plot3(100*xi_traj(1,1:opt_params.num_pts+1), ...
    100*xi_traj(2,1:opt_params.num_pts+1), ...
    180/pi*xi_traj(3,1:opt_params.num_pts+1), ...
    '-o','LineWidth',1.5);
% plot3(100*xi_traj(1,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10), ...
%     100*xi_traj(2,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10), ...
%     180/pi*xi_traj(3,1:opt_params.num_pts+1-(opt_params.num_pts+1)/10), ...
%     '-o','LineWidth',2.5);
grid on;
hold on;
plot3(100*xi_cl(1,1), 100*xi_cl(2,1),180/pi*xi_cl(3,1),'o','LineWidth',3.5);
xlabel('X /cm');
ylabel('Z /cm');
zlabel(' \theta /°');
title('Plot of states X - Z - \theta');
legend('tractory','reference','start point','location','best');
