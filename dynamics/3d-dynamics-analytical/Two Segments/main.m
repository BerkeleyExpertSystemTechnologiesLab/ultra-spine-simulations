% @Author ChanWoo Yang
% UC Berkeley
% BEST Lab
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Created 4/03/2015
% Modified 6/09/2015
% Contact ChanWoo at: chanwoo.yang@berkeley.edu
% Tensegrity Spine Dynamics: Two Stellated Tetrahedron Segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

dt = 0.02;      % Time step [sec]
Tf = 15;         % Final Time Step [sec]

alpha = 109.5;       % angle b/w rods [degree]
R = 0.2;            % Length of Rod [m]
% Distance on 2D projection
l = R*sind(alpha/2);
h = R*cosd(alpha/2);
m = 0.3;              % mass of node [kg]
g = 9.81;                % gravitational acc. [m/s^2]

yy10 = [0 0 0.15 20 0 30 0 0 0 0 0 0];    % Initial Condition


cs = 10;             % Saddle cable damping coefficient [Ns/m]
cv = 10;             % Vertical cable damping coefficient [Ns/m]
ks = 1220;          % Saddle cable spring constant [N/m]
kv = 1220;          % Vertical cable spring constant [N/m]
Ls0 = 0.1;          % Saddle cable initial cable length [m] 
Lv0 = 0.1;          % Vertical cable initial cable length [m]

% Fixed first segment
node1 = [0,0,0;
        -l,0,-h;
        l,0,-h;
        0,-l,h;
        0,l,h];


options = odeset('reltol',1.e-10,'abstol',1.e-10);
[T,Y] = ode45(@(t,yy)eomSolver(t,yy,node1,R,m,g,cv,cs,kv,ks,Ls0,Lv0),0:dt:Tf,yy10,options);

node = getNodeCoord(R,T,Y);   % node([x,y,z],node#,timeStep)

%% Plot Animation

% Uncomment these lines to save a video
% videoObject = VideoWriter('/Users/ChanWoo/Desktop/SpineDynamics.avi');  % Modify this accordingly
% videoObject.Quality = 100;
% videoObject.FrameRate = 5;

for i = 1:length(T)
    figure(1)
    plot3(node(1,1,i),node(2,1,i),node(3,1,i),'*b',...
          node(1,2,i),node(2,2,i),node(3,2,i),'*r',...
          node(1,3,i),node(2,3,i),node(3,3,i),'*k',...
          node(1,4,i),node(2,4,i),node(3,4,i),'*m',...
          node(1,5,i),node(2,5,i),node(3,5,i),'*y',...
          [node(1,1,i) node(1,2,i)],[node(2,1,i) node(2,2,i)],[node(3,1,i) node(3,2,i)],'-b',...
          [node(1,1,i) node(1,3,i)],[node(2,1,i) node(2,3,i)],[node(3,1,i) node(3,3,i)],'-b',...
          [node(1,1,i) node(1,4,i)],[node(2,1,i) node(2,4,i)],[node(3,1,i) node(3,4,i)],'-b',...
          [node(1,1,i) node(1,5,i)],[node(2,1,i) node(2,5,i)],[node(3,1,i) node(3,5,i)],'-b',...
          node1(1,1),node1(1,2),node1(1,3),'*r',...
          node1(2,1),node1(2,2),node1(2,3),'*b',...
          node1(3,1),node1(3,2),node1(3,3),'*k',...
          node1(4,1),node1(4,2),node1(4,3),'*m',...
          node1(5,1),node1(5,2),node1(5,3),'*y',...
          [node1(1,1) node1(2,1)],[node1(1,2) node1(2,2)],[node1(1,3) node1(2,3)],'-r',...
          [node1(1,1) node1(3,1)],[node1(1,2) node1(3,2)],[node1(1,3) node1(3,3)],'-r',...
          [node1(1,1) node1(4,1)],[node1(1,2) node1(4,2)],[node1(1,3) node1(4,3)],'-r',...
          [node1(1,1) node1(5,1)],[node1(1,2) node1(5,2)],[node1(1,3) node1(5,3)],'-r',...
          [node(1,2,i) node1(2,1)],[node(2,2,i) node1(2,2)],[node(3,2,i) node1(2,3)],'--k',...
          [node(1,3,i) node1(3,1)],[node(2,3,i) node1(3,2)],[node(3,3,i) node1(3,3)],'--k',...
          [node(1,4,i) node1(4,1)],[node(2,4,i) node1(4,2)],[node(3,4,i) node1(4,3)],'--k',...
          [node(1,5,i) node1(5,1)],[node(2,5,i) node1(5,2)],[node(3,5,i) node1(5,3)],'--k',...
          [node(1,2,i) node1(4,1)],[node(2,2,i) node1(4,2)],[node(3,2,i) node1(4,3)],'--k',...
          [node(1,2,i) node1(5,1)],[node(2,2,i) node1(5,2)],[node(3,2,i) node1(5,3)],'--k',...
          [node(1,3,i) node1(4,1)],[node(2,3,i) node1(4,2)],[node(3,3,i) node1(4,3)],'--k',...
          [node(1,3,i) node1(5,1)],[node(2,3,i) node1(5,2)],[node(3,3,i) node1(5,3)],'--k');
    grid on
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
%     view([0,90])
    axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
    
    % Uncomment to save a video
%     M(i) = getframe(gcf);
    
    pause(0.001)
    T(i)
end

% Save the movie we generated
% Uncomment these 3 lines to save a video
% open(videoObject);
% writeVideo(videoObject,M);
% close(videoObject);

%% Compute engineering parameters

% Compute Angular Momentum
for i = 1:length(T)
r2(:,i) = node(:,2,i)-node(:,1,i);
r3(:,i) = node(:,3,i)-node(:,1,i);
r4(:,i) = node(:,4,i)-node(:,1,i);
r5(:,i) = node(:,5,i)-node(:,1,i);
% MoI(:,:,i) = m*r2(:,i)*r2(:,i)'+m*r3(:,i)*r3(:,i)'+...
%              m*r4(:,i)*r4(:,i)'+m*r5(:,i)*r5(:,i)';  %Moment of Inertia
% AM(:,i) = MoI(:,:,i)*Y(i,10:12)';
AM(:,i) = m*r2(:,i)*r2(:,i)'*Y(i,10:12)'+m*r3(:,i)*r3(:,i)'*Y(i,10:12)'+...
          m*r4(:,i)*r4(:,i)'*Y(i,10:12)'+m*r5(:,i)*r5(:,i)'*Y(i,10:12)';
end

% Compute spring forces
[Fv22,Fv33,Fv44,Fv55,Fs24,Fs25,Fs34,Fs35] = springForceChecker(node1,R,kv,ks,Ls0,Lv0,T,Y);
for i=1:length(T)
    Fv22Mag(i) = norm(Fv22(:,i)); 
    Fv33Mag(i) = norm(Fv33(:,i)); Fv44Mag(i) = norm(Fv44(:,i)); Fv55Mag(i) = norm(Fv55(:,i));
    Fs24Mag(i) = norm(Fs24(:,i)); Fs25Mag(i) = norm(Fs25(:,i)); Fs34Mag(i) = norm(Fs34(:,i)); Fs35Mag(i) = norm(Fs35(:,i));
end

%% Plot
%Vertical Cables
figure()
plot(T(:),Fv22Mag(:),'*k',T(:),Fv33Mag(:),'<r',T(:),Fv44Mag(:),'^b',T(:),Fv55Mag(:),'g')
legend('Vertical 12-22','Vertical 13-23','Vertical 14-24','Vertical 15-25','Location','best')
xlabel('Simulation Timestep [sec]','FontSize',14)
ylabel('Force [N]','FontSize',14)
grid on
%Saddle Cables
figure()
plot(T(:),Fs24Mag(:),'*k',T(:),Fs25Mag(:),'<r',T(:),Fs34Mag(:),'^b',T(:),Fs35Mag(:),'g')
legend('Saddle 22-14','Saddle 22-15','Saddle 23-14','Saddle 23-15','Location','best')
xlabel('Simulation Timestep [sec]','FontSize',14)
ylabel('Force [N]','FontSize',14)
grid on

% Plot total energy over time step
[totalE,ke,pe] = energyConsvChecker(node1,m,g,R,kv,ks,Ls0,Lv0,T,Y);
figure()
plot(T(:),totalE(:),'k',T(:),ke(:),'r',T(:),pe(:),'b')
legend('Total Energy','Kinetic Energy','Potential Energy','Location','best')
xlabel('Simulation Timestep [sec]','FontSize',14)
ylabel('Energy [J]','FontSize',14)
grid on

% Plot changes in x,y,z coordinate of the center node
figure()
plot(T(:),Y(:,1),'r',T(:),Y(:,2),'b',T(:),Y(:,3),'k')
legend('x1','y1','z1','Location','best')
xlabel('Simulation Timestep [sec]','FontSize',14)
ylabel('Distance Covered [m]','FontSize',14)
grid on

% Plot changes in theta,phi,psi of the structure
figure()
plot(T(:),Y(:,4),'r',T(:),Y(:,5),'b',T(:),Y(:,6),'k')
legend('theta','phi','psi','Location','best')
xlabel('Simulation Timestep [sec]','FontSize',14)
ylabel('Angle [deg]','FontSize',14)
grid on

% Plot Angular Momentum
figure()
plot(T(:),AM(1,:),'r',T(:),AM(2,:),'b',T(:),AM(3,:),'k')
legend('theta','phi','psi','Location','best')
xlabel('Simulation Timestep [sec]','FontSize',14)
ylabel('Angular Momentum','FontSize',14)
grid on

% %Saddle Cables z-direction Force
% figure()
% plot(T(:),Fs24(3,:),'*-k',T(:),Fs25(3,:),'<r',T(:),Fs34(3,:),'^b',T(:),Fs35(3,:),'g')
% legend('Saddle 22-14','Saddle 22-15','Saddle 23-14','Saddle 23-15','Location','best')
% xlabel('Simulation Timestep [sec]','FontSize',14)
% ylabel('Force [N]','FontSize',14)
% grid on
% 
% %Vertical Cables z-direction Force
% figure()
% plot(T(:),Fv22(3,:),'*-k',T(:),Fv33(3,:),'<r',T(:),Fv44(3,:),'^b',T(:),Fv55(3,:),'g')
% legend('Vertical 12-22','Vertical 13-23','Vertical 14-24','Vertical 15-25','Location','best')
% xlabel('Simulation Timestep [sec]','FontSize',14)
% ylabel('Force [N]','FontSize',14)
% grid on