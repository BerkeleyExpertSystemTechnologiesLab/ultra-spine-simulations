% @Author ChanWoo Yang
% UC Berkeley
% BEST Lab
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Created 4/03/2015
% Modified 6/04/2015
% Contact ChanWoo at: chanwoo.yang@berkeley.edu
% Tensegrity Spine Dynamics: Two Stellated Tetrahedron Segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ F ] = appliedForce(node1,R,m,g,cv,cs,kv,ks,Ls0,Lv0,Fu1,Fu2,Fu3,Fu4,Fu5,yy)
%INPUT:
% node1 = Coordinates of fixed first segment
% l = R*sind(alpha/2);
% h = R*cosd(alpha/2);
% Fu1,Fu2,Fu3,Fu4,Fu5 = Applied force at each node
% x,y,z,theta,phi,psi = 6 states of a segment
%OUTPUT:
% F = Force vector; 
%     Summation of forces acting on each node respect to each state

x = yy(1); y = yy(2); z = yy(3); theta = yy(4); phi = yy(5); psi = yy(6);
xd = yy(7); yd = yy(8); zd = yy(9);
thetad = deg2rad(yy(10)); phid = deg2rad(yy(11)); psid = deg2rad(yy(12)); % Angular Velocity in [rad/s]
    

alpha = 109.5;    % angle b/w rods [degree]
% Distance on 2D projection
l = R*sind(alpha/2);
h = R*cosd(alpha/2);

F = zeros(6,1);
node2 = getNodeCoord(R,1,[x,y,z,theta,phi,psi])';     % Coordinate of moving segment

%Gravitational Force
Fm = [0 0 -m*g];

% Spring Forces at each node
Fs2 = kv*(signChecker(norm(node1(2,:)-node2(2,:))-Lv0))*((node1(2,:)-node2(2,:))/norm(node1(2,:)-node2(2,:)))+...
      ks*(signChecker(norm(node1(4,:)-node2(2,:))-Ls0))*((node1(4,:)-node2(2,:))/norm(node1(4,:)-node2(2,:)))+...
      ks*(signChecker(norm(node1(5,:)-node2(2,:))-Ls0))*((node1(5,:)-node2(2,:))/norm(node1(5,:)-node2(2,:)));
Fs3 = kv*(signChecker(norm(node1(3,:)-node2(3,:))-Lv0))*((node1(3,:)-node2(3,:))/norm(node1(3,:)-node2(3,:)))+...
      ks*(signChecker(norm(node1(4,:)-node2(3,:))-Ls0))*((node1(4,:)-node2(3,:))/norm(node1(4,:)-node2(3,:)))+...
      ks*(signChecker(norm(node1(5,:)-node2(3,:))-Ls0))*((node1(5,:)-node2(3,:))/norm(node1(5,:)-node2(3,:)));
Fs4 = kv*(signChecker(norm(node1(4,:)-node2(4,:))-Lv0))*((node1(4,:)-node2(4,:))/norm(node1(4,:)-node2(4,:)));
Fs5 = kv*(signChecker(norm(node1(5,:)-node2(5,:))-Lv0))*((node1(5,:)-node2(5,:))/norm(node1(5,:)-node2(5,:)));

% For no damping, comment below 
% Damping Forces at each node ---------------------------------------------
v1 = [xd yd zd]';

t2 = [-l 0 -h]';
t3 = [l 0 -h]';
t4 = [0 -l h]';
t5 = [0 l h]';
Tx = [1 0 0;
      0 cosd(theta) sind(theta);
      0 -sind(theta) cosd(theta)];
Ty = [cosd(phi) 0 sind(phi);
      0 1 0;
      -sind(phi) 0 cosd(phi)];
Tz = [cosd(psi) sind(psi) 0;
      -sind(psi) cosd(psi) 0;
      0 0 1];
Txd = [0 0 0; 
       0 -sind(theta)*thetad cosd(theta)*thetad;
       0 -cosd(theta)*thetad -sind(theta)*thetad];
Tyd = [-sind(phi)*phid 0 cosd(phi)*phid;
       0 0 0;
       -cosd(phi)*phid 0 -sind(phi)*phid];
Tzd = [-sind(psi)*psid cosd(psi)*psid 0;
       -cosd(psi)*psid -sind(psi)*psid 0;
       0 0 0];

v2 = v1 + (Txd*Ty*Tz + Tx*Tyd*Tz + Tx*Ty*Tzd)*t2;
v3 = v1 + (Txd*Ty*Tz + Tx*Tyd*Tz + Tx*Ty*Tzd)*t3;
v4 = v1 + (Txd*Ty*Tz + Tx*Tyd*Tz + Tx*Ty*Tzd)*t4;
v5 = v1 + (Txd*Ty*Tz + Tx*Tyd*Tz + Tx*Ty*Tzd)*t5;

%With Damping
Fd2 = cv*(((0-v2)'*(node1(2,:)-node2(2,:))')/norm(node1(2,:)-node2(2,:)))*((node1(2,:)-node2(2,:))/norm(node1(2,:)-node2(2,:)))+...
      cs*(((0-v2)'*(node1(4,:)-node2(2,:))')/norm(node1(4,:)-node2(2,:)))*((node1(4,:)-node2(2,:))/norm(node1(4,:)-node2(2,:)))+...
      cs*(((0-v2)'*(node1(5,:)-node2(2,:))')/norm(node1(5,:)-node2(2,:)))*((node1(5,:)-node2(2,:))/norm(node1(5,:)-node2(2,:)));
Fd3 = cv*(((0-v3)'*(node1(3,:)-node2(3,:))')/norm(node1(3,:)-node2(3,:)))*((node1(3,:)-node2(3,:))/norm(node1(3,:)-node2(3,:)))+...
      cs*(((0-v3)'*(node1(4,:)-node2(3,:))')/norm(node1(4,:)-node2(3,:)))*((node1(4,:)-node2(3,:))/norm(node1(4,:)-node2(3,:)))+...
      cs*(((0-v3)'*(node1(5,:)-node2(3,:))')/norm(node1(5,:)-node2(3,:)))*((node1(5,:)-node2(3,:))/norm(node1(5,:)-node2(3,:)));
Fd4 = cv*(((0-v4)'*(node1(4,:)-node2(4,:))')/norm(node1(4,:)-node2(4,:)))*((node1(4,:)-node2(4,:))/norm(node1(4,:)-node2(4,:)));
Fd5 = cv*(((0-v5)'*(node1(5,:)-node2(5,:))')/norm(node1(5,:)-node2(5,:)))*((node1(5,:)-node2(5,:))/norm(node1(5,:)-node2(5,:)));

%No Damping
% Fd2 = [0 0 0];
% Fd3 = [0 0 0];
% Fd4 = [0 0 0];
% Fd5 = [0 0 0];
% -----------------------------------------------

F1 = Fu1 + Fm;
F2 = Fu2 + Fs2 + Fm + Fd2;
F3 = Fu3 + Fs3 + Fm + Fd3;
F4 = Fu4 + Fs4 + Fm + Fd4;
F5 = Fu5 + Fs5 + Fm + Fd5;

% ------------------------------------------------
drdx = [1 0 0]';
drdy = [0 1 0]';
drdz = [0 0 1]';
dr1dtheta = [0 0 0]';
dr2dtheta = [0,... 
            cosd(theta)*(-h*cosd(phi)+l*cosd(psi)*sind(phi))-l*sind(psi)*sind(theta),...
            -l*cosd(theta)*sind(psi)+(h*cosd(phi)-l*cosd(psi)*sind(phi))*sind(theta)]';
dr3dtheta = [0,...
            -cosd(theta)*(h*cosd(phi)+l*cosd(psi)*sind(phi))+l*sind(psi)*sind(theta),...
            l*cosd(theta)*sind(psi)+(h*cosd(phi)+l*cosd(psi)*sind(phi))*sind(theta)]';
dr4dtheta = [0,...
            cosd(theta)*(h*cosd(phi)+l*sind(phi)*sind(psi))+l*cosd(psi)*sind(theta),...
            l*cosd(psi)*cosd(theta)-(h*cosd(phi)+l*sind(phi)*sind(psi))*sind(theta)]';
dr5dtheta = [0,...
            h*cosd(phi)*cosd(theta)-l*(cosd(theta)*sind(phi)*sind(psi)+cosd(psi)*sind(theta)),...
            -l*cosd(psi)*cosd(theta)+(-h*cosd(phi)+l*sind(phi)*sind(psi))*sind(theta)]';
dr1dphi = [0 0 0]';
dr2dphi = [-h*cosd(phi)+l*cosd(psi)*sind(phi),...
          (l*cosd(phi)*cosd(psi)+h*sind(phi))*sind(theta),...
          (l*cosd(phi)*cosd(psi)+h*sind(phi))*cosd(theta)]';
dr3dphi = [-h*cosd(phi)-l*cosd(psi)*sind(phi),...
           -(l*cosd(phi)*cosd(psi)-h*sind(phi))*sind(theta),...
           -(l*cosd(phi)*cosd(psi)-h*sind(phi))*cosd(theta)]';
dr4dphi = [h*cosd(phi)+l*sind(phi)*sind(psi),...
          (-h*sind(phi)+l*cosd(phi)*sind(psi))*sind(theta),...
          (-h*sind(phi)+l*cosd(phi)*sind(psi))*cosd(theta)]';
dr5dphi = [h*cosd(phi)-l*sind(phi)*sind(psi),...
          (-h*sind(phi)-l*cosd(phi)*sind(psi))*sind(theta),...
          (-h*sind(phi)-l*cosd(phi)*sind(psi))*cosd(theta)]';
dr1dpsi = [0 0 0]';
dr2dpsi = [l*cosd(phi)*sind(psi),...
           l*(cosd(psi)*cosd(theta)-sind(phi)*sind(psi)*sind(theta)),...
           -l*(cosd(theta)*sind(phi)*sind(psi)+cosd(psi)*sind(theta))]';
dr3dpsi = [-l*cosd(phi)*sind(psi),...
           -l*cosd(psi)*cosd(theta)+l*sind(phi)*sind(psi)*sind(theta),...
           l*(cosd(theta)*sind(phi)*sind(psi)+cosd(psi)*sind(theta))]';
dr4dpsi = [-l*cosd(phi)*cosd(psi),...
            l*(cosd(theta)*sind(psi)+cosd(psi)*sind(phi)*sind(theta)),...
            l*(cosd(psi)*cosd(theta)*sind(phi)-sind(psi)*sind(theta))]';
dr5dpsi = [l*cosd(phi)*cosd(psi),...
          -l*(cosd(theta)*sind(psi)+cosd(psi)*sind(phi)*sind(theta)),...
          -l*cosd(psi)*cosd(theta)*sind(phi)+l*sind(psi)*sind(theta)]';

F(1) = F1*drdx + F2*drdx + F3*drdx + F4*drdx + F5*drdx;
F(2) = F1*drdy + F2*drdy + F3*drdy + F4*drdy + F5*drdy;
F(3) = F1*drdz + F2*drdz + F3*drdz + F4*drdz + F5*drdz;
F(4) = F1*dr1dtheta + F2*dr2dtheta + F3*dr3dtheta + F4*dr4dtheta + F5*dr5dtheta;
F(5) = F1*dr1dphi + F2*dr2dphi + F3*dr3dphi + F4*dr4dphi + F5*dr5dphi;
F(6) = F1*dr1dpsi + F2*dr2dpsi + F3*dr3dpsi + F4*dr4dpsi + F5*dr5dpsi;

end

