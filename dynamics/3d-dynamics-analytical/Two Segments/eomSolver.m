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

function dqdt = eomSolver(t,yy,node1,R,m,g,cv,cs,kv,ks,Ls0,Lv0)
%INPUT: yy = 12 states
%OUTPUT:
    
    % Force at each node
    if t<10
    F1 = [0 0 0];
    F2 = [0 0 0];
    F3 = [0 0 0];
    F4 = [0 0 0];
    F5 = [0 0 0];
    F0 = zeros(1,3);    % No force
    else
    F1 = [0 0 0];
    F2 = [0 0 0];
    F3 = [0 0 0];
    F4 = [0 0 0];
    F5 = [0 0 0];
    F0 = zeros(1,3);    % No force            
    end
    alpha = 109.5;    % angle b/w rods [degree]
    % Distance on 2D projection
    l = R*sind(alpha/2);
    h = R*cosd(alpha/2);

    x = yy(1); y = yy(2); z = yy(3); theta = yy(4); phi = yy(5); psi = yy(6);
    xd = yy(7); yd = yy(8); zd = yy(9);
    thetad = deg2rad(yy(10)); phid = deg2rad(yy(11)); psid = deg2rad(yy(12)); % Angular Velocity in [rad/s]
    
%     fprintf('%f\n',z);
    
    dqdt = zeros(12,1);
    F = appliedForce(node1,R,m,g,cv,cs,kv,ks,Ls0,Lv0,F1,F2,F3,F4,F5,yy);
    
    dqdt(1:6) = [xd,yd,zd,thetad,phid,psid]';
    dqdt(7:12) = [F(1)/(5*m),F(2)/(5*m),F(3)/(5*m),...
            (secd(phi)*(F(4)*secd(phi)-F(6)*tand(phi)-4*m*phid*((l^2)*psid-2*(h^2)*sind(phi)*thetad)))/...
            (2*m*(2*(h^2)+(l^2))),...
            (F(5)+4*(l^2)*m*cosd(phi)*psid*thetad-(2*(h^2)-(l^2))*m*sind(2*phi)*(thetad^2))/...
            (2*m*(2*(h^2)+(l^2))),...
            ((secd(phi)^2)/(8*m*(l^2)*(2*(h^2)+(l^2))))*...
            (2*F(6)*(h^2)+3*F(6)*(l^2)+F(6)*(2*(h^2)-(l^2))*cosd(2*phi)-4*F(4)*(l^2)*sind(phi)-...
             4*(l^2)*m*cosd(phi)*phid*(-4*(l^2)*sind(phi)*psid+(6*(h^2)+(l^2)+(-2*(h^2)+(l^2))*cos(2*phi))*thetad))]';
end