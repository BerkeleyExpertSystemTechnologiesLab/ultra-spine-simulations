% @Author ChanWoo Yang
% UC Berkeley
% BEST Lab
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Created 4/03/2015
% Modified 4/12/2015
% Contact ChanWoo at: chanwoo.yang@berkeley.edu
% Tensegrity Spine Dynamics: Stellated Tetrahedron Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dqdt = eomSolver(t,yy,l,h,m,F1,F2,F3,F4,F5)
%INPUT: yy = 12 states
%OUTPUT:
    x = yy(1); y = yy(2); z = yy(3); theta = yy(4); phi = yy(5); psi = yy(6);
    xd = yy(7); yd = yy(8); zd = yy(9); 
    thetad = deg2rad(yy(10)); phid = deg2rad(yy(11)); psid = deg2rad(yy(12)); % Angular Velocity in [rad/s]
    
    dqdt = zeros(12,1);
    F = appliedForce(l,h,F1,F2,F3,F4,F5,theta,phi,psi);
    
    dqdt(1:6) = [xd,yd,zd,thetad,phid,psid]';
    dqdt(7:12) = [F(1)/(5*m),F(2)/(5*m),F(3)/(5*m),...
                (F(4)-F(6)*sind(phi)+4*m*cosd(phi)*phid*(-(l^2)*psid+2*(h^2)*sind(phi)*thetad))/...
                (2*(h^2)*m+3*(l^2)*m+(2*(h^2)-(l^2))*m*cosd(2*phi)-4*m*(l^2)*(sind(phi)^2)),...
                (F(5)+4*(l^2)*m*cosd(phi)*psid*thetad-(2*(h^2)-(l^2))*m*sind(2*phi)*(thetad^2))/(2*m*(2*(h^2)+(l^2))),...
                (((-m*(2*(h^2)+3*(l^2)+(2*(h^2)-(l^2))*cosd(2*phi))*(F(6)-4*(l^2)*m*cosd(phi)*phid*thetad))/(m*(l^2)))+...
                (4*sind(phi)*(F(4)-4*m*cosd(phi)*phid*((l^2)*psid+(-2*(h^2)+(l^2))*sind(phi)*thetad))))/...
                (-4*m*(2*(h^2)+3*(l^2)+(2*(h^2)-(l^2))*cosd(2*phi)-4*(l^2)*((sind(phi))^2)))]';
end