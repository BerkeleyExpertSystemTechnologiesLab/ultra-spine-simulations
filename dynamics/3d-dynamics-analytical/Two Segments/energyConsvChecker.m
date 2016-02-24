% @Author ChanWoo Yang
% UC Berkeley
% BEST Lab
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Created 4/10/2015
% Modified 4/12/2015
% Contact ChanWoo at: chanwoo.yang@berkeley.edu
% Tensegrity Spine Dynamics: Two Stellated Tetrahedron Segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [totalE,ke,pe] = energyConsvChecker(node1,m,l,h,kv,ks,Ls0,Lv0,T,Y)
g = 9.81;   % Gravitational acc.
x = Y(:,1); y = Y(:,2); z = Y(:,3); theta = Y(:,4); phi = Y(:,5); psi = Y(:,6);
xd = Y(:,7); yd = Y(:,8); zd = Y(:,9); 
thetad = deg2rad(Y(:,10)); phid = deg2rad(Y(:,11)); psid = deg2rad(Y(:,12)); % Angular Velocity in [rad/s]

node2 = getNodeCoord(T,Y);      % node([x,y,z],node#,timeStep)

for i = 1:length(T)
    ke(i) = (1/2)*m*(2*(2*(h^2)+(l^2))*(phid(i)^2)+4*(l^2)*(psid(i)^2)+8*(l^2)*sind(phi(i))*psid(i)*thetad(i)+... 
            2*(h^2)*(thetad(i)^2)+3*(l^2)*(thetad(i)^2)+2*(h^2)*cosd(2*phi(i))*(thetad(i)^2)-...
            (l^2)*cosd(2*phi(i))*(thetad(i)^2)+5*(xd(i)^2)+5*(yd(i)^2)+5*(zd(i)^2));
    peHeight = m*g*(node2(3,1,i)+node2(3,2,i)+node2(3,3,i)+node2(3,4,i)+node2(3,5,i));
    peSpring = (1/2)*(kv*((signChecker(norm(node1(2,:)'-squeeze(node2(:,2,i)))-Lv0))^2)+...
                      ks*((signChecker(norm(node1(4,:)'-squeeze(node2(:,2,i)))-Ls0))^2)+...
                      ks*((signChecker(norm(node1(5,:)'-squeeze(node2(:,2,i)))-Ls0))^2)+...
                      kv*((signChecker(norm(node1(3,:)'-squeeze(node2(:,3,i)))-Lv0))^2)+...
                      ks*((signChecker(norm(node1(4,:)'-squeeze(node2(:,3,i)))-Ls0))^2)+...
                      ks*((signChecker(norm(node1(5,:)'-squeeze(node2(:,3,i)))-Ls0))^2)+...
                      kv*((signChecker(norm(node1(4,:)'-squeeze(node2(:,4,i)))-Lv0))^2)+...
                      kv*((signChecker(norm(node1(5,:)'-squeeze(node2(:,5,i)))-Lv0))^2));
    pe(i) = peHeight+peSpring;

    totalE(i) = ke(i)+pe(i);
end
    
end