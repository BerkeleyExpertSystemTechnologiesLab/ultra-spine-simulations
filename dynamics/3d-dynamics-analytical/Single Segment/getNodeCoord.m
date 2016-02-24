% @Author ChanWoo Yang
% UC Berkeley
% BEST Lab
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Created 4/03/2015
% Modified 4/04/2015
% Contact ChanWoo at: chanwoo.yang@berkeley.edu
% Tensegrity Spine Dynamics: Stellated Tetrahedron Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = getNodeCoord(T,Y)

r = zeros(3,5,length(T));
% r([x,y,z],node#,timeStep)

alpha = 109.5;    %angle b/w rods

%Distance from the center of a regular tetrahedron to a node
%Length of Rod
R = 1;

% Distance on 2D projection
l = R*sind(alpha/2);
h = R*cosd(alpha/2);

t2 = [-l 0 -h]';
t3 = [l 0 -h]';
t4 = [0 -l h]';
t5 = [0 l h]';

    for i = 1:length(T)
        r(:,1,i) = Y(i,1:3)';
        theta = Y(i,4); %about x-axis
        phi = Y(i,5);   %about y-axis
        psi = Y(i,6);  %about z-axis

        % Rotational Matrix
        Tx = [1 0 0;
              0 cosd(theta) sind(theta);
              0 -sind(theta) cosd(theta)];
        Ty = [cosd(phi) 0 sind(phi);
              0 1 0;
              -sind(phi) 0 cosd(phi)];
        Tz = [cosd(psi) sind(psi) 0;
              -sind(psi) cosd(psi) 0;
              0 0 1];

        % Nodal Position after Rotation
        e2 = Tx*Ty*Tz*t2;
        e3 = Tx*Ty*Tz*t3;
        e4 = Tx*Ty*Tz*t4;
        e5 = Tx*Ty*Tz*t5;

        r(:,2,i) = r(:,1,i) + e2;
        r(:,3,i) = r(:,1,i) + e3;
        r(:,4,i) = r(:,1,i) + e4;
        r(:,5,i) = r(:,1,i) + e5;
    end
end
