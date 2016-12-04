% @Author Jeff Friesen, minor edits and cleanup by Drew Sabelhaus
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Modified 8/01/2014
% Contact Jeff at: jfriesen222@gmail.com
% Contact Drew at: apsabelhaus@berkeley.edu

% Instructions:
% - before starting, be sure that the getLengths function is in your MATLAB path: 
%   an easy way to do this is to make sure that getLengths.m is in the same folder as this file
% - if you get errors about missing functions, such as quadprog (for the optimization) you might be missing a required MATLAB toolbox.
% 1) set N on line 41, the number of spine nodes
% 2) set num_frames_to_render on line 101. this is "how long the movie will be."
% 3) adjust the pause at the end of this file, if things aren't rendering quickly or properly
% 4) adjust the axis limits also at the end of this file to scale thingsdifferently if so desired

clc;
clear all;
close all;


%w is a vector used to navigate the nullspace of cable force densities
w = ones(20,1)*10;
%This is the bar connectivity matrix for one tetrahedron
%          1  2  3  4  5  6
B_base = [ 1  1  1  0  0  0;   %A
          -1  0  0  1  1  0;   %B
           0  0 -1 -1  0  1;   %C
           0 -1  0  0 -1 -1]'; %D
% DREW: Why tetrahedron edges instead of bars?
     
pretension = 200; %force density N/m

%This is the string connetivity matrix for Any two tetrahedra
%          1  2  3  4  5  6  7  8   
C_base = [ 1  0  0  0  0  0  0  0 ;  %A
           0  1  0  0  0  0  0  0 ;  %B
           0  0  1  0  0  1  1  0 ;  %C
           0  0  0  1  1  0  0  1 ;  %D
          -1  0  0  0 -1 -1  0  0 ;  %E
           0 -1  0  0  0  0 -1 -1 ;  %F
           0  0 -1  0  0  0  0  0 ;  %G
           0  0  0 -1  0  0  0  0]'; %H
% Note B_base and C_base are transposed to match Schek's definitions

%Number of tetrahedra (spine nodes)
N = 5;
% B = zeros(N*6,N*4);
% C = zeros((N-1)*8,N*4);
B = zeros(N*size(B_base,1),N*size(B_base,2));
C = zeros((N-1)*size(C_base,1),N*size(B_base,2));
% B and C need to be the same width

%Connectivity Matrices superimposed for number of tetrahedra
for i=0:N-1
    B = B + [zeros(6*i,N*4);
             zeros(6,4*i), B_base,  zeros(6,4*(N-1-i))
             zeros(6*(N-1-i),N*4)];
end
for i=0:N-2
    C = C +[zeros(8*i,N*4);
            zeros(8,4*i),  C_base,  zeros(8,4*(N-2-i));
            zeros(8*(N-2-i),N*4)];
end

K_scale=eye((N-1)*8);
% DREW: Why do you have this scaling factor?

Mass = 1; %kg, mass of single tetra
L = 0.2; %m, edge length of tetra 

%These are coordinates for the first tetrahedron. All subsequent
%tetrahedron are multiplied by the same translations repeatedly... this way
%coupled actuation is straightforward if you kinematically constrain
%yourself to these 4 degrees of freedom
% z = 0.075; %vertical height
xR = 0;    %x rotation bending
z = 0.75; %vertical height
yR = 0;    %y rotation bending
zR = 0;    %torsional rotation

%Cable stiffnesses, only used for back calculating cable length from force
%densities
K = 1000*ones((N-1)*8,1); %N/m

hFig = figure(1);
% set(hFig, 'Position', [800 50 1000 1200])
% -----------Comment/Uncomment below accordingly-------------
options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off'); %2014 ver.
% DREW: Why these optimization options?

% options = optimset('Algorithm','interior-point-convex','Display','off'); %2012 ver.
% -----------------------------------------------------------
state = 1; % DREW: What does this boolean do?
h = sqrt(L^2-1/2*L^2); % height from origin to top of tetra (DOUBLE CHECK)
% h = L*sqrt(0.5); % equivalent, clearer version of height

%Standard Tetrahedron Coord
tetraNodesPreTransform = [ L/2     0  -h/2  1;  %A
                          -L/2     0  -h/2  1;  %B                 
                             0  -L/2   h/2  1;  %C
                             0   L/2   h/2  1]; %D                      

holder1 = tetraNodesPreTransform;

%for plotting purposes
centersPreTransform = [0 0 0 1];
holder2 = centersPreTransform;
     
     
for i = 1:N-1
    tetraNodesPreTransform = [tetraNodesPreTransform; holder1];
    centersPreTransform = [centersPreTransform; holder2];
end
% Right now I have five tetrahedrons stacked on top of each other

% 8 cables x 2 [points/string] (4 verticals + 4 saddles)
% X,Y,Z
cablePts = zeros(16*(N-1),3);

% Iterate over number of frames to render
frame = 0;
num_frames_to_render = 20;

% Uncomment these lines to save a video
%videoObject = VideoWriter('videos/SpineExample.avi');
%videoObject.Quality = 100;
%videoObject.FrameRate = 5;

% Keep track of cable lengths over successive iterations. Column is time index.
% Note that I just copied Jeff's sizing for his stringLengths variable here,
% I haven't actually counted and justified to myself that we really have 16*(N-1) cables.
cableLengthsOverTime = zeros(16*(N-1)-1, num_frames_to_render);
% NOTE: This size doesn't make sense.

% Maybe we just need every other element
lengthsOverTime = zeros((N-1)*8, num_frames_to_render);
% NOTE: This size makes sense.

% ----------------Record changes in cable length (CHANWOO)---------------------------
% stringLengthHistoryVert = zeros(number of connections, 4verticals , time scale);
stringLengthHistoryVert = zeros(N-1,4,num_frames_to_render);
% stringLengthHistorySadd = zeros(number of connections, 4saddles , time scale);
stringLengthHistorySadd = zeros(N-1,4,num_frames_to_render);
%------------------------------------------------------------------------------------

% Record changes in cable force (Drew) 1-23-15
% hard-coded constant comes from the observation that there are 8 cables connecting two
% adjacent tetrahedra
stringForceHistory = zeros((N-1)*8,num_frames_to_render);

% Record the actual position of the robot (Drew) 2016-04-22
% Each 2D 'center' array is for 5 tetrahedra, with 4 dimensions (x,y,z, 1).
centersHistory = zeros(N, 4, num_frames_to_render);
% The angles of rotation for each tetra at each timepoint will be recorded
% There are 5 tetras.
rotationHistory = zeros(5, num_frames_to_render);

% Main loop
while frame < num_frames_to_render
    frame=frame+1; % =counter
    
    % Define the incremental rotations/bending of successive tetra nodes
    % Uncomment zR to see torsion, and yR for bending
    if zR<3.14/8 && yR>-3.14/8 && state==1 % DREW: How did you determine these bounds?
%         zR = zR + 0.02;
        yR = yR - 0.01;
    else   
        if zR>-3.14/8 && yR<3.14/8 && state==-1
%             zR = zR - 0.02;
            yR = yR + 0.01;
        else
            state=-state;
        end
    end

    %Reset the nodes for the next translation... it's easier this way
    tetraNodes=tetraNodesPreTransform;
    centers=centersPreTransform;
    %These two transforms aren't dependent upon coordinates as I've defined
    %them so leave outside of the for loop
    zTranslation = makehgtform('translate',[0 0 z]);
    zRotation = makehgtform('zrotate',zR);
    % Help explaining 4x4 transformation matrices: http://www.euclideanspace.com/maths/geometry/affine/matrix4x4/

    for i=1:N-1    
        tetraNodes((4*i+1):end,:)=(zRotation*zTranslation*tetraNodes((4*i+1):end,:)')'; 
        centers((i+1):N,:)=(zRotation*zTranslation*centers((i+1):N,:)')';
    end
    % I now have a tower, possible spinning like DNA if zR was used
    
    %After you translate the tetrahedrons vertically and spin them about
    %the z-axis you next create an appropriate rotation matrix to spin
    %them about the bar of the tetrahedron to maintain the same saddle
    %string positions
    for i=1:N-1 
        xRotTranslate = makehgtform('translate',-tetraNodes(4*i,1:3));
        xRotation = makehgtform('axisrotate',tetraNodes(4*i-1,1:3)- tetraNodes(4*i,1:3),xR);
        xRotTranslateInv = makehgtform('translate',tetraNodes(4*i,1:3));
        tetraNodes((4*i+1):4*N,:)=(xRotTranslateInv*xRotation*xRotTranslate*tetraNodes((4*i+1):4*N,:)')';
        centers((i+1):N,:)=(xRotTranslateInv*xRotation*xRotTranslate*centers((i+1):N,:)')';
    % NOTE: I'm translating it down to the origin, rotating it, and then
    % translating it back up
    end 
    for i=1:N-1
        yRotTranslate = makehgtform('translate',-tetraNodes(4*i+2,1:3));
        yRotation = makehgtform('axisrotate',tetraNodes(4*i+1,1:3)- tetraNodes(4*i+2,1:3),yR);
        yRotTranslateInv = makehgtform('translate',tetraNodes(4*i+2,1:3));
        tetraNodes((4*i+1):4*N,:)=(yRotTranslateInv*yRotation*yRotTranslate*tetraNodes((4*i+1):4*N,:)')';
        centers((i+1):N,:)=(yRotTranslateInv*yRotation*yRotTranslate*centers((i+1):N,:)')';  
    end

    for i=0:N-2
        cablePts(((1:16)+16*i),:)=tetraNodes(([1,5,2,6,3,7,4,8,4,5,3,5,3,6,4,6]+i*4),1:3); %cable 1
    end
            
    %Newtons Loading of the tetrahedrons, currently configured just to hold up their own mass
    F = [zeros(N*8,1); (N-1)*Mass/4*9.8*ones(4,1); -Mass/4*9.8*ones((N-1)*4,1)]; 
    % DREW: How does this work?
    
    % Set up the optimization that will solve for the inverse kinematics at
    % this timestep
    C_A=[C; B];
    A= [C_A' *diag(C_A*tetraNodes(:,1));
        C_A' *diag(C_A*tetraNodes(:,2));
        C_A' *diag(C_A*tetraNodes(:,3))];
    A_g = pinv(A);
    A_g_A = A_g*A;
    V=(eye(length(A_g_A))-A_g_A);
    [Q,R,E] = qr(V);
    [m , n] = size(R);
    j=1;
    i=1;
    while i<=m
        if norm(R(i,:))>10^-12
            R_new(j,:)=R(i,:);
            j=j+1;
        else
            i=m;
        end
        i=i+1;
    end
    V=Q(:,1:j-1);
    % DREW: Why do you do this decomposition and change V?
    
    % Run the actual optimization for the inverse kinematics
    w = quadprog(V(1:(N-1)*8,:)'*K_scale*V(1:(N-1)*8,:), V(1:(N-1)*8,:)'*K_scale*A_g(1:(N-1)*8,:)*F, -V(1:(N-1)*8,:), A_g(1:(N-1)*8,:)*F-pretension,[],[],[],[],[],options);
    q=A_g*F + V*w;
    
    stringLengths=getLengths(cablePts(:,1),cablePts(:,2),cablePts(:,3));
    Lengths= stringLengths(1:2:end);
    L0=Lengths-Lengths.*q(1:(N-1)*8)./K;
    
    % Record the string lengths
    cableLengthsOverTime(:,frame) = stringLengths;
    lengthsOverTime(:,frame) = Lengths;
    
    % Calculate the forces in the cables, and record the color differently
    % if they're in tension or compression
    % only search the first 's' number of elements because those are the cables
    % this is for example 4*8 = 32
    Force=Lengths.*q(1:(N-1)*8);
    for i=1:(N-1)*8
        if Force(i)>0
            color(i)='b';
        else
            color(i)='r';
        end
    end
    
    %----------- CABLE LENGTH RECORDER (CHANWOO)-------------------------
    for i = 1:N-1       % Connection between tetrahedrons
        for j = 1:4     % 4verticals & 4saddles per 1 connection
            stringLengthHistoryVert(i,j,frame) = getLengths(cablePts((j*2-1:j*2)+16*(i-1),1),cablePts((j*2-1:j*2)+16*(i-1),2),cablePts((j*2-1:j*2)+16*(i-1),3));

            stringLengthHistorySadd(i,j,frame) = getLengths(cablePts(((j*2-1)+8:j*2+8)+16*(i-1),1),cablePts(((j*2-1)+8:j*2+8)+16*(i-1),2),cablePts(((j*2-1)+8:j*2+8)+16*(i-1),3));
        end
    end
    %---------------------------------------------------------------------
    
    % Cable Force Recorder (Drew)
    stringForceHistory(:,frame) = Force;  
    
    % Tetrahedra center nodes history (Drew)
    centersHistory(:,:, frame) = centers;
        
    % Rotation for each tetrahedron can be calculated by the slope of the line connecting
    % its top two nodes. These are the 3rd and 4th points for each tetra.
    for i=1:5
        % Note that we're ignoring the X dimension here. TO-DO: generalize this. (2016-04-22).
        % THIS IS NOT VALID FOR 3D ROTATIONS.
        % Z is the third element, and we want to calculate e.g. rows 4-3, rows 8-7, 12-11, 16-15, 20-19.
        % Orientation here is in
        tetraZdist = tetraNodes( i*4, 3) - tetraNodes( i*4 - 1, 3);
        tetraYdist = tetraNodes( i*4, 2) - tetraNodes( i*4 - 1, 2);
        tetraRotation = atan(tetraZdist/tetraYdist);
        % save this angle (in radians) for the i-th tetrahedron, for this timestep.
        rotationHistory(i, frame) = tetraRotation;
    end
    
    hold off;

    % Plot the tetrahedron structures
    for i=0:N-1
        tetraPlot = [tetraNodes(1+i*4,1:3);
                     centers(1+i,1:3);
                     tetraNodes(2+i*4,1:3);
                     centers(1+i,1:3);
                     tetraNodes(3+i*4,1:3);
                     centers(1+i,1:3);
                     tetraNodes(4+i*4,1:3);
                     centers(1+i,1:3)];
        plot3(tetraPlot(:,1),tetraPlot(:,2),tetraPlot(:,3),'k','LineWidth',3);
        hold on;
    end
    grid on;
    
    % Save tetrahedra locations (Drew)
    %tetraPlotHistory
    
    %axis equal;
    % For 2D perspective:
    %view([0 0]);
    view([90 0]);
    
    % Plot the strings, with the thickness proportional to the force
    for i=0:((N-1)*8-1);
       plot3(cablePts((1:2)+2*i,1),cablePts((1:2)+2*i,2),cablePts((1:2)+2*i,3),color(1+i),'LineWidth',abs(0.1*Force(1+i)))
    end
    
    % Rescale this plot's axes
    xlim([-0.3 0.3]);
    ylim([-0.3 0.3]);
    zlim([0 0.4]);
    
    % Rescale the plot's absolute size
    %set(gcf, 'Units', 'pixels');
    %set(gcf, 'Position', [0, 0, 1600, 800]);
    
    % Save this frame
    M(frame) = getframe(gcf);
    
    % Pause in between successive plots - this makes it video-like
    % To slow down the video, set a longer pause here. Units are seconds.
    pause(0.04);
      
end

% Plot the length trajectory of one cable
figure;
hold on;
plot(cableLengthsOverTime(24,:));

% Plot changes in length of each cable over time (CHANWOO)
% Refer to plotLengthChange.m file; gearRatioFinder.m

%plotLengthChange
%gearRatioFinder



% Save the movie we generated
% Uncomment these 3 lines to save a video
%open(videoObject);
%writeVideo(videoObject, M);
%close(videoObject);
    
    