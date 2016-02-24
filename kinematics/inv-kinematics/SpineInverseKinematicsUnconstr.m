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
w=ones(20,1)*10;
%This is the bar connectivity matrix for one tetrahedron
     %   1  2  3  4  5  6  7  8  9 10 11 12  
B_base=[ 1  1  1  0  0  0;   %A
        -1  0  0  1  1  0;   %B
         0  0 -1 -1  0  1;   %C
         0 -1  0  0 -1 -1]';  %D
       
pretension = 200; %force density N/m

%This is the string connetivity matrix for Any two tetrahedra
     %   1  2  3  4  5  6  7  8   
C_base=[ 1  0  0  0  0  0  0  0 ; %A
         0  1  0  0  0  0  0  0 ; %B
         0  0  1  0  0  1  1  0 ; %C
         0  0  0  1  1  0  0  1 ; %D
        -1  0  0  0 -1 -1  0  0 ; %E
         0 -1  0  0  0  0 -1 -1 ; %F
         0  0 -1  0  0  0  0  0 ; %G
         0  0  0 -1  0  0  0  0]';%H


%Number of tetrahedra (spine nodes)
N=5;
C=zeros((N-1)*8,N*4);
B=zeros(N*6,N*4);

%Connectivity Matrices superimposed for number of tetrahedra
for i=0:(N-2)
    C= C +[zeros(8*i,N*4);
           zeros(8,4*i),  C_base     ,  zeros(8,4*(N-2-i));
           zeros(8*(N-2-i),N*4)];
end
for i=0:N-1
    B= B + [zeros(6*i,N*4);
            zeros(6,4*i), B_base,  zeros(6,4*(N-1-i))
            zeros(6*(N-1-i),N*4)];
end

K_scale=eye((N-1)*8);


Mass= 1; %mass of single tetra kg
L=0.2; %m, Edge of tetra 

%These are coordinates for the first tetrahedron. All subsequent
%tetrahedron are multiplied by the same translations repeatedly... this way
%coupled actuation is straightforward if you kinematically constrain
%yourself to these 4 degrees of freedom
z = 0.075; %vertical height
xR = 0;    %x rotation bending
yR = 0;    %y rotation bending
zR = 0;    %torsional rotation

%Cable stiffnesses, only used for back calculating cable length from force
%densities
K=1000*ones((N-1)*8,1); %N/m
 
hFig = figure(1);
set(hFig, 'Position', [800 50 1000 1200])
% -----------Comment/Uncomment below accordingly-------------
options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');       %2014 ver.
% options = optimset('Algorithm','interior-point-convex','Display','off');                        %2012 ver.
% -----------------------------------------------------------
state=1;
h=sqrt(L^2-1/2*L^2); % height of tetrahedron... I made this equation a while ago I should double check it. 


%Standard Tetrahedron Coord
tetraNodesPreTransform = [L/2   0     -h/2  1; %A
                          -L/2   0     -h/2  1; %B                 
                          0    -L/2    h/2  1; %C
                          0     L/2    h/2  1];%D

holder1=tetraNodesPreTransform;
                           
%for plotting purposes
centersPreTransform=[0 0 0 1];
holder2=centersPreTransform;
     
     
for i=1:N-1
    tetraNodesPreTransform = [tetraNodesPreTransform; holder1];
    centersPreTransform=[centersPreTransform; holder2];
end

%8 strings x 2 [points/string] (4 verticals + 4 saddles)
%X,Y,Z
stringPts=zeros(16*(N-1),3);
temp = zeros((N-1)*8,3);


% iterate over number of frames to render
frame=0;
num_frames_to_render = 40;

% Uncomment these lines to save a video
%videoObject = VideoWriter('videos/SpineExample.avi');
%videoObject.Quality = 100;
%videoObject.FrameRate = 5;


% Keep track of cable lengths over successive iterations. Column is time index.
% Note that I just copied Jeff's sizing for his stringLengths variable here,
% I haven't actually counted and justified to myself that we really have 16*(N-1) cables.
stringLengthsOverTime = zeros(16*(N-1)-1, num_frames_to_render);

% ----------------Record changes in cable length (CHANWOO)---------------------------
% stringLengthHistoryVert = zeros(number of connections, 4verticals , time scale);
stringLengthHistoryVert = zeros(N-1,4,num_frames_to_render);
% stringLengthHistorySadd = zeros(number of connections, 4saddles , time scale);
stringLengthHistorySadd = zeros(N-1,4,num_frames_to_render);
%------------------------------------------------------------------------------------


% Main loop
while frame < num_frames_to_render
    frame=frame+1;      %=counter
    
    % Define the incremental rotations/bending of successive tetra nodes
    % Uncomment zR to see torsion, and yR for bending
    if zR<3.14/8 && yR>-3.14/8 && state==1
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

    for i=1:N-1    
        tetraNodes((4*i+1):end,:)=(zRotation*zTranslation *tetraNodes((4*i+1):end,:)')'; 
        centers((i+1):N,:)=(zRotation*zTranslation*centers((i+1):N,:)')';
    end
    
    %After you translate the tetrahedrons vertically and spin them about
    %the z-axis you next create an appropriate rotation matrix to spin
    %them about the bar of the tetrahedron to maintain the same saddle
    %string positions
    for i=1:N-1 
        xRotTranslate = makehgtform('translate',-tetraNodes(4*i,1:3));
        xRotation = makehgtform('axisrotate',tetraNodes(4*i-1,1:3)- tetraNodes(4*i,1:3),xR);
        xRotTranslateInv = makehgtform('translate',tetraNodes(4*i,1:3));
        tetraNodes((4*i+1):4*N,:)=(xRotTranslateInv*xRotation*xRotTranslate*tetraNodes((4*i+1):4*N,:)')'; centers((i+1):N,:)=(xRotTranslateInv*xRotation*xRotTranslate*centers((i+1):N,:)')';
    end 
    for i=1:N-1
        yRotTranslate = makehgtform('translate',-tetraNodes(4*i+2,1:3));
        yRotation = makehgtform('axisrotate',tetraNodes(4*i+1,1:3)- tetraNodes(4*i+2,1:3),yR);
        yRotTranslateInv = makehgtform('translate',tetraNodes(4*i+2,1:3));
        tetraNodes((4*i+1):4*N,:)=(yRotTranslateInv*yRotation*yRotTranslate*tetraNodes((4*i+1):4*N,:)')'; centers((i+1):N,:)=(yRotTranslateInv*yRotation*yRotTranslate*centers((i+1):N,:)')';  
    end

    for i=0:N-2
        stringPts(((1:16)+16*i),:)=tetraNodes(([1,5,2,6,3,7,4,8,4,5,3,5,3,6,4,6]+i*4),1:3); %string 1
    end
            
    %Newtons Loading of the tetrahedrons, currently configured just to hold up their own mass
    F = [zeros(N*8,1); (N-1)*Mass/4*9.8*ones(4,1); -Mass/4*9.8*ones((N-1)*4,1)]; 
    
    % Set up the optimization that will solve for the inverse kinematics at
    % this timestep
    C_A=[C; B];
    A= [C_A' *diag(C_A*tetraNodes(:,1));
        C_A' *diag(C_A*tetraNodes(:,2));
        C_A' *diag(C_A*tetraNodes(:,3))];
    A_g = pinv( A);
    A_g_A=A_g*A;
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
    
    % Run the actual optimization for the inverse kinematics
    w = quadprog(V(1:(N-1)*8,:)'*K_scale*V(1:(N-1)*8,:),V(1:(N-1)*8,:)'*K_scale*A_g(1:(N-1)*8,:)*F,-V(1:(N-1)*8,:),A_g(1:(N-1)*8,:)*F-pretension,[],[],[],[],[],options);
    q=A_g*F + V*w;
    
    stringLengths=getLengths(stringPts(:,1),stringPts(:,2),stringPts(:,3));
    Lengths= stringLengths(1:2:end);
    L0=Lengths-Lengths.*q(1:(N-1)*8)./K;
    
    % Record the string lengths
    stringLengthsOverTime(:,frame) = stringLengths;
    

% Storing string points for plotting & get length
%     1 3 5 7
%     2 4 6 8
%     18 20 22 24
%     34 36 38 40
%     50 52 54 56
    
    for i = 1:(N-1)
        temp(1+8*(i-1),:) = stringPts(1,:);
        temp(3+8*(i-1),:) = stringPts(3,:);
        temp(5+8*(i-1),:) = stringPts(5,:);
        temp(7+8*(i-1),:) = stringPts(7,:);
        temp(2+8*(i-1),:) = stringPts(2+(i-1)*16,:);
        temp(4+8*(i-1),:) = stringPts(4+(i-1)*16,:);
        temp(6+8*(i-1),:) = stringPts(6+(i-1)*16,:);
        temp(8+8*(i-1),:) = stringPts(8+(i-1)*16,:);
    end
    
    %----------- CABLE LENGTH RECORDER (CHANWOO)-------------------------
    for i = 1:N-1       % Connection between tetrahedrons
        for j = 1:4     % 4verticals & 4saddles per 1 connection
            stringLengthHistoryVert(i,j,frame) = getLengths(temp((j*2-1:j*2)+(i-1)*8,1),temp((j*2-1:j*2)+(i-1)*8,2),temp((j*2-1:j*2)+(i-1)*8,3));

            stringLengthHistorySadd(i,j,frame) = getLengths(stringPts(((j*2-1)+8:j*2+8)+16*(i-1),1),stringPts(((j*2-1)+8:j*2+8)+16*(i-1),2),stringPts(((j*2-1)+8:j*2+8)+16*(i-1),3));
        end
    end
    %---------------------------------------------------------------------

    % Calculate the forces in the cables, and record the color differently
    % if they're in tension or compression
    Force=Lengths.*q(1:(N-1)*8);
    for i=1:(N-1)*8
        if Force(i)>0
            color(i)='b';
        else
            color(i)='r';
        end
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
    
    %axis equal;
    % For 2D perspective:
    %view([0 0]);
    view([90 0]);
    
    % Plot the strings, with the thickness proportional to the force
    
    for i=1:(N-1)*4    
        plot3(temp((1:2)+2*(i-1),1),temp((1:2)+2*(i-1),2),temp((1:2)+2*(i-1),3),color(1+i),'LineWidth',abs(0.1*Force(1+i)))
    end
    for i=1:N-1;
       plot3(stringPts((9:10)+16*(i-1),1),stringPts((9:10)+16*(i-1),2),stringPts((9:10)+16*(i-1),3),color(1+i),'LineWidth',abs(0.1*Force(1+i)))
       plot3(stringPts((11:12)+16*(i-1),1),stringPts((11:12)+16*(i-1),2),stringPts((11:12)+16*(i-1),3),color(1+i),'LineWidth',abs(0.1*Force(1+i)))
       plot3(stringPts((13:14)+16*(i-1),1),stringPts((13:14)+16*(i-1),2),stringPts((13:14)+16*(i-1),3),color(1+i),'LineWidth',abs(0.1*Force(1+i)))
       plot3(stringPts((15:16)+16*(i-1),1),stringPts((15:16)+16*(i-1),2),stringPts((15:16)+16*(i-1),3),color(1+i),'LineWidth',abs(0.1*Force(1+i)))
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
plot(stringLengthsOverTime(24,:));

% Plot changes in length of each cable over time (CHANWOO)
plotLengthChange
gearRatioFinder


% Save the movie we generated
% Uncomment these 3 lines to save a video
%open(videoObject);
%writeVideo(videoObject, M);
%close(videoObject);
    
    