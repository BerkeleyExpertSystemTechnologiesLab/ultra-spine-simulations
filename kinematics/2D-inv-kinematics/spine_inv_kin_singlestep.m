% spine_inv_kin_singlestep.m
% @Author Jeff Friesen, minor edits and cleanup by Drew Sabelhaus
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Modified 8/01/2014
% Contact Jeff at: jfriesen222@gmail.com
% Contact Drew at: apsabelhaus@berkeley.edu

function [node_distances, final_lengths] = spine_inv_kin_singlestep(tetra_vertical_spacing, spine_geometric_parameters)

% This function calculates the required cable lengths in the spine robot for a given kinematic position: it does a single step of inverse kinematics.
% Inputs:
% states = the positions of the spine. This should be a 12 * num_vertebrae vector.
% spine_geometric_parameters = a struct containing the geometry of the spine vertebrae. See, for example, ultra_spine_lqr for its creation and use.
% Outputs:
% lengths = the lengths that would keep this spine in equilibrium

% NOTE THAT AS OF 2016-04-18, THIS SCRIPT ONLY TAKES IN VERTICAL SPINES.
% We're using this to initialize the ultra-spine-lqr script, which will start vertical.
% TO-DO: actually implement rotation of the tetras.

% Instructions:
% - before starting, be sure that the getLengths function is in your MATLAB path: 
%   an easy way to do this is to make sure that getLengths.m is in the same folder as this file
% - if you get errors about missing functions, such as quadprog (for the optimization) you might be missing a required MATLAB toolbox.
% 1) set N on line 41, the number of spine nodes

%% Define constants and other variables to be used later
% Unroll the spine_geometric_parameters struct into individual variables. See the dynamics generation script for more information.
% Gravitational force
g = spine_geometric_parameters.g;
% Total number of spine vertebrae. As of 2016-04-18 for use in ultra_spine_lqr, this is 4.
num_vertebrae = spine_geometric_parameters.N;
% Length of one "leg" of the tetrahedron (straight-line distance from center to outer point)
l = spine_geometric_parameters.l;
% Height of one tetrahedtron
h = spine_geometric_parameters.h;
% total mass of one whole tetrahedron
m_t = spine_geometric_parameters.m_t;



% Hard code the spring constant from the dynamics scripts. TO-DO: roll this into another struct or something.
% Note that this is from getTensions.m in the dynamics folder.
%K = 2000; % N/m

%w is a vector used to navigate the nullspace of cable force densities
w=ones(20,1)*10;

%This is the bar connectivity matrix for one tetrahedron
     %   1  2  3  4  5  6  7  8  9 10 11 12  
B_base=[ 1  1  1  0  0  0;   %A
        -1  0  0  1  1  0;   %B
         0  0 -1 -1  0  1;   %C
         0 -1  0  0 -1 -1]';  %D
       
% TO-DO: Fully include this variable in the control scripts.
% At the moment, is is an independent variable: "pretension" will proportionally change the cable lengths output from this script,
% BUT, changing pretension won't change the fact that the spine will be in equilibrium... maybe?
% Test out: will we have the same behavior (equilibrium), but different cable forces, for different pretension?
%pretension = 200; %force density N/m

% On 2016-04-18, ultra-spine-lqr had cable rest lengths of 0.1 for vertical and 0.187 for horizontal.
% That corresponds to a pretension of zero: the vertical cables are exactly as long as the z_translation between tetrahedra.
pretension = 400;

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
%N=5;
N = num_vertebrae;

% C is the connectivity matrix, as is B.
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

% This is that scaling parameter that Jeff included for the optimization below.
% Set it to 1 here and forget it.
K_scale=eye((N-1)*8);

%% Calculate the physical parameters of each tetrahedron
%Mass= 1; %mass of single tetra kg
Mass = m_t;


tetraNodesPreTransform = [(l^2 - (h/2)^2)^.5, 0, -h/2; ...
             -(l^2 - (h/2)^2)^.5, 0, -h/2; ...
             0, (l^2 - (h/2)^2)^.5, h/2; ...
             0, -(l^2 - (h/2)^2)^.5, h/2];
         
% Make augment this with a 1 vector for the rotations.
tetraNodesPreTransform = [tetraNodesPreTransform, ones(4,1)]; 

%state=1;

% Store a copy of tetraNodesPreTransform to use in copying later.
single_tetra_nodes_atzero=tetraNodesPreTransform;                           


%Cable stiffnesses, only used for back calculating cable length from force densities
%K=1000*ones((N-1)*8,1); %N/m
K=2000*ones((N-1)*8,1); %N/m
 
% The options for the optimization solver
options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');     
     
for i=1:N-1
    tetraNodesPreTransform = [tetraNodesPreTransform; single_tetra_nodes_atzero];
end

%8 strings x 2 [points/string] (4 verticals + 4 saddles)
%X,Y,Z
stringPts=zeros(16*(N-1),3);        

% Loop used to start here
% There was no index for the loop: it was just while frame < frame_max or something like that.

% AS OF 2016-04-18: just move the tetras by tetra_vertical_spacing.
zTranslation = makehgtform('translate',[0 0 tetra_vertical_spacing]);

tetraNodes = tetraNodesPreTransform;

for i=1:N-1    
    tetraNodes((4*i+1):end,:)=(zTranslation*tetraNodes((4*i+1):end,:)')'; 
end



%% Perform the actual optimization.
% This section of the code takes the following parameters to do the calculation:
% tetraNodes
% Mass
% C
% B

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
A_g = pinv( A); % Pseudoinverse
A_g_A=A_g*A;
V=(eye(length(A_g_A))-A_g_A);
[Q,R,E] = qr(V); % orthogonal-triangular decomposition.
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

% Some deciphering by Drew 2016-04-18. K_scale is unused (set to I above.)
% So, this looks like V'V as the first term, and V'A_g as the second term.
% So, the corresponding problem solved by quadprog is min over w of
% 1/2 w' V'V w + (V' A_g F)' w, subject to -V w <= A_g F - pretension.
% Correlating to the ASME IDETC 2015 paper: V is same, w is same, A_g is A^+_s, F is p.
% But, Jeff didn't include pretension in the paper's explanation. The constraint is there to prevent slack cables,
% so I suppose that pretension is just an adjustment to A_g F that truly means "pretension", BUT we need to check units here:
% Since A_g F has units of force density (see the calculation of q below), that means "pretension" is also force density: N/m.
% This is NOT the spring constant.
w = quadprog(V(1:(N-1)*8,:)'*K_scale*V(1:(N-1)*8,:), ...
    V(1:(N-1)*8,:)'*K_scale*A_g(1:(N-1)*8,:)*F, ...
    -V(1:(N-1)*8,:), ...
    A_g(1:(N-1)*8,:)*F-pretension, ...
    [],[],[],[],[],options);

% Given the w from optimization, we can calculate the force density matrix q.
% Jeff's equation in the paper is, equivalently, q_s = A^+_s p + (something) w. He talks about how (something) ends up being V.
q=A_g*F + V*w;

stringLengths=getLengths(stringPts(:,1),stringPts(:,2),stringPts(:,3));
Lengths= stringLengths(1:2:end);
L0=Lengths-(Lengths.*q(1:(N-1)*8))./K;
%L0 = Lengths.*q(1:(N-1)*8)./K;


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

node_distances = Lengths;
final_lengths = Lengths.*q(1:(N-1)*8)./K;

end
    