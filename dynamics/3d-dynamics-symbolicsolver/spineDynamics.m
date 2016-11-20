% Jeff Friesen, Abishek Akella, Andrew P. Sabelhaus
% Copyright 2016
% This script generates multiple m-files:
%   duct_accel
%   dlengths_dt
%   lengths


clc
clear variables
close all

%% Geometric and physical parameters of the spine robot

% There are 5 total point masses in this model: one at the center, and one
% at each of the 4 outer points at the end of the "legs."
% The geometry of the spine is defined by the bounding box for the 4 outer
% points of the tetrahedron. The center of the tetra is located at (0,0,0).
% So, for example, the z location of each of the outer 4 spine nodes will
% be at (height/2).
% All units are in SI.
% NOTE that this script is not currently parameterized by N. Future
% work will need to change all of the hardcoded N=4 usages below.
% NOTE that since the output functions from this script do not involve the
% springs or the ways the forces are applied between tetra, those constants
% are not included here. (e.g., see other scripts for k, and damping.)
% Parameters are, in order:
%   g, gravitational constant used in the dynamics. (m^3 / (kg s^2))
%   N, number of spine links. (unitless)
%   l, length of spine "leg". The straight-line dist from center to node.
%       in (m).
%   h, total height of one spine tetrahedron. (m) 
%   m_t, total mass of one spine tetrahedron. (kg)
%   FoS, factor-of-safety in the mass of the model. Because we use this
%       model for control, amplifying the mass of the tetrahedron nodes
%       lends some sense of robustness (informally.) (unitless)
%   m, mass of one node of the tetrahedron. = (m_t / 5)*FoS, because there
%       are 5 point masses in this model. Includes FoS too. (kg).

g = 9.81;
N = 4;
l = 0.15;
h = 0.15;
m_t = 0.142;
FoS = 1.2;
m = (m_t/5) * FoS;

% Store all these parameters as a struct for later use.
spine_geometric_parameters.g = g;
spine_geometric_parameters.N = N;
spine_geometric_parameters.l = l;
spine_geometric_parameters.h = h;
spine_geometric_parameters.m_t = m_t;
spine_geometric_parameters.FoS = FoS;
spine_geometric_parameters.m = m;
% The path where we want to save these parameters as a .mat file:
spine_geometric_parameters_path = 'spine_geometric_parameters';

% Creating Parallel Pools for Computing
pools = gcp;

%% Building symbolic matrices for computation
% X, Y, Z, Theta, Gamma, Phi of tetrahedrons for N links
x = sym('x',[N 1]);
y = sym('y',[N 1]);
z = sym('z',[N 1]);
T = sym('T',[N 1]);
G = sym('G',[N 1]);
P = sym('P',[N 1]);

% Velocities for each of above variables
dx = sym('dx',[N 1]);
dy = sym('dy',[N 1]);
dz = sym('dz',[N 1]);
dT = sym('dT',[N 1]);
dG = sym('dG',[N 1]);
dP = sym('dP',[N 1]);

% Accelerations for each of above variables
d2x = sym('d2x',[N 1]);
d2y = sym('d2y',[N 1]);
d2z = sym('d2z',[N 1]);
d2T = sym('d2T',[N 1]);
d2G = sym('d2G',[N 1]);
d2P = sym('d2P',[N 1]);

% Solved accelerations (EOM)
D2x = sym('D2x',[N 1]);
D2y = sym('D2y',[N 1]);
D2z = sym('D2z',[N 1]);
D2T = sym('D2T',[N 1]);
D2G = sym('D2G',[N 1]);
D2P = sym('D2P',[N 1]);

% Lagrangian derivatives for each variable
fx = sym('fx',[N 1]);
fy = sym('fy',[N 1]);
fz = sym('fz',[N 1]);
fT = sym('fT',[N 1]);
fG = sym('fG',[N 1]);
fP = sym('fP',[N 1]);

% Tensions of Cables
tensions = sym('tensions',[8 N]);
lengths = sym('lengths', [8 N]);
dlengths_dt = sym('dlengths_dt', [8 N]);
dlengths_dt(:, 1) = zeros(8, 1);

% Nodal positions of each link from fixed frame
t1 = [0, 0, 0]';
t2 = [(l^2 - (h/2)^2)^.5, 0, -h/2]';
t3 = [-(l^2 - (h/2)^2)^.5, 0, -h/2]';
t4 = [0, (l^2 - (h/2)^2)^.5, h/2]';
t5 = [0, -(l^2 - (h/2)^2)^.5, h/2]';

% Nodal positons of each vertex (N by 3)
r1 = sym('r1',[3 N]);
r2 = sym('r2',[3 N]);
r3 = sym('r3',[3 N]);
r4 = sym('r4',[3 N]);
r5 = sym('r5',[3 N]);

% Nodal velocities of each node for each spine link
dr1 = sym('dr1',[3 N]);
dr2 = sym('dr2',[3 N]);
dr3 = sym('dr3',[3 N]);
dr4 = sym('dr4',[3 N]);
dr5 = sym('dr5',[3 N]);

% First link is fixed in place
r1(:,1) = [0; 0; 0]; % Base link centered around (0,0,0)
r2(:,1) = r1(:,1) + t2; % No rotations
r3(:,1) = r1(:,1) + t3; 
r4(:,1) = r1(:,1) + t4;
r5(:,1) = r1(:,1) + t5;

dr1(:,1) = 0; % Nodal velocities are zero for fixed link
dr2(:,1) = 0;
dr3(:,1) = 0;
dr4(:,1) = 0;
dr5(:,1) = 0;
disp('1%')

% Kinetic energy
T_L = 0;
% Potential energy
V_L=0;

% Following loop to calculate kinetic and potential energies of each link
assume(x,'real')
assume(y,'real')
assume(z,'real')
assume(T,'real')
assume(G,'real')
assume(P,'real')
assume(dx,'real')
assume(dy,'real')
assume(dz,'real')
assume(dT,'real')
assume(dG,'real')
assume(dP,'real')

% Anchor points for first link - these don't matter they just simplify
% indexing
anch_1_bot{1} = [0; 0; 0];
anch_1_top{1} = [0; 0; 0];
anch_2_bot{1} = [0; 0; 0];
anch_2_top{1} = [0; 0; 0];
anch_3_bot{1} = [0; 0; 0];
anch_3_top{1} = [0; 0; 0];
anch_4_bot{1} = [0; 0; 0];
anch_4_top{1} = [0; 0; 0];
angle = pi/6;
anch_5_bot{1} = [0; 0; 0];  
anch_5_top{1} = [0; 0; 0];
anch_6_bot{1} = [0; 0; 0];
anch_6_top{1} = [0; 0; 0];
anch_7_bot{1} = [0; 0; 0];
anch_7_top{1} = [0; 0; 0];
anch_8_bot{1} = [0; 0; 0];
anch_8_top{1} = [0; 0; 0];

prevRot = eye(3);

Lagr = sym('Lagr',[1 N]);
for i=2:N
    % Rotation matrices for theta gamma and phi (euler rotations)

    R1 = [ 1          0           0         ;
           0          cos(T(i))   sin(T(i)) ;
           0         -sin(T(i))   cos(T(i))];

    R2 = [cos(G(i))   0           sin(G(i)) ;
          0           1           0         ;
         -sin(G(i))   0           cos(G(i))];

    R3 = [cos(P(i))   sin(P(i))   0         ;
         -sin(P(i))   cos(P(i))   0         ;
          0           0           1        ];
      disp('2%')

    % Building full Rotation Matrix
    Rot=R3*R2*R1;

    t2 = [(l^2 - (h/2)^2)^.5, 0, -h/2]';
    t3 = [-(l^2 - (h/2)^2)^.5, 0, -h/2]';
    t4 = [0, (l^2 - (h/2)^2)^.5, h/2]';
    t5 = [0, -(l^2 - (h/2)^2)^.5, h/2]';


    % Compute nodal positions after a given rotation
    e2 = Rot*t2;
    e3 = Rot*t3;
    e4 = Rot*t4;
    e5 = Rot*t5;

    disp('3%')
    % Nodal Positions
    r1(:,i) = [x(i); y(i); z(i)];
    r2(:,i) = r1(:,i) + e2;
    r3(:,i) = r1(:,i) + e3;
    r4(:,i) = r1(:,i) + e4;
    r5(:,i) = r1(:,i) + e5;

    % Defining anchor points for each cable mechanism (Slightly simplified -
    % dealing with exit points of cables)
    anch_1_bot{i} = r1(:, i-1)+prevRot*t2;
    anch_1_top{i} = r1(:,i) + Rot*t2;
    anch_2_bot{i} = r1(:, i-1)+prevRot*t3;
    anch_2_top{i} = r1(:,i) + Rot*t3;
    anch_3_bot{i} = r1(:, i-1)+prevRot*t4;
    anch_3_top{i} = r1(:,i) + Rot*t4;
    anch_4_bot{i} = r1(:, i-1)+prevRot*t5;
    anch_4_top{i} = r1(:,i) + Rot*t5;
    angle = pi/6;
    anch_5_bot{i} = r1(:, i-1)+prevRot*t4;
    anch_5_top{i} = r1(:,i) + Rot*t2; % t4 to t2
    anch_6_bot{i} = r1(:, i-1)+prevRot*t5;
    anch_6_top{i} = r1(:,i) + Rot*t2; % t5 to t2
    anch_7_bot{i} = r1(:, i-1)+prevRot*t4;
    anch_7_top{i} = r1(:,i) + Rot*t3; % t4 to t3
    anch_8_bot{i} = r1(:, i-1)+prevRot*t5;
    anch_8_top{i} = r1(:,i) + Rot*t3; % t5 to t3

    % Nodal Velocities
    dr1(:,i) = fulldiff(r1(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
    dr2(:,i) = fulldiff(r2(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
    dr3(:,i) = fulldiff(r3(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
    dr4(:,i) = fulldiff(r4(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
    dr5(:,i) = fulldiff(r5(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});

    % Kinetic energy each tetrahedron
    T_L = simplify(1/2*m*(dr1(:,i).'*dr1(:,i) + dr2(:,i).'*dr2(:,i)+ ...
    dr3(:,i).'*dr3(:,i)+ dr4(:,i).'*dr4(:,i)+ dr5(:,i).'*dr5(:,i)));

    % Potential energy each tetrahedron
    V_L = simplify(m*g*[0 0 1]*(r1(:,i)+r2(:,i)+r3(:,i)+r4(:,i)+r5(:,i)));
    prevRot = Rot;
    
    % Lagrangian
    Lagr(i)=simplify(T_L-V_L);
end

disp('4%')

for i = 1:N % Build time derivatives for each vertebra
    %Calculate time derivitaves etc. according to lagrangian dynamics
    fx(i) =  fulldiff(diff(Lagr(i) , dx(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                    - diff(Lagr(i) , x(i));
    fy(i) = fulldiff(diff(Lagr(i) , dy(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                    - diff(Lagr(i) , y(i));
    fz(i) = fulldiff(diff(Lagr(i) , dz(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                    - diff(Lagr(i) , z(i));
    fT(i) = fulldiff(diff(Lagr(i) , dT(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                    - diff(Lagr(i) , T(i));
    fG(i) = fulldiff(diff(Lagr(i) , dG(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                    - diff(Lagr(i) , G(i));
    fP(i) =  fulldiff(diff(Lagr(i) , dP(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                    - diff(Lagr(i) , P(i));

    % Compute in parallel to speed things up
    pf1 = parfeval(pools,@simplify,1,fx(i),'Steps',100);
    pf2 = parfeval(pools,@simplify,1,fy(i),'Steps',100);
    pf3 = parfeval(pools,@simplify,1,fz(i),'Steps',100);
    pf4 = parfeval(pools,@simplify,1,fP(i),'Steps',100);
    pf5 = parfeval(pools,@simplify,1,fG(i),'Steps',100);
    pf6 = parfeval(pools,@simplify,1,fT(i),'Steps',100);

    % Get Ouptuts
    fx(i) = fetchOutputs(pf1);
    disp('5%')
    fy(i) = fetchOutputs(pf2);
    disp('6%')
    fz(i) = fetchOutputs(pf3);
    disp('7%')
    fP(i) = fetchOutputs(pf4);
    disp('8%')
    fG(i) = fetchOutputs(pf5);
    disp('9%')
    fT(i) = fetchOutputs(pf6);      
    disp('10%') 
end

% Again, building useless variables for link 1 to help with indexing
lengths(1, 1) = 0;
lengths(2, 1) = 0;
lengths(3, 1) = 0;
lengths(4, 1) = 0;
lengths(5, 1) = 0;
lengths(6, 1) = 0;
lengths(7, 1) = 0;
lengths(8, 1) = 0;
tensions(1, 1) = 0;
tensions(2, 1) = 0;
tensions(3, 1) = 0;
tensions(4, 1) = 0;
tensions(5, 1) = 0;
tensions(6, 1) = 0;
tensions(7, 1) = 0;
tensions(8, 1) = 0;

for i = 2:(N-1)

    % Spring lengths for springs below current tetrahedron  
    lengths(1, i) = norm(anch_1_bot{i}-anch_1_top{i});
    lengths(2, i) = norm(anch_2_bot{i}-anch_2_top{i});
    lengths(3, i) = norm(anch_3_bot{i}-anch_3_top{i});
    lengths(4, i) = norm(anch_4_bot{i}-anch_4_top{i});
    lengths(5, i) = norm(anch_5_bot{i}-anch_5_top{i});
    lengths(6, i) = norm(anch_6_bot{i}-anch_6_top{i});
    lengths(7, i) = norm(anch_7_bot{i}-anch_7_top{i});
    lengths(8, i) = norm(anch_8_bot{i}-anch_8_top{i});

    dlengths_dt(:, i) = fulldiff(lengths(:, i),{x(i),y(i),z(i),T(i),G(i),P(i)});
    parfor(j=1:8)
        dlengths_dt(j, i) = simplify(dlengths_dt(j, i),50);
    end

    % Compute in Parallel
    disp('11%')
    pf1 = parfeval(pools,@simplify,1,lengths(1, i),'Steps',100);
    pf2 = parfeval(pools,@simplify,1,lengths(2, i),'Steps',100);
    pf3 = parfeval(pools,@simplify,1,lengths(3, i),'Steps',100);
    pf4 = parfeval(pools,@simplify,1,lengths(4, i),'Steps',100);
    pf5 = parfeval(pools,@simplify,1,lengths(5, i),'Steps',100);
    pf6 = parfeval(pools,@simplify,1,lengths(6, i),'Steps',100);
    pf7 = parfeval(pools,@simplify,1,lengths(7, i),'Steps',100);
    pf8 = parfeval(pools,@simplify,1,lengths(8, i),'Steps',100);
    disp('12%')
    lengths(1, i) = fetchOutputs(pf1);
    disp('13%')
    lengths(2, i) = fetchOutputs(pf2);
    disp('14%')
    lengths(3, i) = fetchOutputs(pf3);
    disp('15%')
    lengths(4, i) = fetchOutputs(pf4);
    disp('16%')
    lengths(5, i) = fetchOutputs(pf5);
    disp('17%')
    lengths(6, i) = fetchOutputs(pf6);
    disp('18%')
    lengths(7, i) = fetchOutputs(pf7);
    disp('19%')
    lengths(8, i) = fetchOutputs(pf8);

    disp('20%')

    F1 = tensions(1, i)*(anch_1_bot{i}-anch_1_top{i}) - tensions(1, i+1)*(anch_1_bot{i+1}-anch_1_top{i+1});
    F2 = tensions(2, i)*(anch_2_bot{i}-anch_2_top{i}) - tensions(2, i+1)*(anch_2_bot{i+1}-anch_2_top{i+1});
    F3 = tensions(3, i)*(anch_3_bot{i}-anch_3_top{i}) - tensions(3, i+1)*(anch_3_bot{i+1}-anch_3_top{i+1});
    F4 = tensions(4, i)*(anch_4_bot{i}-anch_4_top{i}) - tensions(4, i+1)*(anch_4_bot{i+1}-anch_4_top{i+1});
    
    % Saddle Cables
    disp('21%')
    F5 = tensions(5, i)*(anch_5_bot{i}-anch_5_top{i});
    F6 = tensions(6, i)*(anch_6_bot{i}-anch_6_top{i});
    F7 = tensions(7, i)*(anch_7_bot{i}-anch_7_top{i});
    F8 = tensions(8, i)*(anch_8_bot{i}-anch_8_top{i});
    disp('22%')
    disp('23%')
    disp('24%')
    pf1 = parfeval(pools,@simplify,1,F1,'Steps',100);
    pf2 = parfeval(pools,@simplify,1,F2,'Steps',100);
    pf3 = parfeval(pools,@simplify,1,F3,'Steps',100);
    pf4 = parfeval(pools,@simplify,1,F4,'Steps',100);
    pf5 = parfeval(pools,@simplify,1,F5,'Steps',100);
    pf6 = parfeval(pools,@simplify,1,F6,'Steps',100);
    pf7 = parfeval(pools,@simplify,1,F7,'Steps',100);
    pf8 = parfeval(pools,@simplify,1,F8,'Steps',100);
    disp('25%')
    F1 = fetchOutputs(pf1);
    disp('26%')
    disp('27%')
    F2 = fetchOutputs(pf2);
    disp('28%')
    disp('29%')
    F3 = fetchOutputs(pf3);
    disp('30%')
    disp('31%')
    F4 = fetchOutputs(pf4);
    disp('32%')
    disp('33%')
    F5 = fetchOutputs(pf5);
    disp('34%')
    disp('35%')
    F6 = fetchOutputs(pf6);
    disp('36%')
    disp('37%')
    F7 = fetchOutputs(pf7);
    disp('38%')
    disp('39%')
    F8 = fetchOutputs(pf8);
    disp('40%')    
    % Vertical string forces

    % Forces in global coordinates
    Fx(i) = diff(anch_1_top{i},x(i)).'*(F1)...
        +diff(anch_2_top{i},x(i)).'*(F2)...
        +diff(anch_3_top{i},x(i)).'*(F3)...
        +diff(anch_4_top{i},x(i)).'*(F4)...
        +diff(anch_5_top{i},x(i)).'*(F5)...
        +diff(anch_6_top{i},x(i)).'*(F6)...
        +diff(anch_7_top{i},x(i)).'*(F7)...
        +diff(anch_8_top{i},x(i)).'*(F8);
    disp('42%')
    disp('43%')
    Fy(i) =  diff(anch_1_top{i},y(i)).'*(F1)...
        +diff(anch_2_top{i},y(i)).'*(F2)...
        +diff(anch_3_top{i},y(i)).'*(F3)...
        +diff(anch_4_top{i},y(i)).'*(F4)...
        +diff(anch_5_top{i},y(i)).'*(F5)...
        +diff(anch_6_top{i},y(i)).'*(F6)...
        +diff(anch_7_top{i},y(i)).'*(F7)...
        +diff(anch_8_top{i},y(i)).'*(F8);
    disp('44%')
    Fz(i) =  diff(anch_1_top{i},z(i)).'*(F1)...
        +diff(anch_2_top{i},z(i)).'*(F2)...
        +diff(anch_3_top{i},z(i)).'*(F3)...
        +diff(anch_4_top{i},z(i)).'*(F4)...
        +diff(anch_5_top{i},z(i)).'*(F5)...
        +diff(anch_6_top{i},z(i)).'*(F6)...
        +diff(anch_7_top{i},z(i)).'*(F7)...
        +diff(anch_8_top{i},z(i)).'*(F8);
    disp('45%')
    FT(i) =  diff(anch_1_top{i},T(i)).'*(F1)...
        +diff(anch_2_top{i},T(i)).'*(F2)...
        +diff(anch_3_top{i},T(i)).'*(F3)...
        +diff(anch_4_top{i},T(i)).'*(F4)...
        +diff(anch_5_top{i},T(i)).'*(F5)...
        +diff(anch_6_top{i},T(i)).'*(F6)...
        +diff(anch_7_top{i},T(i)).'*(F7)...
        +diff(anch_8_top{i},T(i)).'*(F8);
    disp('46%')
    FG(i) =  diff(anch_1_top{i},G(i)).'*(F1)...
        +diff(anch_2_top{i},G(i)).'*(F2)...
        +diff(anch_3_top{i},G(i)).'*(F3)...
        +diff(anch_4_top{i},G(i)).'*(F4)...
        +diff(anch_5_top{i},G(i)).'*(F5)...
        +diff(anch_6_top{i},G(i)).'*(F6)...
        +diff(anch_7_top{i},G(i)).'*(F7)...
        +diff(anch_8_top{i},G(i)).'*(F8);
    disp('47%')
    FP(i) =  diff(anch_1_top{i},P(i)).'*(F1)...
        +diff(anch_2_top{i},P(i)).'*(F2)...
        +diff(anch_3_top{i},P(i)).'*(F3)...
        +diff(anch_4_top{i},P(i)).'*(F4)...
        +diff(anch_5_top{i},P(i)).'*(F5)...
        +diff(anch_6_top{i},P(i)).'*(F6)...
        +diff(anch_7_top{i},P(i)).'*(F7)...
        +diff(anch_8_top{i},P(i)).'*(F8);
    disp('48%')
    pf1 = parfeval(pools,@simplify,1,Fx(i), 'Criterion','preferReal','Steps',100);
    pf2 = parfeval(pools,@simplify,1,Fy(i), 'Criterion','preferReal','Steps',100);
    pf3 = parfeval(pools,@simplify,1,Fz(i), 'Criterion','preferReal','Steps',100);
    pf4 = parfeval(pools,@simplify,1,FP(i), 'Criterion','preferReal','Steps',100);
    pf5 = parfeval(pools,@simplify,1,FG(i), 'Criterion','preferReal','Steps',100);
    pf6 = parfeval(pools,@simplify,1,FT(i), 'Criterion','preferReal','Steps',100);
    disp('49%')
    Fx(i) = fetchOutputs(pf1);
    disp('50%')
    Fy(i) = fetchOutputs(pf2);
    disp('52%')
    Fz(i) = fetchOutputs(pf3);
    disp('54%')
    FP(i) = fetchOutputs(pf4);
    disp('56%')
    FG(i) = fetchOutputs(pf5);
    disp('58%')
    FT(i) = fetchOutputs(pf6);
    disp('60%')
    [D2x(i), D2y(i), D2z(i), D2T(i), D2G(i), D2P(i)] = solve(Fx(i)==fx(i), Fy(i)==fy(i),Fz(i) == fz(i), FT(i)==fT(i),...
                                       FG(i)==fG(i),FP(i)==fP(i),d2x(i),d2y(i),d2z(i),d2T(i),d2G(i),d2P(i));
end

% Code for final tetrahedron at top of structure
% Spring lengths for springs below current tetrahedron  
lengths(1, N) = norm(anch_1_bot{N}-anch_1_top{N});
lengths(2, N) = norm(anch_2_bot{N}-anch_2_top{N});
lengths(3, N) = norm(anch_3_bot{N}-anch_3_top{N});
lengths(4, N) = norm(anch_4_bot{N}-anch_4_top{N});
lengths(5, N) = norm(anch_5_bot{N}-anch_5_top{N});
lengths(6, N) = norm(anch_6_bot{N}-anch_6_top{N});
lengths(7, N) = norm(anch_7_bot{N}-anch_7_top{N});
lengths(8, N) = norm(anch_8_bot{N}-anch_8_top{N});
    
dlengths_dt(:, N) = fulldiff(lengths(:, N),{x(N),y(N),z(N),T(N),G(N),P(N)});
parfor(j=1:8)
    dlengths_dt(j, N) = simplify(dlengths_dt(j, N),50);
end
    
% Compute in Parallel
disp('11%')
pf1 = parfeval(pools,@simplify,1,lengths(1, N),'Steps',100);
pf2 = parfeval(pools,@simplify,1,lengths(2, N),'Steps',100);
pf3 = parfeval(pools,@simplify,1,lengths(3, N),'Steps',100);
pf4 = parfeval(pools,@simplify,1,lengths(4, N),'Steps',100);
pf5 = parfeval(pools,@simplify,1,lengths(5, N),'Steps',100);
pf6 = parfeval(pools,@simplify,1,lengths(6, N),'Steps',100);
pf7 = parfeval(pools,@simplify,1,lengths(7, N),'Steps',100);
pf8 = parfeval(pools,@simplify,1,lengths(8, N),'Steps',100);
disp('12%')
lengths(1, N) = fetchOutputs(pf1);
disp('13%')
lengths(2, N) = fetchOutputs(pf2);
disp('14%')
lengths(3, N) = fetchOutputs(pf3);
disp('15%')
lengths(4, N) = fetchOutputs(pf4);
disp('16%')
lengths(5, N) = fetchOutputs(pf5);
disp('17%')
lengths(6, N) = fetchOutputs(pf6);
disp('18%')
lengths(7, N) = fetchOutputs(pf7);
disp('19%')
lengths(8, N) = fetchOutputs(pf8);

disp('20%')

F1 = tensions(1, N)*(anch_1_bot{N}-anch_1_top{N});
F2 = tensions(2, N)*(anch_2_bot{N}-anch_2_top{N});
F3 = tensions(3, N)*(anch_3_bot{N}-anch_3_top{N});
F4 = tensions(4, N)*(anch_4_bot{N}-anch_4_top{N});
% Saddle cable forces
disp('21%')
F5 = tensions(5, N)*(anch_5_bot{N}-anch_5_top{N});
F6 = tensions(6, N)*(anch_6_bot{N}-anch_6_top{N});
F7 = tensions(7, N)*(anch_7_bot{N}-anch_7_top{N});
F8 = tensions(8, N)*(anch_8_bot{N}-anch_8_top{N});
disp('22%')
disp('23%')
disp('24%')
pf1 = parfeval(pools,@simplify,1,F1,'Steps',100);
pf2 = parfeval(pools,@simplify,1,F2,'Steps',100);
pf3 = parfeval(pools,@simplify,1,F3,'Steps',100);
pf4 = parfeval(pools,@simplify,1,F4,'Steps',100);
pf5 = parfeval(pools,@simplify,1,F5,'Steps',100);
pf6 = parfeval(pools,@simplify,1,F6,'Steps',100);
pf7 = parfeval(pools,@simplify,1,F7,'Steps',100);
pf8 = parfeval(pools,@simplify,1,F8,'Steps',100);
disp('25%')
F1 = fetchOutputs(pf1);
disp('26%')
disp('27%')
F2 = fetchOutputs(pf2);
disp('28%')
disp('29%')
F3 = fetchOutputs(pf3);
disp('30%')
disp('31%')
F4 = fetchOutputs(pf4);
disp('32%')
disp('33%')
F5 = fetchOutputs(pf5);
disp('34%')
disp('35%')
F6 = fetchOutputs(pf6);
disp('36%')
disp('37%')
F7 = fetchOutputs(pf7);
disp('38%')
disp('39%')
F8 = fetchOutputs(pf8);
disp('40%')    
% Vertical string forces


% Forces in global coordinates
Fx(N) = diff(anch_1_top{N},x(N)).'*(F1)...
    +diff(anch_2_top{N},x(N)).'*(F2)...
    +diff(anch_3_top{N},x(N)).'*(F3)...
    +diff(anch_4_top{N},x(N)).'*(F4)...
    +diff(anch_5_top{N},x(N)).'*(F5)...
    +diff(anch_6_top{N},x(N)).'*(F6)...
    +diff(anch_7_top{N},x(N)).'*(F7)...
    +diff(anch_8_top{N},x(N)).'*(F8);
disp('42%')
disp('43%')
Fy(N) =  diff(anch_1_top{N},y(N)).'*(F1)...
    +diff(anch_2_top{N},y(N)).'*(F2)...
    +diff(anch_3_top{N},y(N)).'*(F3)...
    +diff(anch_4_top{N},y(N)).'*(F4)...
    +diff(anch_5_top{N},y(N)).'*(F5)...
    +diff(anch_6_top{N},y(N)).'*(F6)...
    +diff(anch_7_top{N},y(N)).'*(F7)...
    +diff(anch_8_top{N},y(N)).'*(F8);
disp('44%')
Fz(N) =  diff(anch_1_top{N},z(N)).'*(F1)...
    +diff(anch_2_top{N},z(N)).'*(F2)...
    +diff(anch_3_top{N},z(N)).'*(F3)...
    +diff(anch_4_top{N},z(N)).'*(F4)...
    +diff(anch_5_top{N},z(N)).'*(F5)...
    +diff(anch_6_top{N},z(N)).'*(F6)...
    +diff(anch_7_top{N},z(N)).'*(F7)...
    +diff(anch_8_top{N},z(N)).'*(F8);
disp('45%')
FT(N) =  diff(anch_1_top{N},T(N)).'*(F1)...
    +diff(anch_2_top{N},T(N)).'*(F2)...
    +diff(anch_3_top{N},T(N)).'*(F3)...
    +diff(anch_4_top{N},T(N)).'*(F4)...
    +diff(anch_5_top{N},T(N)).'*(F5)...
    +diff(anch_6_top{N},T(N)).'*(F6)...
    +diff(anch_7_top{N},T(N)).'*(F7)...
    +diff(anch_8_top{N},T(N)).'*(F8);
disp('46%')
FG(N) =  diff(anch_1_top{N},G(N)).'*(F1)...
    +diff(anch_2_top{N},G(N)).'*(F2)...
    +diff(anch_3_top{N},G(N)).'*(F3)...
    +diff(anch_4_top{N},G(N)).'*(F4)...
    +diff(anch_5_top{N},G(N)).'*(F5)...
    +diff(anch_6_top{N},G(N)).'*(F6)...
    +diff(anch_7_top{N},G(N)).'*(F7)...
    +diff(anch_8_top{N},G(N)).'*(F8);
disp('47%')
FP(N) =  diff(anch_1_top{N},P(N)).'*(F1)...
    +diff(anch_2_top{N},P(N)).'*(F2)...
    +diff(anch_3_top{N},P(N)).'*(F3)...
    +diff(anch_4_top{N},P(N)).'*(F4)...
    +diff(anch_5_top{N},P(N)).'*(F5)...
    +diff(anch_6_top{N},P(N)).'*(F6)...
    +diff(anch_7_top{N},P(N)).'*(F7)...
    +diff(anch_8_top{N},P(N)).'*(F8);
disp('48%')
pf1 = parfeval(pools,@simplify,1,Fx, 'Criterion','preferReal','Steps',100);
pf2 = parfeval(pools,@simplify,1,Fy, 'Criterion','preferReal','Steps',100);
pf3 = parfeval(pools,@simplify,1,Fz, 'Criterion','preferReal','Steps',100);
pf4 = parfeval(pools,@simplify,1,FP, 'Criterion','preferReal','Steps',100);
pf5 = parfeval(pools,@simplify,1,FG, 'Criterion','preferReal','Steps',100);
pf6 = parfeval(pools,@simplify,1,FT, 'Criterion','preferReal','Steps',100);
disp('49%')
Fx = fetchOutputs(pf1);
disp('50%')
Fy = fetchOutputs(pf2);
disp('52%')
Fz = fetchOutputs(pf3);
disp('54%')
FP = fetchOutputs(pf4);
disp('56%')
FG = fetchOutputs(pf5);
disp('58%')
FT = fetchOutputs(pf6);
disp('60%')
[D2x(N), D2y(N), D2z(N), D2T(N), D2G(N), D2P(N)] = solve(Fx(N)==fx(N), Fy(N)==fy(N),Fz(N) == fz(N), FT(N)==fT(N),...
                                   FG(N)==fG(N),FP(N)==fP(N),d2x(N),d2y(N),d2z(N),d2T(N),d2G(N),d2P(N));

%% EOM Computation
disp('61%')
disp('62%')
disp('63%')
disp('64%')
disp('65%')
disp('66%')
disp('67%')
disp('68%')
disp('69%')
pf1 = parfeval(pools,@simplify,1,D2x, 'Criterion','preferReal','Steps',100);
pf2 = parfeval(pools,@simplify,1,D2y, 'Criterion','preferReal','Steps',100);
pf3 = parfeval(pools,@simplify,1,D2z, 'Criterion','preferReal','Steps',100);
pf4 = parfeval(pools,@simplify,1,D2P, 'Criterion','preferReal','Steps',100);
pf5 = parfeval(pools,@simplify,1,D2G, 'Criterion','preferReal','Steps',100);
pf6 = parfeval(pools,@simplify,1,D2T, 'Criterion','preferReal','Steps',100);

D2x = fetchOutputs(pf1);
disp('70%')
D2y = fetchOutputs(pf2);
disp('74%')
D2z = fetchOutputs(pf3);
disp('78%')
D2P = fetchOutputs(pf4);
disp('82%')
D2G = fetchOutputs(pf5);
disp('86%')
D2T = fetchOutputs(pf6);
disp('90%')

%% Save the spine geometric parameters to a .mat file for use later
save(spine_geometric_parameters_path, 'spine_geometric_parameters');

%% Generate MATLAB functions
Dyn_eqn = [D2x(2); D2y(2); D2z(2); D2T(2); D2G(2); D2P(2); D2x(3); D2y(3); D2z(3); D2T(3); D2G(3); D2P(3); D2x(4); D2y(4); D2z(4); D2T(4); D2G(4); D2P(4)];

matlabFunction(Dyn_eqn,'file','duct_accel','Vars',[x(2); y(2); z(2); T(2); G(2); P(2); dx(2); dy(2); dz(2); dT(2); dG(2); dP(2); tensions(:, 2); ...
    x(3); y(3); z(3); T(3); G(3); P(3); dx(3); dy(3); dz(3); dT(3); dG(3); dP(3); tensions(:, 3); ...
    x(4); y(4); z(4); T(4); G(4); P(4); dx(4); dy(4); dz(4); dT(4); dG(4); dP(4); tensions(:, 4)]);
matlabFunction(lengths,'file','lengths','Vars',[x(2); y(2); z(2); T(2); G(2); P(2); x(3); y(3); z(3); T(3); G(3); P(3); x(4); y(4); z(4); T(4); G(4); P(4)]);
matlabFunction(dlengths_dt,'file','dlengths_dt','Vars',[x(2); y(2); z(2); T(2); G(2); P(2); dx(2); dy(2); dz(2); dT(2); dG(2); dP(2); ...
    x(3); y(3); z(3); T(3); G(3); P(3); dx(3); dy(3); dz(3); dT(3); dG(3); dP(3);
    x(4); y(4); z(4); T(4); G(4); P(4); dx(4); dy(4); dz(4); dT(4); dG(4); dP(4)]);
disp('100%')

 
