% two_d_dynamics_JeffMethod.m
% Andrew P. Sabelhaus (adapted from Jeff Friesen)
% Copyright 2016
% This script generates multiple m-files:
%   duct_accel
%   dlengths_dt
%   lengths


clc;
clear variables;
close all;

% As of 2016-11-13, fulldiff uses some functionality
% that apparently will be deprecated in a future release of MATLAB.
% For now, turn off that warning.
% TO-DO: re-write fulldiff!
warning off symbolic:sym:DeprecateFindsym

disp('Running the 2D version of Jeff Friesens dynamics derivation...');

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
N = 2;
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
spine_geometric_parameters_path = 'spine_geometric_parameters_2D';

% Creating Parallel Pools for Computing
pools = gcp;

%% Building symbolic matrices for computation
% X, Z, Theta of tetrahedrons for N links
x = sym('x',[N 1]);
z = sym('z',[N 1]);
T = sym('T',[N 1]);

% Velocities for each of above variables
dx = sym('dx',[N 1]);
dz = sym('dz',[N 1]);
dT = sym('dT',[N 1]);

% Accelerations for each of above variables
d2x = sym('d2x',[N 1]);
d2z = sym('d2z',[N 1]);
d2T = sym('d2T',[N 1]);

% Solved accelerations (EOM)
D2x = sym('D2x',[N 1]);
D2z = sym('D2z',[N 1]);
D2T = sym('D2T',[N 1]);

% Lagrangian derivatives for each variable
fx = sym('fx',[N 1]);
fz = sym('fz',[N 1]);
fT = sym('fT',[N 1]);

% Tensions of Cables
tensions = sym('tensions',[4 N]);
lengths = sym('lengths', [4 N]);
dlengths_dt = sym('dlengths_dt', [4 N]);
dlengths_dt(:, 1) = zeros(4, 1);

% Nodal positions of each link from fixed frame
% t1 = [0, 0, 0]';
% t2 = [(l^2 - (h/2)^2)^.5, 0, -h/2]';
% t3 = [-(l^2 - (h/2)^2)^.5, 0, -h/2]';
% t4 = [0, (l^2 - (h/2)^2)^.5, h/2]';
% t5 = [0, -(l^2 - (h/2)^2)^.5, h/2]';

% Now, with the 2D version: these are in (x,z).
t1 = [0, 0]';
t2 = [ -(l^2 - (h/2)^2)^.5, -h/2]';
t3 = [ (l^2 - (h/2)^2)^.5, -h/2]';
t4 = [0, h/2]';

% Nodal positons of each vertex (N by 2)
r1 = sym('r1',[2 N]);
r2 = sym('r2',[2 N]);
r3 = sym('r3',[2 N]);
r4 = sym('r4',[2 N]);

% Nodal velocities of each node for each spine link
dr1 = sym('dr1',[2 N]);
dr2 = sym('dr2',[2 N]);
dr3 = sym('dr3',[2 N]);
dr4 = sym('dr4',[2 N]);

% First link is fixed in place
r1(:,1) = [0; 0]; % Base link centered around (0,0,0)
r2(:,1) = r1(:,1) + t2; % No rotations
r3(:,1) = r1(:,1) + t3; 
r4(:,1) = r1(:,1) + t4;

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
assume(z,'real')
assume(T,'real')
assume(dx,'real')
assume(dz,'real')
assume(dT,'real')

% Anchor points for first link - these don't matter they just simplify
% indexing
anch_1_bot{1} = [0; 0];
anch_1_top{1} = [0; 0];
anch_2_bot{1} = [0; 0];
anch_2_top{1} = [0; 0];
anch_3_bot{1} = [0; 0];
anch_3_top{1} = [0; 0];
anch_4_bot{1} = [0; 0];
anch_4_top{1} = [0; 0];
angle = pi/6;

prevRot = eye(2);

Lagr = sym('Lagr',[1 N]);
for i=2:N
    % Rotation matrices for theta gamma and phi (euler rotations)

    R1 = [ cos(T(i)),    -sin(T(i));
           sin(T(i)),     cos(T(i))];

    disp('2%')

    % Building full Rotation Matrix
    % For the 2D case, there is only one rotation matrix.
    Rot=R1;

    % Re-initialize the t_i points (why???)
    t1 = [0, 0]';
    t2 = [ -(l^2 - (h/2)^2)^.5, -h/2]';
    t3 = [ (l^2 - (h/2)^2)^.5, -h/2]';
    t4 = [0, h/2]';


    % Compute nodal positions after a given rotation
    e2 = Rot*t2;
    e3 = Rot*t3;
    e4 = Rot*t4;

    disp('3%')
    % Nodal Positions
    r1(:,i) = [x(i); z(i)];
    r2(:,i) = r1(:,i) + e2;
    r3(:,i) = r1(:,i) + e3;
    r4(:,i) = r1(:,i) + e4;

    % Defining anchor points for each cable mechanism (Slightly simplified -
    % dealing with exit points of cables)
    anch_1_bot{i} = r1(:, i-1)+prevRot*t2; % vertical from t2 to t2 next
    anch_1_top{i} = r1(:,i) + Rot*t2;
    anch_2_bot{i} = r1(:, i-1)+prevRot*t3; % vertical from t3 to t3 next
    anch_2_top{i} = r1(:,i) + Rot*t3;
    anch_3_bot{i} = r1(:, i-1)+prevRot*t4; % saddle from t4 to t2 next
    anch_3_top{i} = r1(:,i) + Rot*t2;
    anch_4_bot{i} = r1(:, i-1)+prevRot*t4; % saddle from t4 to t3 next
    anch_4_top{i} = r1(:,i) + Rot*t3;
    angle = pi/6;

    % Nodal Velocities
    dr1(:,i) = fulldiff(r1(:,i),{x(i),z(i),T(i)});
    dr2(:,i) = fulldiff(r2(:,i),{x(i),z(i),T(i)});
    dr3(:,i) = fulldiff(r3(:,i),{x(i),z(i),T(i)});
    dr4(:,i) = fulldiff(r4(:,i),{x(i),z(i),T(i)});

    % Kinetic energy each tetrahedron
    T_L = simplify(1/2*m*(dr1(:,i).'*dr1(:,i) + dr2(:,i).'*dr2(:,i)+ ...
    dr3(:,i).'*dr3(:,i)+ dr4(:,i).'*dr4(:,i)));

    % Potential energy each tetrahedron
    V_L = simplify(m*g*[0 1]*(r1(:,i)+r2(:,i)+r3(:,i)+r4(:,i)));
    prevRot = Rot;
    
    % Lagrangian
    Lagr(i)=simplify(T_L-V_L);
end

disp('4%')

for i = 1:N % Build time derivatives for each vertebra
    %Calculate time derivitaves etc. according to lagrangian dynamics
    fx(i) =  fulldiff(diff(Lagr(i) , dx(i)),{x(i),z(i),T(i)})...
                    - diff(Lagr(i) , x(i));
    fz(i) = fulldiff(diff(Lagr(i) , dz(i)),{x(i),z(i),T(i)})...
                    - diff(Lagr(i) , z(i));
    fT(i) = fulldiff(diff(Lagr(i) , dT(i)),{x(i),z(i),T(i)})...
                    - diff(Lagr(i) , T(i));

    % Compute in parallel to speed things up
    pf1 = parfeval(pools,@simplify,1,fx(i),'Steps',100);
    pf2 = parfeval(pools,@simplify,1,fz(i),'Steps',100);
    pf3 = parfeval(pools,@simplify,1,fT(i),'Steps',100);

    % Get Ouptuts
    fx(i) = fetchOutputs(pf1);
    disp('6%')
    fz(i) = fetchOutputs(pf2);
    disp('7%')
    fT(i) = fetchOutputs(pf3);      
    disp('10%') 
end

% Again, building useless variables for link 1 to help with indexing
lengths(1, 1) = 0;
lengths(2, 1) = 0;
lengths(3, 1) = 0;
lengths(4, 1) = 0;
tensions(1, 1) = 0;
tensions(2, 1) = 0;
tensions(3, 1) = 0;
tensions(4, 1) = 0;


for i = 2:(N-1)

    % Spring lengths for springs below current tetrahedron  
    lengths(1, i) = norm(anch_1_bot{i}-anch_1_top{i});
    lengths(2, i) = norm(anch_2_bot{i}-anch_2_top{i});
    lengths(3, i) = norm(anch_3_bot{i}-anch_3_top{i});
    lengths(4, i) = norm(anch_4_bot{i}-anch_4_top{i});

    dlengths_dt(:, i) = fulldiff(lengths(:, i),{x(i),z(i),T(i)});
    parfor(j=1:4)
        dlengths_dt(j, i) = simplify(dlengths_dt(j, i),50);
    end

    % Compute in Parallel
    disp('11%')
    pf1 = parfeval(pools,@simplify,1,lengths(1, i),'Steps',100);
    pf2 = parfeval(pools,@simplify,1,lengths(2, i),'Steps',100);
    pf3 = parfeval(pools,@simplify,1,lengths(3, i),'Steps',100);
    pf4 = parfeval(pools,@simplify,1,lengths(4, i),'Steps',100);
    disp('12%')
    lengths(1, i) = fetchOutputs(pf1);
    disp('13%')
    lengths(2, i) = fetchOutputs(pf2);
    disp('14%')
    lengths(3, i) = fetchOutputs(pf3);
    disp('15%')
    lengths(4, i) = fetchOutputs(pf4);

    disp('20%')

    % Forces due to the vertical cables
    F1 = tensions(1, i)*(anch_1_bot{i}-anch_1_top{i}) - tensions(1, i+1)*(anch_1_bot{i+1}-anch_1_top{i+1});
    F2 = tensions(2, i)*(anch_2_bot{i}-anch_2_top{i}) - tensions(2, i+1)*(anch_2_bot{i+1}-anch_2_top{i+1});
    
    % Saddle Cables
    disp('21%')
    F3 = tensions(3, i)*(anch_3_bot{i}-anch_3_top{i});
    F4 = tensions(4, i)*(anch_4_bot{i}-anch_4_top{i});

    disp('22%')
    disp('23%')
    disp('24%')
    pf1 = parfeval(pools,@simplify,1,F1,'Steps',100);
    pf2 = parfeval(pools,@simplify,1,F2,'Steps',100);
    pf3 = parfeval(pools,@simplify,1,F3,'Steps',100);
    pf4 = parfeval(pools,@simplify,1,F4,'Steps',100);
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

    % Vertical string forces

    % Forces in global coordinates
    Fx(i) = diff(anch_1_top{i},x(i)).'*(F1)...
        +diff(anch_2_top{i},x(i)).'*(F2)...
        +diff(anch_3_top{i},x(i)).'*(F3)...
        +diff(anch_4_top{i},x(i)).'*(F4);
    disp('42%')

    Fz(i) =  diff(anch_1_top{i},z(i)).'*(F1)...
        +diff(anch_2_top{i},z(i)).'*(F2)...
        +diff(anch_3_top{i},z(i)).'*(F3)...
        +diff(anch_4_top{i},z(i)).'*(F4);
    disp('45%')
    FT(i) =  diff(anch_1_top{i},T(i)).'*(F1)...
        +diff(anch_2_top{i},T(i)).'*(F2)...
        +diff(anch_3_top{i},T(i)).'*(F3)...
        +diff(anch_4_top{i},T(i)).'*(F4);

    disp('48%')
    pf1 = parfeval(pools,@simplify,1,Fx(i), 'Criterion','preferReal','Steps',100);
    pf3 = parfeval(pools,@simplify,1,Fz(i), 'Criterion','preferReal','Steps',100);
    pf6 = parfeval(pools,@simplify,1,FT(i), 'Criterion','preferReal','Steps',100);
    disp('49%')
    Fx(i) = fetchOutputs(pf1);
    disp('50%')
    Fz(i) = fetchOutputs(pf3);
    disp('54%')
    FT(i) = fetchOutputs(pf6);
    disp('60%')
    [D2x(i), D2z(i), D2T(i)] = solve(Fx(i)==fx(i), Fz(i) == fz(i), FT(i)==fT(i),...
                                      d2x(i),d2z(i),d2T(i));
end

% Code for final tetrahedron at top of structure
% Spring lengths for springs below current tetrahedron  
lengths(1, N) = norm(anch_1_bot{N}-anch_1_top{N});
lengths(2, N) = norm(anch_2_bot{N}-anch_2_top{N});
lengths(3, N) = norm(anch_3_bot{N}-anch_3_top{N});
lengths(4, N) = norm(anch_4_bot{N}-anch_4_top{N});
    
dlengths_dt(:, N) = fulldiff(lengths(:, N),{x(N),z(N),T(N)});
parfor(j=1:4)
    dlengths_dt(j, N) = simplify(dlengths_dt(j, N),50);
end
    
% Compute in Parallel
disp('11%')
pf1 = parfeval(pools,@simplify,1,lengths(1, N),'Steps',100);
pf2 = parfeval(pools,@simplify,1,lengths(2, N),'Steps',100);
pf3 = parfeval(pools,@simplify,1,lengths(3, N),'Steps',100);
pf4 = parfeval(pools,@simplify,1,lengths(4, N),'Steps',100);
disp('12%')
lengths(1, N) = fetchOutputs(pf1);
disp('13%')
lengths(2, N) = fetchOutputs(pf2);
disp('14%')
lengths(3, N) = fetchOutputs(pf3);
disp('15%')
lengths(4, N) = fetchOutputs(pf4);

disp('20%')

F1 = tensions(1, N)*(anch_1_bot{N}-anch_1_top{N});
F2 = tensions(2, N)*(anch_2_bot{N}-anch_2_top{N});
% Saddle cable forces
disp('21%')
F3 = tensions(3, N)*(anch_3_bot{N}-anch_3_top{N});
F4 = tensions(4, N)*(anch_4_bot{N}-anch_4_top{N});


disp('22%')
disp('23%')
disp('24%')
pf1 = parfeval(pools,@simplify,1,F1,'Steps',100);
pf2 = parfeval(pools,@simplify,1,F2,'Steps',100);
pf3 = parfeval(pools,@simplify,1,F3,'Steps',100);
pf4 = parfeval(pools,@simplify,1,F4,'Steps',100);
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
disp('40%')    
% Vertical string forces

% Forces in global coordinates
Fx(N) = diff(anch_1_top{N},x(N)).'*(F1)...
    +diff(anch_2_top{N},x(N)).'*(F2)...
    +diff(anch_3_top{N},x(N)).'*(F3)...
    +diff(anch_4_top{N},x(N)).'*(F4);
disp('42%')

Fz(N) =  diff(anch_1_top{N},z(N)).'*(F1)...
    +diff(anch_2_top{N},z(N)).'*(F2)...
    +diff(anch_3_top{N},z(N)).'*(F3)...
    +diff(anch_4_top{N},z(N)).'*(F4);
disp('45%')
FT(N) =  diff(anch_1_top{N},T(N)).'*(F1)...
    +diff(anch_2_top{N},T(N)).'*(F2)...
    +diff(anch_3_top{N},T(N)).'*(F3)...
    +diff(anch_4_top{N},T(N)).'*(F4);

disp('48%')
pf1 = parfeval(pools,@simplify,1,Fx(N), 'Criterion','preferReal','Steps',100);
pf3 = parfeval(pools,@simplify,1,Fz(N), 'Criterion','preferReal','Steps',100);
pf6 = parfeval(pools,@simplify,1,FT(N), 'Criterion','preferReal','Steps',100);
disp('49%')
Fx(N) = fetchOutputs(pf1);
disp('50%')
Fz(N) = fetchOutputs(pf3);
disp('54%')
FT(N) = fetchOutputs(pf6);
disp('60%')
[D2x(N), D2z(N), D2T(N)] = solve(Fx(N)==fx(N), Fz(N) == fz(N), FT(N)==fT(N),...
                                  d2x(N),d2z(N),d2T(N));


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
pf3 = parfeval(pools,@simplify,1,D2z, 'Criterion','preferReal','Steps',100);
pf6 = parfeval(pools,@simplify,1,D2T, 'Criterion','preferReal','Steps',100);

D2x = fetchOutputs(pf1);
disp('74%')
D2z = fetchOutputs(pf3);
disp('82%')
D2T = fetchOutputs(pf6);
disp('90%')

%% Save the spine geometric parameters to a .mat file for use later
save(spine_geometric_parameters_path, 'spine_geometric_parameters');

%% Generate MATLAB functions

% NOTE: NEED TO ADJUST THE INDICES HERE IF N IS CHANGED.

% For 3 vertebrae:
% Dyn_eqn = [D2x(2); D2z(2); D2T(2); D2x(3); D2z(3); D2T(3)];
% matlabFunction(Dyn_eqn,'file','two_d_spine_accel','Vars',[x(2); z(2); T(2); dx(2); dz(2); dT(2); tensions(:, 2); ...
%     x(3); z(3); T(3); dx(3); dz(3); dT(3); tensions(:, 3)]);
% matlabFunction(lengths,'file','two_d_spine_lengths','Vars',[x(2); z(2); T(2); x(3); z(3); T(3)]);
% matlabFunction(dlengths_dt,'file','two_d_spine_dlengths_dt','Vars',[x(2); z(2); T(2); dx(2); dz(2); dT(2); ...
%     x(3); z(3); T(3); dx(3); dz(3); dT(3)]);

% For 2 vertebrae:
Dyn_eqn = [D2x(2); D2z(2); D2T(2)];
matlabFunction(Dyn_eqn,'file','two_d_spine_accel','Vars',[x(2); z(2); T(2); dx(2); dz(2); dT(2); tensions(:, 2)]);
matlabFunction(lengths,'file','two_d_spine_lengths','Vars',[x(2); z(2); T(2)]);
matlabFunction(dlengths_dt,'file','two_d_spine_dlengths_dt','Vars',[x(2); z(2); T(2); dx(2); dz(2); dT(2)]);
disp('100%')

 








