% @Author Jeff Friesen
% Dynamic Tensegrity Robotics Lab
% Intelligent Robotics Group, NASA Ames Research Center
% Modified 5/20/2015
% Contact Jeff at: jfriesen222@gmail.com
%Dynamics for a flemons spine model with 2 links

clc
clear variables
close all
%number of tetrahedra (still working on making this work for more than 2
%tetrahedra)
N=2;
% l = length of tetrahedron side
% h = height tetrahedon
% m is the mass of one particle, 5 particles per tetrahedron one at each
% corner node and one in the center. 

rad = 0.015;
g = 9.81;
h = 23/100;
l = 15.9*2/100;
lA1 = l;
lA2 = l;
m = (0.4536 + 4*0.12 + 0.6)/5;

%create parrallel pools for fast computing
pools = gcp;

%x y z theta gamma and phi of tetrahedrons, ith entry corresponds to ith
%tetrahedron
x = sym('x',[N 1]);
y = sym('y',[N 1]);
z = sym('z',[N 1]);
T = sym('T',[N 1]);
G = sym('G',[N 1]);
P = sym('P',[N 1]);

%respective velocities of the above variables
dx = sym('dx',[N 1]);
dy = sym('dy',[N 1]);
dz = sym('dz',[N 1]);
dT = sym('dT',[N 1]);
dG = sym('dG',[N 1]);
dP = sym('dP',[N 1]);

%respective accelerations of the above variables
d2x = sym('d2x',[N 1]);
d2y = sym('d2y',[N 1]);
d2z = sym('d2z',[N 1]);
d2T = sym('d2T',[N 1]);
d2G = sym('d2G',[N 1]);
d2P = sym('d2P',[N 1]);

%accelerations once solved
D2x = sym('D2x',[N 1]);
D2y = sym('D2y',[N 1]);
D2z = sym('D2z',[N 1]);
D2T = sym('D2T',[N 1]);
D2G = sym('D2G',[N 1]);
D2P = sym('D2P',[N 1]);

%lagrangian derivative equations for the respective variables
fx = sym('fx',[N 1]);
fy = sym('fy',[N 1]);
fz = sym('fz',[N 1]);
fT = sym('fT',[N 1]);
fG = sym('fG',[N 1]);
fP = sym('fP',[N 1]);


%tensions of cables
tensions = sym('tensions',[8 1]);

%Nodal positions of a tetrahedron from a fixed frame
t2 = [ 0;  -l/2; -h/2];
t3 = [ 0;   l/2; -h/2];
ta2 = [ 0;  -lA1/2; -h/2];
ta3 = [ 0;   lA1/2; -h/2];

t4 = [-l/2;   0;  h/2];
t5 = [ l/2;   0;  h/2];
ta4 = [-lA2/2;   0;  h/2];
ta5 = [ lA2/2;   0;  h/2];



%Nodal positons of the tetrahedra, ith column correspons to ith cartesian
%coordinates of a node
r1 = sym('r1',[3 N]);
r2 = sym('r2',[3 N]);
r3 = sym('r3',[3 N]);
r4 = sym('r4',[3 N]);
r5 = sym('r5',[3 N]);

%Coordinates for shifting location of cable attachments

top_bar_prerot = [0; 0;  rad];
bot_bar_prerot = [0; 0; -rad];

%Nodal velocities of the tetrahedra, ith column correspons to ith cartesian
%coordinates of a node
dr1 = sym('dr1',[3 N]);
dr2 = sym('dr2',[3 N]);
dr3 = sym('dr3',[3 N]);
dr4 = sym('dr4',[3 N]);
dr5 = sym('dr5',[3 N]);

%The first tetrahedron is fixed in place below, otherwise
%comment this out and make i go from 1:N instead of 2:N and add constraint
%forces as desired later
r1(:,1) = [0; 0; 0]; %Base tetrahedron centered around 0,0,0
r2(:,1) = r1(:,1) + ta2;%No rotations
r3(:,1) = r1(:,1) + ta3; %actuator is on bottom
r4(:,1) = r1(:,1) + t4;
r5(:,1) = r1(:,1) + t5;
dr1(:,1) = 0;%Nodal velocities are zero for fixed tetrahedron
dr2(:,1) = 0;
dr3(:,1) = 0;
dr4(:,1) = 0;
dr5(:,1) = 0;
disp('1%')
%kinetic energy
T_L = 0;
%potential energy
V_L=0;

%following loop is used to calculate potential and kinetic energy of each
%tetrahedron
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

for i=2:N
%Rotation matrices for theta gamma and phi (euler rotations)    
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
%Combine rotations by multiplying   
Rot=R3*R2*R1;

t2 = [ 0;  -l/2; -h/2];
t3 = [ 0;   l/2; -h/2];
t4 = [-l/2;   0;  h/2];
t5 = [ l/2;   0;  h/2];


%Compute nodal positions after a given rotation
e2 = Rot*t2;
e3 = Rot*t3;
e4 = Rot*ta4;
e5 = Rot*ta5;
% top_bar(:,i) = Rot*top_bar_prerot;
% bot_bar(:,i) = Rot*bot_bar_prerot;

disp('3%')
%Nodal Positions
r1(:,i) = [x(i); y(i); z(i)];
r2(:,i) = r1(:,i) + e2;
r3(:,i) = r1(:,i) + e3;
r4(:,i) = r1(:,i) + e4;
r5(:,i) = r1(:,i) + e5;

%define string anchor positions
anch_1_bot = t2 + [0; 0; 3/4]*rad;
anch_1_top = r1(:,i) + Rot*(t2 + [0; 0; -1]*rad);
anch_2_bot = t3 + [0; 0;  3/4]*rad;
anch_2_top = r1(:,i) + Rot*(t3 + [0; 0; -1]*rad);
anch_3_bot = t4 + [0; 0;  1]*rad;
anch_3_top = r1(:,i) + Rot*(t4 + [0; 0; -3/4]*rad);
anch_4_bot = t5 + [0; 0;  1]*rad;
anch_4_top = r1(:,i) + Rot*(t5 + [0; 0; -3/4]*rad);
angle = pi/6;
anch_5_bot = t4 +                [0;          -sin(angle); -cos(angle)]*rad;
anch_5_top = r1(:,i) + Rot*(t2 + [ -sin(angle);         0;  cos(angle)]*rad);
anch_6_bot = t5 +                [0;          -sin(angle); -cos(angle)]*rad;
anch_6_top = r1(:,i) + Rot*(t2 + [sin(angle);           0;  cos(angle)]*rad);
anch_7_bot = t4 +                [0;           sin(angle); -cos(angle)]*rad;
anch_7_top = r1(:,i) + Rot*(t3 + [ -sin(angle);         0;  cos(angle)]*rad);
anch_8_bot = t5 +                [0;           sin(angle); -cos(angle)]*rad;
anch_8_top = r1(:,i) + Rot*(t3 + [sin(angle);           0;  cos(angle)]*rad);

%Nodal Velocities
dr1(:,i) = fulldiff(r1(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
dr2(:,i) = fulldiff(r2(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
dr3(:,i) = fulldiff(r3(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
dr4(:,i) = fulldiff(r4(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
dr5(:,i) = fulldiff(r5(:,i),{x(i),y(i),z(i),T(i),G(i),P(i)});
%Kinetic energy each tetrahedron
T_L = simplify(T_L + 1/2*m*(dr1(:,i).'*dr1(:,i) + dr2(:,i).'*dr2(:,i)+ ...
dr3(:,i).'*dr3(:,i)+ dr4(:,i).'*dr4(:,i)+ dr5(:,i).'*dr5(:,i)));
%Potential energy each tetrahedron
V_L = simplify(V_L + m*g*[0 0 1]*(r1(:,i)+r2(:,i)+r3(:,i)+r4(:,i)+r5(:,i)));
end

top_bar(:,1) = top_bar_prerot;
bot_bar(:,1) = bot_bar_prerot;
dtop_bar(:,1) = [0;0;0];
dbot_bar(:,1) = [0;0;0];

%Lagrangian
Lagr=simplify(T_L-V_L);
disp('4%')


    %Calculate time derivitaves etc. according to lagrangian dynamics
fx(i) =  fulldiff(diff(Lagr , dx(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                - diff(Lagr , x(i));
fy(i) = fulldiff(diff(Lagr , dy(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                - diff(Lagr , y(i));
fz(i) = fulldiff(diff(Lagr , dz(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                - diff(Lagr , z(i));
fT(i) = fulldiff(diff(Lagr , dT(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                - diff(Lagr , T(i));
fG(i) = fulldiff(diff(Lagr , dG(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                - diff(Lagr , G(i));
fP(i) =  fulldiff(diff(Lagr , dP(i)),{x(i),y(i),z(i),T(i),G(i),P(i)})...
                - diff(Lagr , P(i));
            %Compute in parrallel to speed things up
pf1 = parfeval(pools,@simplify,1,fx(i),'Steps',100);
pf2 = parfeval(pools,@simplify,1,fy(i),'Steps',100);
pf3 = parfeval(pools,@simplify,1,fz(i),'Steps',100);
pf4 = parfeval(pools,@simplify,1,fP(i),'Steps',100);
pf5 = parfeval(pools,@simplify,1,fG(i),'Steps',100);
pf6 = parfeval(pools,@simplify,1,fT(i),'Steps',100);
%get outputs
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

%Spring lengths for springs below current tetrahedron  
    lengths(1) = norm(anch_1_bot-anch_1_top);
    lengths(2) = norm(anch_2_bot-anch_2_top);
    lengths(3) = norm(anch_3_bot-anch_3_top);
    lengths(4) = norm(anch_4_bot-anch_4_top);
    lengths(5) = norm(anch_5_bot-anch_5_top);
    lengths(6) = norm(anch_6_bot-anch_6_top);
    lengths(7) = norm(anch_7_bot-anch_7_top);
    lengths(8) = norm(anch_8_bot-anch_8_top);
    
    dlengths_dt = fulldiff(lengths,{x(i),y(i),z(i),T(i),G(i),P(i)});
    parfor(j=1:8)
        
        dlengths_dt(j) = simplify(dlengths_dt(j),50);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %again compute parrallel
    disp('11%')
    pf1 = parfeval(pools,@simplify,1,lengths(1),'Steps',100);
    pf2 = parfeval(pools,@simplify,1,lengths(2),'Steps',100);
    pf3 = parfeval(pools,@simplify,1,lengths(3),'Steps',100);
    pf4 = parfeval(pools,@simplify,1,lengths(4),'Steps',100);
    pf5 = parfeval(pools,@simplify,1,lengths(5),'Steps',100);
    pf6 = parfeval(pools,@simplify,1,lengths(6),'Steps',100);
    pf7 = parfeval(pools,@simplify,1,lengths(7),'Steps',100);
    pf8 = parfeval(pools,@simplify,1,lengths(8),'Steps',100);
    disp('12%')
    lengths(1) = fetchOutputs(pf1);
    disp('13%')
	lengths(2) = fetchOutputs(pf2);
    disp('14%')
    lengths(3) = fetchOutputs(pf3);
    disp('15%')
    lengths(4) = fetchOutputs(pf4);
    disp('16%')
    lengths(5) = fetchOutputs(pf5);
    disp('17%')
    lengths(6) = fetchOutputs(pf6);
    disp('18%')
    lengths(7) = fetchOutputs(pf7);
    disp('19%')
    lengths(8) = fetchOutputs(pf8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('20%')

    F1 = tensions(1)*(anch_1_bot-anch_1_top);
    F2 = tensions(2)*(anch_2_bot-anch_2_top);
    F3 = tensions(3)*(anch_3_bot-anch_3_top);
    F4 = tensions(4)*(anch_4_bot-anch_4_top);
    %saddle string forces
    disp('21%')
    F5 = tensions(5)*(anch_5_bot-anch_5_top);
    F6 = tensions(6)*(anch_6_bot-anch_6_top);
    F7 = tensions(7)*(anch_7_bot-anch_7_top);
    F8 = tensions(8)*(anch_8_bot-anch_8_top);
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
%vertical string forces


%Forces in global coordinates
Fx = diff(anch_1_top,x(i)).'*(F1)...
    +diff(anch_2_top,x(i)).'*(F2)...
    +diff(anch_3_top,x(i)).'*(F3)...
    +diff(anch_4_top,x(i)).'*(F4)...
    +diff(anch_5_top,x(i)).'*(F5)...
    +diff(anch_6_top,x(i)).'*(F6)...
    +diff(anch_7_top,x(i)).'*(F7)...
    +diff(anch_8_top,x(i)).'*(F8);
disp('42%')
disp('43%')
Fy =  diff(anch_1_top,y(i)).'*(F1)...
    +diff(anch_2_top,y(i)).'*(F2)...
    +diff(anch_3_top,y(i)).'*(F3)...
    +diff(anch_4_top,y(i)).'*(F4)...
    +diff(anch_5_top,y(i)).'*(F5)...
    +diff(anch_6_top,y(i)).'*(F6)...
    +diff(anch_7_top,y(i)).'*(F7)...
    +diff(anch_8_top,y(i)).'*(F8);
disp('44%')
Fz =  diff(anch_1_top,z(i)).'*(F1)...
    +diff(anch_2_top,z(i)).'*(F2)...
    +diff(anch_3_top,z(i)).'*(F3)...
    +diff(anch_4_top,z(i)).'*(F4)...
    +diff(anch_5_top,z(i)).'*(F5)...
    +diff(anch_6_top,z(i)).'*(F6)...
    +diff(anch_7_top,z(i)).'*(F7)...
    +diff(anch_8_top,z(i)).'*(F8);
disp('45%')
FT =  diff(anch_1_top,T(i)).'*(F1)...
    +diff(anch_2_top,T(i)).'*(F2)...
    +diff(anch_3_top,T(i)).'*(F3)...
    +diff(anch_4_top,T(i)).'*(F4)...
    +diff(anch_5_top,T(i)).'*(F5)...
    +diff(anch_6_top,T(i)).'*(F6)...
    +diff(anch_7_top,T(i)).'*(F7)...
    +diff(anch_8_top,T(i)).'*(F8);
disp('46%')
FG =  diff(anch_1_top,G(i)).'*(F1)...
    +diff(anch_2_top,G(i)).'*(F2)...
    +diff(anch_3_top,G(i)).'*(F3)...
    +diff(anch_4_top,G(i)).'*(F4)...
    +diff(anch_5_top,G(i)).'*(F5)...
    +diff(anch_6_top,G(i)).'*(F6)...
    +diff(anch_7_top,G(i)).'*(F7)...
    +diff(anch_8_top,G(i)).'*(F8);
disp('47%')
FP =  diff(anch_1_top,P(i)).'*(F1)...
    +diff(anch_2_top,P(i)).'*(F2)...
    +diff(anch_3_top,P(i)).'*(F3)...
    +diff(anch_4_top,P(i)).'*(F4)...
    +diff(anch_5_top,P(i)).'*(F5)...
    +diff(anch_6_top,P(i)).'*(F6)...
    +diff(anch_7_top,P(i)).'*(F7)...
    +diff(anch_8_top,P(i)).'*(F8);
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
[D2x(i), D2y(i), D2z(i), D2T(i), D2G(i), D2P(i)] = solve(Fx==fx(i), Fy==fy(i),Fz == fz(i), FT==fT(i),...
                                   FG==fG(i),FP==fP(i),d2x(i),d2y(i),d2z(i),d2T(i),d2G(i),d2P(i));

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
Dyn_eqn = [D2x; D2y; D2z; D2T; D2G; D2P];
matlabFunction(Dyn_eqn(2:2:12),'file','duct_accel','Vars',[x(2); y(2); z(2); T(2); G(2); P(2); dx(2); dy(2); dz(2); dT(2); dG(2); dP(2); tensions]);
matlabFunction(lengths,'file','lengths','Vars',[x(2); y(2); z(2); T(2); G(2); P(2)]);
matlabFunction(dlengths_dt,'file','dlengths_dt','Vars',[x(2); y(2); z(2); T(2); G(2); P(2); dx(2); dy(2); dz(2); dT(2); dG(2); dP(2)]);
disp('100%')

