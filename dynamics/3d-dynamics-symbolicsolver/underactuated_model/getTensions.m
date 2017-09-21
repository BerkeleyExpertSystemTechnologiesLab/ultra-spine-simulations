function [ Te ] = getTensions(x2,y2,z2,T2,G2,P2,dx2,dy2,dz2,dT2,dG2,dP2,x3,y3,z3,T3,G3,P3,dx3,dy3,dz3,dT3,dG3,dP3,...
    x4,y4,z4,T4,G4,P4,dx4,dy4,dz4,dT4,dG4,dP4,inp1,inp2,inp3,inp4,inp5,inp6,inp7,inp8,rest1,rest2,rest3,rest4,rest5,rest6,rest7,rest8,link)
% getTensions.m
% Copyright 2016-2017 Jeff Friesen, Andrew Sabelhaus

% This function calculates the forces due to each of the cables in our tensegrity spine robot 
% as a function of the system states (vertebrae positions) and the rest lengths of each cable.

% The spring constant K and damping constant c are hard-coded here:
K = 2000;
c = -100;

% Do a little bit of correction on the positions and rotations of each vertebra
x2 = x2 + (abs(x2)<1e-8)*1e-6;
y2 = y2 + (abs(y2)<1e-8)*1e-6;
z2 = z2 + (abs(z2)<1e-8)*1e-6;
T2 = T2 + (abs(T2)<1e-8)*1e-6;
G2 = G2 + (abs(G2)<1e-8)*1e-6;
P2 = P2 + (abs(P2)<1e-8)*1e-6;

x3 = x3 + (abs(x3)<1e-8)*1e-6;
y3 = y3 + (abs(y3)<1e-8)*1e-6;
z3 = z3 + (abs(z3)<1e-8)*1e-6;
T3 = T3 + (abs(T3)<1e-8)*1e-6;
G3 = G3 + (abs(G3)<1e-8)*1e-6;
P3 = P3 + (abs(P3)<1e-8)*1e-6;

% To calculate the cable forces, we need the lengths and velocities of the cable
% (the vectors of the cable positions and velocities)
L = lengths(x2,y2,z2,T2,G2,P2,x3,y3,z3,T3,G3,P3,x4,y4,z4,T4,G4,P4);
dlengths = dlengths_dt(x2,y2,z2,T2,G2,P2,dx2,dy2,dz2,dT2,dG2,dP2,x3,y3,z3,T3,G3,P3,dx3,dy3,dz3,dT3,dG3,dP3,x4,y4,z4,T4,G4,P4,dx4,dy4,dz4,dT4,dG4,dP4);

% reshape a bit
L = L(:, link);
dlengths = dlengths(:, link);

% The 'stretch length' is like delta x in F = K delta x
stretch_length = [L(1) - rest1 + inp1;
                  L(2) - rest2 + inp2;
                  L(3) - rest3 + inp3;
                  L(4) - rest4 + inp4;
                  L(5) - rest5 + inp5;
                  L(6) - rest6 + inp6;
                  L(7) - rest7 + inp7;
                  L(8) - rest8 + inp8];
              
% Calculate the tension force (I think this is a scalar?) F = K delta X - C x_dot
Te = (stretch_length>=0).*(K*stretch_length-c*dlengths);

% Correct for slack cables. Cables cannot apply a negative force (no pushing.)
Te = (Te>=0).*Te;    
end

