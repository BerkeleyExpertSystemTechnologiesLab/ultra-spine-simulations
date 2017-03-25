clear
close all

load('datav4.mat')

% Geometric parameters
minCableTension = 0;

% Mass and force parameters
g = gravity; % m/s^2, acceleration due to gravity
m = m_rod;

forcesZ = [-m*g -m*g -m*g -m*g -m*g -m*g -m_payload*g]';
momentsX = [0 0 0 0 0 0 0]';
momentsY = [0 0 0 0 0 0 0]';
momentsZ = [0 0 0 0 0 0 0]';
coms = [15 16 17 18 19 20 21]';

fixed = [1 5 9]';
% Bodies (1,3,6)

[q,A,p] = InvKin(C,x,y,z,forcesZ,momentsX,momentsY,momentsZ,coms,fixed,minCableTension);


%% Sample force output check

load('forceOutputv5.mat')

bodies = length(coms);

BodyForceReader(q,A,p)

