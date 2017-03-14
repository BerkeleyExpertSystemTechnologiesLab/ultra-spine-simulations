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

[ q, A, p] = InvKin( C , x, y, z, forcesZ, momentsX, momentsY, momentsZ, coms, fixed, minCableTension );


%% Sample force output check

load('forceOutputv5.mat')

cableOuts = A*q;
disp(' ')
for a = 1:7
    bodyOuts = cableOuts(6*a-5:6*a);
    disp(['For body ' num2str(a) ' cable forces produce'])
    disp(['X Force:' num2str(-bodyOuts(1))])
    disp(['Y Force:' num2str(-bodyOuts(2))])
    disp(['Z Force:' num2str(-bodyOuts(3))])
    disp(['X Moment:' num2str(bodyOuts(4))])
    disp(['Y Moment:' num2str(bodyOuts(5))])
    disp(['Z Moment:' num2str(bodyOuts(6))])
    disp(' ')
end

