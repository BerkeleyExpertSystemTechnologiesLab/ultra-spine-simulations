
C = [...
    0 1 0 0 -1 0;...
    0 0 1 0 0 -1;...
    1 -1 0 0 0 0;...
    1 0 -1 0 0 0;...
    0 0 0 1 -1 0;...
    0 0 0 1 0 -1];

x = [0 1.5 -1.5 -.5 1.5 -1.5]';
z = [0 0 0 -1 -1 -1]';
y = [0 0 0 0 0 0]';
forcesZ = [-1 -1]';
momentsX = [0 0]';
momentsY = [0 2]';
momentsZ = [0 0]';
coms = [1 4]';
fixed = [2 3]'; 
minCableTension = 0;


[q,A,p] = InvKin( C, x, y, z, forcesZ, momentsX, momentsY, momentsZ, coms, fixed, minCableTension  );

BodyForceReader(q,A,p)
