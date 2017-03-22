%% Schek Example
% 12/1/16
clear; clc

% 2D system
%      *1
%      |
% 4*---*---*2
%      |
%      *3
% Center node is 5
% Bars are numbered 1-4 starting at top and rotating clockwise
% Node 4 is fixed

% Connectivity matrices

    %    1  2  3  4    
    C = [1  0  0  0;  % 1
         0  1  0  0;  % 2
         0  0  1  0;  % 3
         0  0  0  1]; % 4
     
    %      5  
    Cf = [-1;  % 1
          -1;  % 2
          -1;  % 3
          -1]; % 4
      
    Cs = [C Cf];

% Free nodes

    x1 = sym('x1','real');
    x2 = sym('x2','real');
    x3 = sym('x3','real');
    x4 = sym('x4','real');
    x = [x1 x2 x3 x4]';

    y1 = sym('y1','real');
    y2 = sym('y2','real');
    y3 = sym('y3','real');
    y4 = sym('y4','real');
    y = [y1 y2 y3 y4]';

    z1 = sym('z1','real');
    z2 = sym('z2','real');
    z3 = sym('z3','real');
    z4 = sym('z4','real');
    z = [z1 z2 z3 z4]';

% Fixed nodes    
xf = 1;
yf = 0;
zf = 0;

% Coordinate differences of connected points
u = C*x + Cf*xf;
v = C*y + Cf*yf;
w = C*z + Cf*zf;

% Diagonal matrices of coordinate differences
U = diag(u);
V = diag(v);
W = diag(w);

% Lengths of bars
l = [1 1 1 1]';

% Length matrix
L = diag(l);

% Forces in bars
s1 = sym('s1','real');
s2 = sym('s2','real');
s3 = sym('s3','real');
s4 = sym('s4','real');
s = [s1 s2 s3 s4]';

% Equilibrium equations
px = C'*U*inv(L)*s;
py = C'*V*inv(L)*s;
pz = C'*W*inv(L)*s;

% Define q
q = inv(L)*s;
Q = diag(q);

% Redefine equlibrium equations
px = C'*U*q;
py = C'*V*q;
pz = C'*W*q;