%% Spine Parameters
clear
close all

% Geometric parameters
ll = .15; % m, length of long bars
h = .15; % m, height from top to bottom of tetra
ls = h/2; % m, length of short bar
w = sqrt(ll^2-(h/2)^2); % m, width from center of tetra
minCableTension = 0;


% Mass and force parameters
g = 9.8; % m/s^2, acceleration due to gravity
m = .136; % kg/tetra

% xi = [0,.1,0]; % Forget what this does.

%     1  2     3     4     5    6      7       8       9       10
x = [ 0  w    -w     0     0    0      w      -w       0        0]';
% x = [ 0  w     -w    0    0    0+.01  w+.01  -w+.01   0+.01   0+.01]';
y = [ 0  0     0     w    -w   0      0       0       w       -w]';
% y = [ 0  0     0     w    -w   .01      .01      .01      .01+w      .01-w]';
z = [ 0  -h/2  -h/2  h/2  h/2  .1     .1-h/2  .1-h/2  .1+h/2  .1+h/2]';

forcesZ = [-m*g -m*g]';
momentsX = [0 0];
momentsY = [0 0];
momentsZ = [0 0];
coms = [1 6]';
fixed = [2 3 4]';


%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."

% Full connectivity matrix
% Rows 1-4 are cables
% Rows 5-10 are bars
% Columns 1-4 are bottom tetra nodes
% Columns 5-8 are top tetra nodes
%    1  2  3  4  5  6  7  8  9 10
C = [0  1  0  0  0  0 -1  0  0  0;  %  1
    0  0  1  0  0  0  0 -1  0  0;  %  2
    0  0  0  1  0  0 -1  0  0  0;  %  3
    0  0  0  1  0  0  0 -1  0  0;  %  4
    0  0  0  1  0  0  0  0 -1  0;  %  5
    0  0  0  0  1  0 -1  0  0  0;  %  6
    0  0  0  0  1  0  0 -1  0  0;  %  7
    0  0  0  0  1  0  0  0  0 -1;  %  8
    1 -1  0  0  0  0  0  0  0  0;  %  9
    1  0 -1  0  0  0  0  0  0  0;  % 10
    1  0  0 -1  0  0  0  0  0  0;  % 11
    1  0  0  0 -1  0  0  0  0  0;  % 12
    0  0  0  0  0  1 -1  0  0  0;  % 13
    0  0  0  0  0  1  0 -1  0  0;  % 14
    0  0  0  0  0  1  0  0 -1  0;  % 15
    0  0  0  0  0  1  0  0  0 -1]; % 16


[ q, A, p ] = InvKin( C , x, y, z, forcesZ, momentsX, momentsY, momentsZ, coms, fixed, minCableTension  );


%% Sample force output check

load('forceOutputv5.mat')

bodies = length(coms);

BodyForceReader(q,A,p)