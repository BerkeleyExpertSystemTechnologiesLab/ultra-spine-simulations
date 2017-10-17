% InvKin_Horizontal_Spine.m
% Copyright Andrew P. Sabelhaus and the Berkeley Emergent Space
% Tensegrities Lab, 2017

% This script sets up an inverse kinematics problem for a horizontal,
% tetreahedral, tensegrity spine.
% Then, the inverse kinematics solver is called, and the results are
% plotted.

% See InvKin.m for the inverse kinematics solver itself.

% As of the current commit, this script doesn't provide a correct answer,
% but does at least give something. Should allow for some calculations of
% the forces the spine will experience, plus a safety factor.

%% Setup: Spine Parameters
clear
close all
clc

% Geometric parameters

% % The non-symmetric spine was:
% ll = .15; % m, length of long bars
% h = .15; % m, height from top to bottom of tetra
% ls = h/2; % m, length of short bar
% w = sqrt(ll^2-(h/2)^2); % m, width from center of tetra

% The symmetric spine is:
h = 0.15; % m, height from bottom to top of vertebra
w = 0.075; % m, the distance from the center to one edge

% Used to constrain all cables to have some pre-tension
minCableTension = 0.5; % Newtons, I think?

% Mass and force parameters
g = 9.8; % m/s^2, acceleration due to gravity
m = .136; % kg/tetra

% We could specify each of the positions of each node manually:
% %     1  2     3     4     5    6      7       8       9       10
% x = [ 0  w    -w     0     0    0      w      -w       0        0]';
% y = [ 0  0     0     w    -w   0      0       0       w       -w]';
% z = [ 0  -h/2  -h/2  h/2  h/2  .1     .1-h/2  .1-h/2  .1+h/2  .1+h/2]';

% For the 2D spine (I think?) would be:
% x = [ 0  w     -w    0    0    0+.01  w+.01  -w+.01   0+.01   0+.01]';
% y = [ 0  0     0     w    -w   .01      .01      .01      .01+w      .01-w]';

% Or, we could specify one local coordinate frame:
% (all in meters)
x_local = [ 0  w    -w     0     0]';
y_local = [ 0  0     0     w    -w ]';
z_local = [ 0  -h/2  -h/2  h/2  h/2]';

% Let's make this easier, we can use vectors to represent coordinates.
% Each point coordinate is a column vector.
coordinates_local = [x_local'; y_local'; z_local'];

% The above coordinates are correct for a vertical spine.
% However, for the horizontal spine, the vertebra needs to be rotated
% by 90 degrees around the x-axis.
% Noting that sin(pi/2) = 1, cos(pi/2) = 0;
% The rotation matrix here would be:
Rot = [1, 0, 0; ...
       0, 0, -1; ...
       0, 1, 0];
   
% Apply the rotation:
coordinates_local = Rot * coordinates_local;

% Note: by Drew's understanding, after doing the rotation, the nodes are
% the following:
% 1 = center
% 2, 3 = middle
% 4 = top (+z)
% 5 = bottom (-z)

% ...and then translate the local frame to make global frames for each
% body.

% Translation between each body:
%translation = [0, 0, 0.1]'; % translate by 0.1m in the Z-direction.
translation = [0, 0.1, 0]'; % translate by 0.1m in the Y-direction. This is the horizontal spine.

% Expand this vector out so we can add it to each body.
% We can determine the number of nodes based on how many coordinates have
% been specified above (e.g., 5 columns means 5 nodes.)
num_nodes = length(x_local); 
% So the matrix that represents a translation of each nodal coordinate is
translation = repmat(translation, 1, num_nodes);
% ...which we can add to each local frame. 

% Initialize each position to have the first body at the origin
coordinates = coordinates_local;

% Append new bodies' nodal coordinates.
% including the base body:
num_bodies = 2;
% example: with two bodies, run this loop once.
for i=1:(num_bodies-1) 
    % Append new coordinates for nodes.
    % First, the new frame will be 'i' translations away from the base
    % body.
    % We can multiply the scalar i by the matrix 'translation'.
    coordinates_translated = coordinates_local + i * translation;
    % Then append the new frame.
    coordinates = [coordinates, coordinates_translated];
    % ...where we keep each point as a row vector.
end

% Finally, the inverse kinematics function requires separate vectors of
% each of the x, y, z positions.
% Looks like they need to be column vectors.
x = coordinates(1,:)';
y = coordinates(2,:)';
z = coordinates(3,:)';

% To-Do: these vectors need to be made to be the same dimension as the
% number of bodies.
forcesZ = [-m*g -m*g]';
momentsX = [0 0];
momentsY = [0 0];
momentsZ = [0 0];
coms = [1 6]';

% The vertical spine had three points fixed???? This needs to be checked,
% does NOT seem correct.
fixed = [2 3 4]';

% For the horizontal spine, let's fix the center and the two middle
% nodes. This is kind of like supporting the spine by its 'center'.
%fixed = [1 2 3]';
% ...doesn't seem to work. maybe doesn't provide a reaction force in one
% direction?



%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."

% Full connectivity matrix
% Rows 1-8 are cables.
% Rows 9-12 are rods for v1, 13-16 are rods for v2
% Columns 1-5 are bottom tetra nodes
% Columns 6-10 are top tetra nodes

% For the horizontal spine:
%    1  2  3  4  5  6  7  8  9 10
C = [0  1  0  0  0  0 -1  0  0  0;  %  1 horiz
     0  0  1  0  0  0  0 -1  0  0;  %  2 horiz
     0  0  0  1  0  0  0  0 -1  0;  %  3 horiz
     0  0  0  0  1  0  0  0  0  -1; %  4 horiz
    
     0  1  0  0  0  0  0  0 -1  0;  %  5 saddle
     0  0  1  0  0  0  0  0 -1  0;  %  6 saddle
     0  1  0  0  0  0  0  0  0  -1; %  7 saddle
     0  0  1  0  0  0  0  0  0  -1; %  8 saddle
    
    1 -1  0  0  0  0  0  0  0  0;  %  9
    1  0 -1  0  0  0  0  0  0  0;  % 10
    1  0  0 -1  0  0  0  0  0  0;  % 11
    1  0  0  0 -1  0  0  0  0  0;  % 12
    
    0  0  0  0  0  1 -1  0  0  0;  % 13
    0  0  0  0  0  1  0 -1  0  0;  % 14
    0  0  0  0  0  1  0  0 -1  0;  % 15
    0  0  0  0  0  1  0  0  0 -1]; % 16

[ q, A, p, tensions ] = InvKin( C , x, y, z, forcesZ, momentsX, momentsY, momentsZ, coms, fixed, minCableTension  );


%% Sample force output check

BodyForceReader(q,A,p)

disp('');
disp('The cable tensions are, in Newtons,')
tensions
