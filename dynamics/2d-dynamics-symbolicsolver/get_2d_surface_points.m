% Andrew P. Sabelhaus, Feb 2018
% heavily adapted from Jeff Friesen's cylinder2P.

function [x, y, z] = get_2d_surface_points(radius, rad_disc, ...
                        length_disc, start_pt, end_pt)
%get_2d_surface_points Calculate points for a surf of a bar for a 2D
%tensegrity
%
%   This function calculates and returns a set of x, y, z points that can
%   be passed into 'surf' to create a cylinder. In this case, it's used for
%   the bars of a tensegrity, in plot_2d_tensegrity_surfaced.
%
% Inputs:
%   radius = radius of the cylinder
%   rad_disc = discretization along the radius (number of points)
%   length_disc = discretization long the long edge of the cylinder
%   start_pt = a single 3D point in space, the center of the circle at one
%       end of the cylinder
%   end_pt = point at other end of cylinder.
%
% Outputs:
%   x, y, z = points to pass into surf.

% Start a grid of the points along the circle and cylinder lengths.
theta = linspace(0, 2*pi, rad_disc);
dl = linspace(0, 1, length_disc); % since we're moving along a unit vector

% Preallocate the whole array of points.
% There are going to be length_disc x rad_disc number of points, and we're
% storing them in separate vectors since that's easier for surf.
x = zeros(length_disc, rad_disc);
y = zeros(length_disc, rad_disc);
z = zeros(length_disc, rad_disc);

% rename for convenience and brevity:
r1 = start_pt;
r2 = end_pt;

% Unit vector along the axis of the cylinder
v = (r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
%cylinder axis described by: r(dl)=r1+v*dl for 0<dl<1

% Let's get two arbitrary vectors that complete a set of axes for the
% cylinder.
% We can get a random linearly independent vector (of v) to seed the
% computation
R2 = rand(1,3);
% Get an orthogonal vector to v (the second axis)
x2 = v-R2/(R2*v');
% Normalize (make the orthogonal -> orthonormal)
x2 = x2/sqrt(x2*x2');
% Get the third axis, orthonormal to both v and x2.
x3 = cross(v,x2);
% renormalize it too.
x3 = x3/sqrt(x3*x3');

% Pull out each element component-wise (makes the computations a bit less
% messy, less indexing:
% start and end points
r1x = r1(1); r1y = r1(2); r1z = r1(3);
r2x = r2(1); r2y = r2(2); r2z = r2(3);
% axes that direct around the circular edge of the cylinder
x2x = x2(1); x2y = x2(2); x2z = x2(3);
x3x = x3(1); x3y = x3(2); x3z = x3(3);

% Finally, grid out the points along the cylinder surface.
for i = 1 : length_disc
    % pull out the length along cylinder from the length discretization grid
    l = dl(i);
    % for the i-th circle, store all points around the circle.
    x(i, :) = r1x + (r2x-r1x)*l + radius*cos(theta)*x2x + radius*sin(theta)*x3x; 
    y(i, :) = r1y + (r2y-r1y)*l + radius*cos(theta)*x2y + radius*sin(theta)*x3y; 
    z(i, :) = r1z + (r2z-r1z)*l + radius*cos(theta)*x2z + radius*sin(theta)*x3z;
end


end

