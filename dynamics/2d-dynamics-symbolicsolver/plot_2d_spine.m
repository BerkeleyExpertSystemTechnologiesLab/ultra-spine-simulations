function handles = plot_2d_spine( xi, geometry )
%plot_2d_spine.m Plots the location of the single-vertebra tensegrity spine, in 2D
%   This function plots a very simple illustration of the 2D, "inverted-Y",
%   tensegrity spine. It assumes that the plotting window is open and ready to be 
%   drawn on. (So, for simulations, it would be good to clear the window inside a loop 
%   before calling this function - use the 'handles' that are returned.)
%   Inputs:
%       xi, the system state vector, \in R^6. This is (x, z, theta, dx, dz, dtheta).
%       geometry, a struct that has the geometric parameters of the spine
%   Outputs:
%       handles, the handles to each of the elements that was drawn this timestep.
%           These are used in the calling function to clear what was drawn this timestep.

% First, pull out the geometry of the spine.
% (See the two_d_dynamics_* scripts, where this struct is defined.)
% The "leg length" of the vertebra:
l = geometry.l;
% The total height of a vertebra:
h = geometry.h;

% Next, calculate the local coordinate frame of a vertebra.
% Again, see the dynamics derivation scripts, from which this is copied.
% These are the four node locations (transposed, for ease of writing.
a = [0,                     0;
     -(l^2 - (h/2)^2)^.5,   -h/2;
     (l^2 - (h/2)^2)^.5,    -h/2;
     0,                     h/2]';
 
% The handles array can be a cell array:
handles = {};

% Plot the bottom (not-moving) spine
% First, the four nodes, as circles:
for i=1:size(a,2);
    handles{end+1} = plot(a(1,i), a(2,i), 'k.', 'markersize', 40);
end

% Plot the vertebra itself (lines between nodes)
% t1 to t2:
handles{end+1} = line(a(1,1:2), a(2,1:2),'Color','k');
% t1 to t3:
handles{end+1} = line( [a(1,1), a(1,3)], [a(2,1), a(2,3)],'Color','k');
% t1 to t4:
handles{end+1} = line( [a(1,1), a(1,4)], [a(2,1), a(2,4)],'Color','k');

% Next, calculate the nodes of the moving vertebra.
% These are the locations a, translated by (x,z), and rotated by theta.

% The rotation matrix for this given angle of theta is:
rot = [ cos(xi(3)),    -sin(xi(3));
        sin(xi(3)),     cos(xi(3))];
    
% We need to multiply each point in a by rot,
% and add the (x,z) offset from the xi state vector.
% To make this easier, so that we can use matrix algebra,
% copy out the (x,z) states to a 2x4 matrix:
xz_offset = repmat( [xi(1); xi(2)], 1, 4);

% Now, the locations of the nodes of the moving vertebra are
b = rot*a + xz_offset;

% Finally, plot these in the same way as the static vertebra,
% only with a different color.
% First, the four nodes, as circles:
for i=1:size(b,2);
    handles{end+1} = plot(b(1,i), b(2,i), 'b.', 'markersize', 40);
end

% Plot the vertebra itself (lines between nodes)
% t1 to t2:
handles{end+1} = line(b(1,1:2), b(2,1:2),'Color','b');
% t1 to t3:
handles{end+1} = line( [b(1,1), b(1,3)], [b(2,1), b(2,3)],'Color','b');
% t1 to t4:
handles{end+1} = line( [b(1,1), b(1,4)], [b(2,1), b(2,4)],'Color','b');

% Finally, plot the cables, between each node.
% 1) vertical, a2 to b2:
handles{end+1} = line( [a(1,2), b(1,2)], [a(2,2), b(2,2)],'Color','r');
% 2) vertical, a3 to b3:
handles{end+1} = line( [a(1,3), b(1,3)], [a(2,3), b(2,3)],'Color','r');
% 3) saddle, a4 to b2:
handles{end+1} = line( [a(1,4), b(1,2)], [a(2,4), b(2,2)],'Color','r');
% 4) saddle, a4 to b3:
handles{end+1} = line( [a(1,4), b(1,3)], [a(2,4), b(2,3)],'Color','r');

end

