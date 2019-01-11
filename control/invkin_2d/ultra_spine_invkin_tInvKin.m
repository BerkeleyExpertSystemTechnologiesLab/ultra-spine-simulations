%% ultra_spine_invkin_tInvKin.m
% Copyright Andrew P. Sabelhaus 2019

% This script used the tInvKin libraries to calculate the inverse
% kinematics for the 2D, two-vertebra, vertical spine that's part of the
% MPC and MPC w/InvKin tests. 
% It is different than the original (2017/2018) results in that the new
% tInvKin library is used, which has been verified in hardware as of Dec.
% 2018, so it's *known* that the inverse kinematics calculations themselves
% are theoretically correct.

% **** NOTE that the tInvKin library MUST be in the current MATLAB path /
% must be added to the path before calling this function.
% This function specifically calls (from that library):
% 
% get_node_coordinates_2d
% get_reaction_forces_2d
% invkin_core_2d_rb

function [f_opt, u_opt] = ultra_spine_invkin_tInvKin(xi, spineParameters, q_min)

%% Set up the parameters

% Debugging level.
% 0 = no output except for errors
% 1 = Starting message, results from quadprog
% 2 = Verbose output of status.
debugging = 0;

% If appropriate, output a starting message.
if debugging >= 1
    disp('Starting ultra_spine_invkin_tInvKin for 2D single-vert vertical spine...');
end

% minimum cable force density, q_min = 0.5, is minCableTension

% Parameters that are hard-coded (here is where to adapt the script to be
% more general)

s = 4; % cables

% number of rigid bodies
%b = 2;
% When removing the anchor nodes, it's like removing one of the bodies:
b = 1;

% Parameters from the passed-in struct
g = spineParameters.g; % m/s^2, acceleration due to gravity
M = spineParameters.total_m; % kg/vertebra
kappa = spineParameters.k_vert; % structure has vertical and horizontal k, but they're the same, so ignore for now
m_node = spineParameters.m; % kg/node for all nodes within a vertebra (this case, 4.)
% Thus, change it to one node:
m_node = m_node(1); % To-do: check this!!!
% For the local frame of each vertebra,
a = spineParameters.a;
% For creating the system state,
N = spineParameters.N;

if debugging >= 2
    a
end

% Configuration matrix for WHOLE STRUCTURE.

% Full connectivity matrix
% Rows 1-4 are cables
% Rows 5-10 are bars
% Columns 1-4 are bottom tetra nodes
% Columns 5-8 are top tetra nodes
%    1  2  3  4  5  6  7  8  
C = [0  1  0  0  0 -1  0  0;  %  1
     0  0  1  0  0  0 -1  0;  %  2
     0  0  0  1  0 -1  0  0;  %  3
     0  0  0  1  0  0 -1  0;  %  4
     1 -1  0  0  0  0  0  0;  %  5
     1  0 -1  0  0  0  0  0;  %  6
     1  0  0 -1  0  0  0  0;  %  7
     0  0  0  0  1 -1  0  0;  %  8
     0  0  0  0  1  0 -1  0;  %  9
     0  0  0  0  1  0  0 -1]; % 10

% r follows directly, since we've hard-coded s.
r = size(C,1) - s;
% ...because C is \in R^{10 x 8}.

% number of nodes
n = size(C, 2);

if debugging >= 2
    C
end

% Taking from the struct and putting into usable form 
% (TO-DO: how does the struct's specification align with the dynamics??)
m = ones(n, 1) * m_node;

% Example of how to do the 'anchored' analysis.
% Declare a vector w \in R^n, 
% where w(i) == 1 if the node should be 'kept'.
% For this example, want to treat body 1 as the anchored nodes.
% So, we zero-out anchored nodes 1 through 4, and keep nodes 5-8
% (which is vertebra two.)
w = [0; 0; 0; 0; 1; 1; 1; 1];
% Including all nodes:
%w = ones(n,1);

% IMPORTANT! If chosing to remove nodes, must change 'b' also, or else inv
% kin will FAIL.

% We also need to declare which nodes are pinned (with external reaction
% forces) and which are not.
% We're choosing not to have this be the same as w, since there are some
% nodes to "ignore" for w which are not necessarily built in to the ground
% for "pinned". Example, nodes inside the leftmost vertebra, where we're
% deciding to assume that only the tips of its "Y" are supported.

% For the tetrahedral vertical spine, nodes 2 and 3 (pointy parts of Y) are
% pinned.
pinned = zeros(n,1);
pinned(2) = 1;
pinned(3) = 1;

%% Position/state

% all the positions of each rigid body (expressed as their COM positions
% and rotation). that's 3 states: [x; z; \gamma] with
% the angle being an intrinsic rotation.

% The xi that's passed in here is the rigid body state for the vertebra
% (we're doing dynamics now.) We don't want that, we want it to be 3 states
% one for each vertebra including the rigidly-fixed one.

xi_all = zeros(3*N, 1);

% Insert the moving vertebra's states (non-moving is zeros state)
% make indexing better here I'm lazy for now
xi_all(4:6) = xi(1:3);

% if debugging >= 2
%     xi
% end

%% Calculations for the inputs to the core invkin library

% The nodal coordinates (x, z)
% calculate from state
% initialize the results. There are n nodes (n position vectors.)
x = zeros(n,1);
y = zeros(n,1);

% We split up according to rigid body. The fixed body is xi_all(1:3,i)
coordinates_i_fixed = get_node_coordinates_2d(a, ...
    xi_all(1:3), debugging);
% and for the free vertebra
coordinates_i_free = get_node_coordinates_2d(a, ...
    xi_all(4:6), debugging);
% At this point along the trajectory, get the coordinates
% ...and split them into coordinate-wise vectors per node, per
% vertebra.
x(1:4) = coordinates_i_fixed(1,:)';
x(5:8) = coordinates_i_free(1,:)';
y(1:4) = coordinates_i_fixed(2,:)';
y(5:8) = coordinates_i_free(2,:)';    

if debugging >= 2
    %coordinates
    x
    y
end

% Reaction forces can be calculated by this helper, which assumes that only
% gravity is present as an external force.
% Initialize results
% px = zeros(n);
% py = zeros(n);

[px, py] = get_reaction_forces_2d(x, y, pinned, ...
    m, g, debugging);

% for more details, you can look at commits to the library before Nov. 2018
% where this reaction force/moment balance was written out by hand.

% Since this was just a calculation of the reaction forces, we ALSO need to
% add in the external forces themselves (grav forces) for use as the whole
% inverse kinematics problem.

% Add the gravitational reaction forces for each mass.
% a slight abuse of MATLAB's notation: this is vector addition, no indices
% needed, since py and m are \in R^n.
py = py + -m*g;


%% Solve the inverse kinematics problem

% Solve, over each point.
% Let's use a cell array for Ab and pb, since I don't feel like thinking
% over their sizes right now.
%f_opt = zeros(s);
%q_opt = zeros(s);
% we also need the lengths for calculating u from q*.
%lengths = zeros(s);
Ab = {};
pb = {};

% finally, the big function call:
% quadprog is inside this function.
[f_opt, q_opt, lengths, Ab_i, pb_i] = invkin_core_2d_rb(x, ...
    y, px, py, w, C, s, b, q_min, debugging);
% and insert this Ab and pb.
Ab{end+1} = Ab_i;
pb{end+1} = pb_i;

% Seems correct, intuitively!
% Cable 1 is horizontal, below.
% Cable 2 is horizontal, above.
% Cable 3 is saddle, below.
% Cable 4 is saddle, above.

% A quick plot of the cable tensions.
% figure; 
% hold on;
% subplot(4,1,1)
% hold on;
% title('Cable tensions');
% plot(f_opt(1,:))
% ylabel('1 (N)');
% subplot(4,1,2)
% plot(f_opt(2,:))
% ylabel('2 (N)');
% subplot(4,1,3)
% plot(f_opt(3,:))
% ylabel('3 (N)');
% subplot(4,1,4);
% plot(f_opt(4,:));
% ylabel('4 (N)');

%% Convert the optimal forces into optimal rest lengths.
% u_i = l_i - (F_i/kappa_i)

% save in a vector
u_opt = zeros(s,1);
% TO-DO: DOES THIS WORK WITH NUMPOINTS = 1?
% it's more intuitive to iterate for now. At least, we can iterate over
% cables and not over timesteps.
for k=1:s
    % For cable k, divide the row in f_opt by kappa(k)
    % But, now include the length offset term. Accounts for the initial
    % spring length, as well as the little extender we had to use.
    %u_opt(k, :) = lengths(k,:) - init_len_offset(k) - (f_opt(k,:) ./ kappa(k));
    u_opt(k) = lengths(k) - (f_opt(k) ./ kappa);
end

% a bit of debugging, now
if debugging >= 2
    disp('Optimal rest lengths are:');
    u_opt
end















