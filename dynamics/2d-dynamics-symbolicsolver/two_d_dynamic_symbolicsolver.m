%% two_d_dynamics_symbolicsolver.m

% Berkeley Emergent Space Tensegrities Lab and Andrew P. Sabelhaus
% Copyright 2016
% This script symbolically solves Lagrange's equation for 2-dimensional
%   tensegrity systems which involve repeated "units," each of which are
%   a rigid body that is defined by point masses at nodes, and where all the
%   point masses are rigidly connected. Then, connection points between different
%   "units" are defined, as pairs of nodes that have a cable between them.
%   Note that I might use the word "node" to talk about the location of a point mass.
%   This is inspired by, and some of the terminology comes from, 
%   the type-II tensegrity spine under development in the BEST Lab.
% Many thanks to Jeff Friesen and Abishek Akella for earlier versions of this script.
% PLEASE NOTE, VERY IMPORTANT:
%   The dynamics created by this script assume that all tensions are nonnegative.
%   If, for example, the units to the system cause some tensions to be negative,
%   these dynamics will allow that cable to "push."
%   So, if these are used, you MUST include constraints on the tensions such that
%   the dynamics always hold. 
%   You can do this through the 'tensions' symbolic variable that's calculated
%   below, which gives a column vector of cable tensions as a function of
%   the system state xi. This can be treated as some weird, nonconvex state constraint,
%   something like tensions(xi) >= 0.

% This script generates multiple m-files:
%   xi_dot.m
%   tensions.m
%   dlengths_dt.m
%   lengths.m

% TO-DO: validate the use of the whole state vector xi for fulldiff.
% This works if there are no velocity terms in the equations that we're differentiating (which
% is true for the positions), but MAY NOT BE TRUE for the LHS of Lagrange's equations
% when we take the derivative of that partial.

% Huajing Zhao, 6/20/2017
% NOTE: 
% Error occurs at line (roughly) 1169, while N > 2:
% d2xi_solved_un(i) = replace_derivatives(getfield(d2xi_solved,
% char(d2xi(i))), xi, num_states_per_unit, debugging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% While N = 2, the code works:
%
% d2xi_solved_un: 3*1 sym;
% run: replace_derivatives(getfield(d2xi_solved, char(d2xi(i))), ...)
% Replacing derivatives...
% ans = (500*(24*tensions_un3*cos(xi3 + pi/6) + 45*tensions_un1*sin(xi3) + ...
% size(ans) = 1 1;    
% class(ans) = double
%
% getfield(d2xi_solved, char(d2xi(i)))
% ans = (500*(24*tensions_un3*cos(xi3 + pi/6) + 45*tensions_un1*sin(xi3) + ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% While N = 3, error occurs: 
% "In an assignment  A(:) = B, the number of elements in A and B must be
% the same";
%
% d2xi_solved_un: 6*1 sym;
%
% run: replace_derivatives(getfield(d2xi_solved, char(d2xi(i))), ...)
% Replacing derivatives...
% ans = Empty sym: 0-by-1
% getfield(d2xi_solved, char(d2xi(i))):
% ans = Empty sym: 0-by-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% While N = 4, the same error occurs:
% d2xi_solved_un: 1*1 sym;
% replace_derivatives(getfield(d2xi_solved, char(d2xi(i))), ...
% ans = Empty sym: 0-by-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE in line 1062, about solving Lagr_eqns and write into d2xi_solved
% d2xi_solved = solve( Lagr_eqns, accel_vars);
%% 1) Prepare the workspace:
clc;
clear all;
close all;

% A flag to turn debugging on and off.
% This is useful to see the symbolic variables that are created.
debugging = 1;

% Some paths to the files that will be saved.
% These are the names of the functions and files that are generated at the
% end of this script. Change these to whatever you want!
% They are all relative paths with respect to the current folder.
two_d_geometry_path = 'two_d_geometry.mat';

% As of 2016-11-13, fulldiff uses some functionality
% that apparently will be deprecated in a future release of MATLAB.
% For now, turn off that warning.
% TO-DO: re-write fulldiff!
warning off symbolic:sym:DeprecateFindsym

% We'll be calling simplify quite often.
% Specify the number of steps that the simplifier should take:
num_simplify_steps = 100;

% Throughout this script, I output some messages that show the progress
% of the script during its calculations.
% Those are labelled 'PROGRESS_BAR'.

%PROGRESS_BAR
disp('*********************************');
disp('Two-Dimensional Tensegrity Dynamics Symbolic Solver');
disp(' ');

% Some other parameters to look at:
% the 'equal_masses' flag in the following section,
% ...

%% 2) Geometric and physical parameters of the 2D tensegrity system

%PROGRESS_BAR
disp('Defining tensegrity system physical parameters...');

% This script seeks to automatically derive the dynamics
% for an arbitrary 2D tensegrity.
% It's assumed that this system consists of N rigid bodies,
% each of which are translated in some direction in relation
% to each other.
% We do this by first defining the local coordinate system
% for one unit, using a series of point masses which are
% assumed to be rigidly connected.
% Then, for N-1 of these structures, symbolically calculate
% the velocities, kinetic and potential energies, the 
% Lagrangian, and then the left-hand-side of Lagrange's equations.
% After that, calculate the right-hand-side (the cable forces),
% and call 'solve' to get xi_dot.

% This script assumes that the first unit is fixed at the origin,
% and that the first node of this first unit is at exactly (0,0,0).
% There are N units, including the fixed one:
% N = 3;
N= 2;

% Gravitational constant (used for calculating potential energy in Lagrangian):
g = 9.8;

% Next, we define the positions of each of the point masses
% for a single unit's local coordinate system.

% For this example, one unit is an upside-down Y-shape,
% for the 2D "flemons spine" example. Let's define that
% Y-shape by the length of each of the three lines in the Y,
% and then the total height of one Y.
% We've been calling these lengths "legs", so continue that notation here.
% For others: feel free to go ahead and hard-code your coordinates 
% of point masses in the a matrix below, we're just being fancy here.
leg = 0.15;
h = 0.15;
% for this specific design, the "width" of one unit is then
w = 2 * sqrt(leg^2 - (h/2)^2);

% Each mass position is a column vector,
% but let's define it in transpose here, since
% it's easier to write out a vector as a row
% in MATLAB.
% Again, a1 is the (arbitrary) center of the unit, and
% these are positions in (x,z).
% HERE IS WHERE YOU WOULD CHANGE THIS CODE TO MATCH YOUR TENSEGRITY MODEL.
a = [ 0,        0; ...
      -w/2,     -h/2; ...
      w/2,      -h/2; ...
      0,        h/2]';

% Here's a lop-sided upward-C-shape design:
% An interval in the x-direction for determining
% the other distances:
x_interval = 0.05;
% Same thing, a z-interval:
z_interval = 0.05;
% The node positions are then:
% a = [ 0,            0; ...
%       -x_interval,  z_interval; ...
%       x_interval,   z_interval; ...
%       2*x_interval, 0; ...
%       4*x_interval, 2*z_interval]';
  
% again, note the ' transpose above.
% Given the a matrix above, we can pick out how many
% point masses there are per unit:
num_pm_unit = size(a,2);
% ... again noting that a is transposed.

% This script, itself, does not depend on the connections between
% point masses within a unit. They could be connected in any way,
% it doesn't matter. HOWEVER, for plotting purposes,
% it's useful to define which node-node connections should be
% plotted.
% Call these "bars" since the tensegrity unit will consist
% of bars that connect nodes. 
% Note that this matrix must be symmetric:
% If point 2 is connected to point 3, then point 3
% is connected to point 2 by definition.
% THUS, BY CONVENTION, use the lower triangle
% Enforce this by putting "NaN"
% in the places that we won't be using.
% Note also that there should be no "1"s along the
% main diagonal: nodes don't connect to themselves.
% The "NaN" function generates matrices of NaN,
% and the "triu" generates upper triangular matrices.
bars = triu(NaN(num_pm_unit, num_pm_unit));
%bars = zeros(num_pm_unit, num_pm_unit);
% Insert a "1" where two nodes are connected.
% For our specific example, each of the outer nodes
% is connected to the center node, and that's it.
bars(2,1) = 1;
bars(3,1) = 1;
bars(4,1) = 1;

% For the C-shape:
% Bars connect nodes 2,1; 1,3; 3,4; 4,5;
% bars(2,1) = 1;
% bars(3,1) = 1;
% bars(4,3) = 1;
% bars(5,4) = 1;

% Similarly, define the mass of each point mass.
% This could be done either of the following ways, you pick,
% depending on if the mass in each unit is evenly distributed
% (each has the same mass) or if they are unevenly distriuted
% and you would like to assign each of them individually.
% Here's a handy flag to switch between versions:
equal_masses = 1;
% The masses of each point will be stored in a vector:
m = zeros(num_pm_unit, 1);
% First, a version where all the masses are equal
if equal_masses
    % Let's say each unit has a certain total mass:
    unit_total_mass = 0.1360;
    % Then, the mass of each point will be:
    m_each = unit_total_mass / num_pm_unit;
    % Put into our vector of masses:
    m(:) = m_each;
else
    % Assign each mass individually.
    % NOTE that you have to MANUALLY check to be sure
    % you assign the proper number of point masses!
    m = [ 0.35; 0.35; 0.4; 0.5];
    % A quick automatic check to be sure that m
    % has the required dimensions:
    error_msg = 'Your vector of masses m is not the same size as the number of point masses!!';
    assert( (size(m,1) == num_pm_unit) && (size(m,2) == 1), error_msg);
end

% Next, define the cables that connect adjacent units.
% These are modeled as spring-damper systems, F = k \delta x - c dx/dt.
% In order to determine the lengths of these cables, 
% let's create a "connectivity matrix" that represents
% which point mass locations are connected via cables to
% the locations on adjacent units.
% The following relation is assumed to hold for each 
% pair of units in order, e.g., between 1 and 2, between 2 and 3, up to N.

% The row index is the "from" node, and the column index is the "to" node.
% At each location where there is a cable, specify the spring constant
% and damping constant for that cable.
% So, 'connections' will be a cell array, containing vectors of doubles.
% Note that you can assign the constants however you'd like, including writing
% out each one, but for the inverted-Y-spine example, 
% let's declare some sets of constants for ease of use.
% We'll have 'vertical' cables and 'saddle' cables, as defined by me.
k_vert = 2000;
k_saddle = 2000;
% NOTE THAT these damping constants should be POSITIVE.
c_vert = 50;
c_saddle = 50;

connections = cell(num_pm_unit, num_pm_unit);
% NOTE that these are assuming that a 'lower' unit is the 'from', 
% and it is connecting to a 'to' that's 'above' it.
% For example, with the inverted-Y spine vertebra, node 4 from a lower
% vertebra connects to node 3 of the vertebra above it, 
% but it's not true that node 3 from a lower vertebra connects to 
% node 4 of the vertebra above it.
% For our specific example, the following nodes have cable connections
% between them,
% with 2 vertical and 2 saddle cables:
connections{2,2} = [k_vert, c_vert];
connections{3,3} = [k_vert, c_vert];
connections{4,2} = [k_saddle, c_saddle];
connections{4,3} = [k_saddle, c_saddle];

% For the c-shape:
% Connections between 3,3 (vertical)
% Connections between 2,1; 3,1; 3,4; and 4,5 (saddle)
% connections{3,3} = [k_vert, c_vert];
% connections{2,1} = [k_saddle, c_saddle];
% connections{3,1} = [k_saddle, c_saddle];
% connections{3,4} = [k_saddle, c_saddle];
% connections{5,4} = [k_saddle, c_saddle];

% For later below, calculate how many cables we'll expect to have
% in this system.
% Now that 'connections' is a cell matrix, count up the nonzeros by 
% checking 'isempty' on each element. The ~ negates 'isempty',
% since we want to count non-empty elements.
% So, connections_locations is a matrix of integers (0 or 1).
connections_locations = ~cellfun('isempty', connections);
% Then, count the nonzeros.
num_cables_per_unit = nnz(connections_locations);
% Remembering that there is one set of cables between each two units,
% the total number for the whole tensegrity system is then
num_cables = num_cables_per_unit*(N-1);
% Think about it this way: there is no set of cables reaching "upward"
% from the final unit.

% Save all of this geometry information in a struct, for plotting later.
two_d_geometry.N = N;
two_d_geometry.g = g; % even though this isn't geometry, it's useful to save.
two_d_geometry.a = a;
two_d_geometry.num_pm_unit = num_pm_unit;
two_d_geometry.bars = bars;
two_d_geometry.connections = connections;
two_d_geometry.connections_locations = connections_locations;
two_d_geometry.num_cables_per_unit = num_cables_per_unit;
two_d_geometry.num_cables = num_cables;
% Add the spring constants for the cables.
two_d_geometry.k_vert = k_vert;
two_d_geometry.k_saddle = k_saddle;
two_d_geometry.c_vert = c_vert;
two_d_geometry.c_saddle = c_saddle;
% Include the mass of the system.
two_d_geometry.m = m;
two_d_geometry.total_m = sum(m);
% For backwards compatibility: if the leg length and height variables
% are declared, save them too.
% On 2018-06-05: removed these variables just to be sure that the 'a' frame
% is the only variable used in the MPC and plotting. E.g., we should just
% be able to declare 'a' here, and have it used everywhere, without having
% to fiddle with widths and heights.
% if exist('leg','var')
%     % I had called this "l" in the past.
%     two_d_geometry.l = leg;
% end
% if exist('h','var')
%     two_d_geometry.h = h;
% end

% For the logistic barrier function, declare two constants that control
% its form and location.
% See below for more information about what these are and what they do.
% Take a look at the Wikipedia page about analytic approximations to 
% the step function:
% https://en.wikipedia.org/wiki/Heaviside_step_function#Analytic_approximations
logistic_k = 500;
logistic_x0 = 0.01;

%% 3) Create the symbolic variables for the solver to use

%PROGRESS_BAR
disp('Creating symbolic variables...');

% The size of the system state varies with number of units.
% For each unit, there are six state variables: [x; z; theta];
% and their derivatives [dx; dz; dtheta];,
% which are the location of the origin of the coordinate system,
% and the rotation about that center of the coordinate system,
% for each unit.
% Call this state vector xi. (This is the greek letter \xi.)

% Noting that the first unit is fixed at the origin,
% the number of state variables is
num_states_per_unit = 6;
num_states = (N-1)*num_states_per_unit;

% Add this to the geometry struct, too.
two_d_geometry.num_states_per_unit = num_states_per_unit;
two_d_geometry.num_states = num_states;

% NOTE that we constrain these symbolic variables
% to be in the real numbers. Complex-number distance vectors
% don't make sense here.
% The state vector is then
xi = sym('xi', [num_states, 1], 'real');

% We'll be solving for the derviative of the state vector,
% as in the equation xi_dot = f(xi, u), so we need to declare
% that as a symbolic variable also.
% This will be our final result at the end of the script.
%xi_dot = sym('xi_dot', [num_states, 1], 'real');
% In order to get the fulldiff results we want,
% need to declare the derivatives of xi and second derivatives
% of xi to both be real numbers.
dxi = sym('dxi', [num_states, 1]);
d2xi = sym('d2xi', [num_states, 1]);
assume(dxi, 'real');
assume(d2xi, 'real');

% We will also need the position vectors for each point mass,
% for each unit.
% Let's make a big 3D matrix r, where the first two dimensions
% are the (x,z) position of the point mass, and the final dimension
% is the unit to which that point mass belongs.
% This matrix has the same shape as a.
r = sym('r', [2, num_pm_unit, N]);

% In addition to the position vectors, we'll need to store the 
% velocities of each point mass in order to calculate the
% kinetic energy of each mass.
% These are of the same dimension as the position vectors.
r_dot = sym('r_dot', [2, num_pm_unit, N]);

% We'll be storing the Lagrangian, L = T-V, for each unit,
% as a symbolic expression in terms of xi.
% This will be a scalar for each unit, so store it as
% a column vector.
% TO-DO: should we be storing this for the 1st unit? L=0 always, there...
Lagrangian = sym('Lagrangian', [N, 1]);

% When we equate the left-hand side and right-hand side of Lagrange's equations,
% to solve for xi_dot, we'll be equating some functions of each of the 
% position coordinates in xi. E.g, functions of (x,z,theta) for all units.
% Create a set of symbolic variables to store these functions.
% This first one represents d/dt * (partial L / partial xi_dot), where xi_dot
% are the velocities of the units, expressed in terms of the point mass velocities.

% NOTE THAT THE xi_dot SYMBOLIC VARIABLE IS NOT USED HERE. Since we have not solved
% for xi_dot explicitly yet (that's the whole point of this script!), we'll be 
% differentiating with respect to some element of xi itself. See the ACC2017 ULTRA Spine
% paper for eqn (6), which makes this make more sense.

% We need one of these for each of the position variables, or half the state vector:
ddt_L_xi_dot = sym('ddt_L_xi_dot', [size(xi,1)/2, 1]);
% We also need the other left-hand-side term, (partial L / partial xi),
L_xi = sym('L_xi', [size(xi,1)/2, 1]);

% Now, for the cables:
% Create a symbolic variable for the length of each cable.
% Store them all in a single vector:
lengths = sym('lengths', [num_cables, 1]);
% We'll also need to store the rate of change of the cable lengths,
% in order to calculate the damping force.
dlengths_dt = sym('dlengths_dt', [num_cables, 1]);
% Finally, we'll store the tensions due to each cable,
% which are the forces that the cables enact on the rigid bodies.
tensions = sym('tensions', [num_cables, 1]);

% As of 2016-11-19: This didn't work. The MATLAB solve/simplify commands
% don't work if 'tensions' is full of the actual tensions on the cables.
% I guess that makes some sense, maybe the sines and cosines make things cancel.
% So, instead, let's keep a completely separate set of tensions,
% not solving for the tensions as a function of system states,
% but then substitute one for the other later once Lagrange's equations are solved.
% Call it "tensions_un" for tensions unsolved.
tensions_un = sym('tensions_un', size(tensions));

% Each cable will have not just a tension, but also the two locations
% where the tension acts. 
% The first column is the "from" location, and the second
% is the "to" location. 
% Additionally, we can store the two units that a cable acts on
% This could be calculated using some smart indexing, but it's
% easier to do it this way.
% So, the third column here are the "from" and "to" indices
% of a specific cable, thus a 2x3 matrix.
% An example would be:
% [ rx_from, rz_from; rx_to, rz_to; 2, 3]'
% ...for a cable connecting units 2 and 3, with the from node
% located on unit 2, and the to node on unit 3.
tension_points = sym('tension_points', [2, 3, num_cables]);

% When solving for the right-hand side of Lagrange's equations,
% we'll need to store the (global) forces on each unit.
% These are (partial r \ partial q_i)'*F_cable(i), for a force acting
% at location r.
% These will each be one scalar.
% Store each force vector as a column.
% Note that we don't care about the forces on the first unit, so only store N-1 of these.
global_forces = sym('global_forces', [num_states_per_unit/2, N-1]);
% Since we'll be adding to these forces as cables are iterated over,
% initialze all to be zero.
global_forces(:) = 0;

% This system has inputs, too.
% Here, let's have the inputs represent the rest length in the cables,
% e.g., u is the new rest length of the spring-cable mechanism,
% transformed into cable tensions below.
% There is one rest length per cable.
u = sym('u', [num_cables, 1]);

%% TO-DO: do these dynamics need to change now that our spine is not symmetric in 3D?

% Answer: I think not. If we were calculating the rigid body dynamics,
% not using point masses, then yes we'd have to calculate the center of mass
% and express the rotational velocity/acceleration around that COM. But here,
% we're just using point masses. The rotation here is one coordinate frame
% versus another, which can be completely arbitrary, as along as it's consistent.

%% 4) Next, constrain the first unit. 

%PROGRESS_BAR
disp('Adding constraints on the first unit...');

% The positions of its nodes are equal to a in both the
% local and global coordinate frames, since this unit does not move.
r(:,:,1) = a;

% The velocities of each point mass here are exactly zero.
% Make a matrix of zeros that fits exactly.
r_dot(:,:,1) = zeros(size(r_dot(:,:,1)));

%% 5) Then, express the coordinates of each point mass in terms of the system states.

%PROGRESS_BAR
disp('Assigning the point mass locations in terms of system states...');

% We'll need to store rotation matrices for each moving unit,
% each of which will be of size 2x2.
% TO-DO: preallocate this. Should it be of type sym?
%R = zeros(2,2,N-1);

% For each of the 2...N units that are moving,
% the point masses are rotated by the 3rd coordinate in that unit's block
% of the state vector (the angle, theta) and translated by the 1st and 2nd
% coordinates (the (x,z) vector.)
for k=1:N-1
    %PROGRESS_BAR
    disp(strcat('     Assigning point mass locations for unit number: ', num2str(k+1)));
    % At the k-th unit, calculate the index of the angle theta
    % into the state vector, noting that the state vector starts at unit
    % two, since the first unit does not move.
    % Here, k=1 is unit 2, and so on.
    % NOTE: this make assumptions about the number of variables (3 pos, 3 velocity.)
    x_index = 1 + (k-1)*num_states_per_unit;
    z_index = 2 + (k-1)*num_states_per_unit;
    theta_index = 3 + (k-1)*num_states_per_unit;
    % Calculate the rotation matrix for this unit,
    % with respect to the global coordinate system.
    R(:,:,k) = [cos(xi(theta_index)),  -sin(xi(theta_index)); ...
                sin(xi(theta_index)),   cos(xi(theta_index))];
    % Then, the position of each point mass will be
    % that mass's position in the local frame, rotated by R(k),
    % translated by the amount between local/global ref frame.
    % Iterate over each point mass:
    for p=1:num_pm_unit
        % The pm-th point mass is a column vector.
        % Note that the index k is with respect to the moving units,
        % so it's really the k+1th unit in terms of the point mass locations r.
        % r = R*a + e; where e is the position of the origin of the coordinate frame for 
        % this unit, e.g., the value of the state vector for x and z for this unit.
        r(:,p,k+1) = R(:,:,k)*a(:,p) + xi(x_index:z_index);
    end
end

%% 6) Similarly, express the velocities of each point mass in terms of xi.

%PROGRESS_BAR
disp('Assigning point mass velocities in terms of system states...');

% For each of the moving units...
for k=1:N-1
    %PROGRESS_BAR
    disp(strcat('     Calculating symbolic derivatives of point masses for unit number: ', num2str(k+1)));
    % For each of the point masses in this unit:
    for p=1:num_pm_unit
        % The velocity is the full derivative of position (with respect to time.)
        % Credit goes to Tim Jorris for the fulldiff function.
        % Note, however, that the independent variables must be passed in as a cell array.
        r_dot(:,p,k+1) = fulldiff(r(:,p,k+1), sym2cell(xi));
    end
end

% Now, to make the simplification easier on the solver later, substitute
% all the derivatives that fulldiff calculated back into the symbolic variables xi.
% Fulldiff writes these variables as, for example, dxi1, dxi1, etc., which are really
% xi4, xi5, etc.
% Note that this works, since the velocities of each point mass will not be
% dependent on the velocity of any other point mass: e.g., there will not be
% any second derivatives in r_dot.

%PROGRESS_BAR
%disp('Substituting system states back into r_dot...');
% Swap out all the 'dxi(whatever)' for the xi+3 coordinate
%r_dot = replace_derivatives(r_dot, xi, num_states_per_unit, debugging);

%% 7) Calculate the kinetic and potential energy, and the Lagrangians, for the whole system.

%PROGRESS_BAR
disp('Solving for L = T-V for each unit...');

% Use the same indexing as in the previous section.
% Note that Lagrangian(1) doesn't matter, it's never used.
for k=1:N
    % For this unit, calculate the kinetic energy, (1/2)mv^2,
    % as in T = (1/2)* sum( m(i)*r_dot(i)^2), for point masses i.
    % Add in an extra simplify step to make things go faster, later.
    %T = sym;
    syms T;
    for p=1:num_pm_unit
        % Add the kinetic energy for this point mass,
        % noting that the vector m has the masses for each point mass.
        % Here, we use r_dot'*r_dot as the 2-norm.
        T = T + (1/2) * m(p)*(r_dot(:,p,k)'*r_dot(:,p,k));
    end
    % Similarly, calculate the potential energy for this unit.
    % The potential energy is gravity times the z component of each unit.
    % Since the point masses may have different mass, do a loop.
    %V = sym;
    syms V;
    for p=1:num_pm_unit
        % Add the potential energy of this point mass,
        % noting that the z component is the second row of r.
        V = V + m(p)*g*r(2,p,k);
    end
        
    % simplify both of these.
    T = simplify(T, num_simplify_steps);
    V = simplify(V, num_simplify_steps);
    
    % Calculate the Lagrangian, T-V;
    Lagrangian(k) = T-V;
    % Simplify this too.
    Lagrangian(k) = simplify(Lagrangian(k), num_simplify_steps);
    
    %DEBUGGING
    if debugging
        disp(strcat('Kinetic energy for unit: ', num2str(k), ' is'));
        disp(T);
        disp(strcat('Potential energy for unit: ', num2str(k), ' is'));
        disp(V);
        disp(strcat('Lagrangian for unit: ', num2str(k), ' is'));
        disp(Lagrangian(k));
    end
end

%% 8) Finally, we can calculate the left-hand side of Lagrange's equations of motion.

% AS OF COMMIT 093d381, this is consistent with Jeff's work.
% The bug was that I had been differentiating against Lagrangian(k+1)
% where it needed to be ddt_L_xi_dot, and also that since we're using
% dxi now instead of replacing derivatives, the variable of differentiation 
% had to be updated.

%PROGRESS_BAR
disp('Calculating the left-hand-side of Lagranges equations...');

% Lagrange's equation(s) of motion are of the following form:
% (d/dt) * (partial L)/(partial xi_dot) - (partial L)/(partial xi)
% ==
% sum(forces acting on the masses due to the cables).
% Here, the 'xi_dot' term is really just the velocity terms
% inside xi, but I wrote it a bit lazily.

% Let's also start up some parallel pools here for quicker calculation.
%pools = gcp;

% Though we know that there are only 6 variables per unit here, let's still use 
% the n variable to calculate which states are positions and
% which are velocities, so this could be used for 2D and 3D dynamics in the future.
% This variable should be 3: (when n=6)
velocity_start_offset = num_states_per_unit/2;

% Let's store a counter into the symbolic variables, which will index
% into the ddt_L_xi_dot and L_xi variables, which will are the ones
% which will really be solved later to get xi_dot.
count = 1;

% For each unit,
for k=1:N-1    
    % We need to calculate some indices into the state vector xi.
    % This is needed to determine how to take the derivatives for Lagrange's eqns.
    % The state variables for this unit start at intervals of num_states_per_unit apart,
    % and end at the next interval of num_states_per unit.
    % For example, in the 6-state-per-unit spine vertebra, these
    % intervals are 1-6, 7-12, 13-18, ...
    unit_index_start = 1 + (k-1)*num_states_per_unit;
    velocity_index_start = unit_index_start + velocity_start_offset;
    unit_index_end = (k)*num_states_per_unit;
    % Calculate the left-hand-side equations for this unit.
    % Iterate over the velocity coordinates for this unit:
    % (e.g., 4-6, 10-12, ...)
    for p = velocity_index_start:unit_index_end
        %PROGRESS_BAR
        disp(strcat('Calculating the LHS of Lagr. Eqns. for state pair: ', ...
            char(xi(p-velocity_start_offset)), ', ', char(xi(p))));
        % The first term is the full time derivative of (partial L / partial xi_dot).
        % Note that I'm not storing L_xi_dot as a symbolic variable itself,
        % it gets overwritten here. That's just because it is never used directly,
        % we only calculate it in order to calcualate ddt_L_xi_dot.
        % Note that we must index into 'Lagrangrian' via k+1 and not k,
        % since Lagrangian(1) is the constant value for the not-moving unit.
        
        % WHEN THIS IS FIXED, and I use replace_derivatives again
        % We'll differentiate with respect to the correct position in xi,
        % like this:
        %L_xi_dot = diff(Lagrangian(k+1), xi(p));
        % FOR NOW: we need to use a dxi variable.
        L_xi_dot = diff(Lagrangian(k+1), dxi(p - velocity_start_offset));
        
        % Then take the full time derivative:
        ddt_L_xi_dot(count) = fulldiff( L_xi_dot, sym2cell(xi));
        % Replace the dxi* terms in the symbolic variable:
        %ddt_L_xi_dot(count) = replace_derivatives(ddt_L_xi_dot(count), xi, num_states_per_unit, debugging);
        
        % The second term is derivatives with respect to the POSITION variables,
        % so subtract away the +3 offset. 
        % My apologies for all this index-shuffling.
        L_xi(count) = diff(Lagrangian(k+1), xi( p - velocity_start_offset));
        
        % Do a quick simplify. Prior work used a parallel pool here,
        % might as well do that again.
        % ...
        
        count = count+1;
    end
end

% Now, the LHS = ddt_L_xi_dot - L_xi. Done!
% At this point, 'count' should equal
% (num_states)/2.

% In comparison with Jeff's script, where he uses (for example)
% fx(2) to represent the LHS of Lagrange's equations for the x-position
% of vertebra, 2
% fx(2) == ddt_L_xi_dot(1) - L_xi(1)
% fx(3) == ddt_L_xi_dot(4) - L_xi(4)
% ...note that this implies that we could have indexed directly
%    into ddt_L_xi_dot and L_xi, instead of using the 'count'
%    variable.
% This comparison can be seem from his line 491, where he 
% calls D2x(i) = solve(Fx(i) == fx(i)).
% The dimensions of fx(i) are 1x1.

%% 9) Now, start calculating the cable forces.
% Calculate the lengths of each cable, as well as the rate of change
% of those lengths, so that we can calculate F = kx - c dx/dt.

%PROGRESS_BAR
disp('Calculating cable tensions...')

% An important assumption:
% Below, the cable "speeds", the rate-of-length-change of the cables,
% which are scalar quantities, 
% are calculated as a function of the system states. This assumes
% that the rate of length change of the cables can be related directly
% to the velocity of the tensegrity units (e.g., vertebrae).
% This isn't necessarily clear to Drew right now: can we come up 
% with a counterexample where these dlength/dt quantities don't relate
% to the vertebrae velocities in the way we assume? For example, does an
% angular velocity between two vertebrae get captured properly in these cable "speeds"?
% Need to think more about if this is valid or not...

% The length of each cable is the distance between
% the two point masses specified by the 'connections' matrix.
% Remember that we've already calculated the location
% of each of these point masses: they are in the 'r' matrix.
% Keep a counter into the symbolic variables for the cables:
cable_num = 1;
% Let's iterate over the units, noting that the top unit
% does not have any cables above it.
for i=1:N-1
    % The 'from' and 'to' unit for all cables with this pairing
    % are i and i+1.
    from_unit = i;
    to_unit = i+1;
    % Then, let's iterate over the connections matrix,
    % calculating lengths if we find a 1.
    % Iterate row-wise.
    % TO-DO: should these be indexed according to point masses
    % or according to cables? I believe the connections_locations matrix
    % is of size (num_pm_unit)x(num_pm_unit).
    %for k=1:num_cables_per_unit
    for k=1:num_pm_unit
        for p=1:num_pm_unit
            % If there is a cable between these indices:
            % (note that we can use the connections_locations 
            %  matrix here, since that's already calculated,
            %  instead of having to check ~isempty on connections.)
            if connections_locations(k,p) == 1
                %PROGRESS_BAR
                disp(strcat('     Calculating symbolic length of cable number: ', num2str(cable_num)));
                % Then pick out the two locations of the nodes
                from_node = r(:,k,from_unit);
                to_node = r(:,p,to_unit);
                % Finally, the length is the 2-norm
                % of the difference between these points.
                lengths(cable_num) = norm(to_node-from_node, 2);
                % Do a quick simplify step
                lengths(cable_num) = simplify(lengths(cable_num), num_simplify_steps);
                
                % Similarly, we can calculate the change in cable lengths.
                %PROGRESS_BAR
                disp(strcat('     Calculating symbolic dlengths_dt of cable number: ', num2str(cable_num)));
                % Calculate by calling fulldiff again.
                dlengths_dt(cable_num) = fulldiff( lengths(cable_num), sym2cell(xi));
                % Replace out the derivatives and simplify:
                %dlengths_dt(cable_num) = simplify( ...
                %    replace_derivatives(dlengths_dt(cable_num), xi, num_states_per_unit, debugging), num_simplify_steps);
                dlengths_dt(cable_num) = simplify( dlengths_dt(cable_num), num_simplify_steps);

                % Next, calculate the (scalar) tensions in each of these cables.
                %PROGRESS_BAR
                disp(strcat('     Calculating symbolic tension of cable number: ', num2str(cable_num)));
                % The input to the system, u, is the "rest length" of each cable,
                % e.g., emulating a motor retracting a cable attached to a spring.
                
                % NOTE THAT right here is where the "cable pushing" problem
                % occurs. It's fully possible for the calculated tension here to be
                % "negative," making these dynamics equations no longer valid.
                % Remember that the force in a cable is max( tension, 0).
                % But, since including a "max" would make the resulting symbolic equations
                % quite horrible, let's do the constraining later on.
                % TO-DO: maybe make this 'smooth' by doing some type of log-barrier method,
                % as is used in convex optimization, to make the tension go arbitrarily close to 0
                % as the cables become 'slack.'
                
                % F = k \delta x - c * d/dt (lengths)
                % NOTE THAT it seems that dlengths_dt really contains -d/dt(lengths), 
                % so it's a "+" here.
                tensions(cable_num) = ...
                    connections{k,p}(1) * (lengths(cable_num) - u(cable_num)) ...
                    + connections{k,p}(2) * dlengths_dt(cable_num);
                
                % Do a quick simplify step to give the symbolic solver
                % an easier time later.
                tensions(cable_num) = simplify(tensions(cable_num), num_simplify_steps);
                
                % Finally, save the two locations where this cable applies its force.
                tension_points(:,1,cable_num) = from_node;
                tension_points(:,2,cable_num) = to_node;
                % We're also saving the indices of the units that this cable attaches to.
                % These are scalars.
                tension_points(1,3,cable_num) = from_unit;
                tension_points(2,3,cable_num) = to_unit;
                
                % Increment the counter into the cables matrices.
                cable_num = cable_num+1;
            end
        end
    end
end

%% 10) Calculate the forces on each unit.

%PROGRESS_BAR
disp('Calculating global cable forces...');

% For each unit, calculate (in global coordinates)
% the forces that act on it.
% These are (partial r \ partial q_i)'*F_cable(i), for a cable force acting
% at location r, where F is 2x1 and (partial r \ partial q_i) is also 2x1,
% since r is 2x1.

% This corresponds to lines starting at 370 in Jeff's code.

% Perform this iteration over the units.
% Since tension_points is indexed to include 1 as the non-moving unit,
% let's keep with that indexing.
% (sorry, I know this is different than above.)
for i=2:N
    % Iterate over our list of cables, and if a cable
    % attaches to this unit...
    for k=1:num_cables
        % If this cable is a "from" for this unit:
        % (note that tension_points is a symbolic matrix, need
        %  to convert back to integers/doubles.)
        if double(tension_points(1,3,k)) == i
            %PROGRESS_BAR
            disp(strcat('Calculating forces due to cable ', num2str(k), ' on unit ', num2str(i), '...'));
            disp('(note, this is a "from" cable.)');
            % Then calculate the force vector associated with 
            % cable k, noting that this cable will pull this unit
            % "upward" since it's the 'from'.
            % This force is tension * direction,
            % where direction is (to - from).
            % NOTE that since the solver doesn't work when plugging in the actual
            % tension directly, let's use the symbolic tension for now, then substitute
            % back later below.
            %F_cable = tensions(k) * (tension_points(:,2,k) - tension_points(:,1,k));
            F_cable = tensions_un(k) * (tension_points(:,2,k) - tension_points(:,1,k));
            % Run a simplify step
            F_cable = simplify(F_cable, num_simplify_steps);
            % Next, calculate (partial r \ partial q_i) for the
            % 'from' node location = r.
            % The O'Reilly textbook calls these covariant basis vectors.
            % We need to do one per each of the three coordinates.
            for q=1:num_states_per_unit/2
                % ...note that we could just as easily do q=1:3,
                % so each of the (x, z, theta) alignment is implied here.
                % We need to differentiate with respect to the
                % corresponding xi state:
                % Shift up q according to which unit is under consideration at the moment:
                q_shifted = q + num_states_per_unit*(i-2);
                % This will give, for example, 1-3 and 7-9.
                % This cable attaches at the 'from' point, column 1 in tension_points.
                cov_bas_vec = diff( tension_points(:,1,k), xi(q_shifted));
                % Run a simplify step
                cov_bas_vec = simplify(cov_bas_vec, num_simplify_steps);
                % FINALLY, ADD TO THE TOTAL GLOBAL FORCE IN THE q-DIRECTION
                % THIS IS THE ACTUAL, FINAL COMPUTATION FOR THE RIGHT-HAND SIDE
                % OF LAGRANGE'S EQUATIONS!!!
                % (note the transpose: this is a scalar.)
                force_addition = simplify(cov_bas_vec'*F_cable, num_simplify_steps);
                %global_forces(q, i-1) = global_forces(q, i-1) + cov_bas_vec'*F_cable;
                global_forces(q, i-1) = global_forces(q, i-1) + force_addition;
                %DEBUGGING
                if debugging
                    disp(strcat('Contribution to global force for unit ', num2str(i), ' in direction ', num2str(q), ' was: ', char( force_addition )));
                end
            end
        % Now, also check if this cable is, instead, a "to" for this unit:
        elseif double(tension_points(2,3,k)) == i
            %PROGRESS_BAR
            disp(strcat('Calculating forces due to cable ', num2str(k), ' on unit ', num2str(i), '...'));
            disp('(note, this is a "to" cable.)');
            % Then calculate the force vector associated with 
            % cable k, noting that this cable will pull this unit
            % "downward" since it's the 'to'.
            % This force is tension * direction,
            % where direction is (from - to).
            % NOTE that since the solver doesn't work when plugging in the actual
            % tension directly, let's use the symbolic tension for now, then substitute
            % back later below.
            %F_cable = tensions(k) * (tension_points(:,1,k) - tension_points(:,2,k));
            F_cable = tensions_un(k) * (tension_points(:,1,k) - tension_points(:,2,k));
            % Run a simplify step
            F_cable = simplify(F_cable, num_simplify_steps);
            % Next, calculate (partial r \ partial q_i) for the
            % 'from' node location = r.
            % The O'Reilly textbook calls these covariant basis vectors.
            % We need to do one per each of the three coordinates.
            for q=1:num_states_per_unit/2
                % ...note that we could just as easily do q=1:3,
                % so each of the (x, z, theta) alignment is implied here.
                % We need to differentiate with respect to the
                % corresponding xi state:
                % Shift up q according to which unit is under consideration at the moment:
                q_shifted = q + num_states_per_unit*(i-2);
                % This will give, for example, 1-3 and 7-9.
                % Note that the cable attaches at the 'to' point, column 2 in tension_points.
                cov_bas_vec = diff( tension_points(:,2,k), xi(q_shifted));
                % Run a simplify step
                cov_bas_vec = simplify(cov_bas_vec, num_simplify_steps);
                % FINALLY, ADD TO THE TOTAL GLOBAL FORCE IN THE q-DIRECTION
                % THIS IS THE ACTUAL, FINAL COMPUTATION FOR THE RIGHT-HAND SIDE
                % OF LAGRANGE'S EQUATIONS!!!
                % (note the transpose: this is a scalar.)
                force_addition = simplify(cov_bas_vec'*F_cable, num_simplify_steps);
                %global_forces(q, i-1) = global_forces(q, i-1) + cov_bas_vec'*F_cable;
                global_forces(q, i-1) = global_forces(q, i-1) + force_addition;
                %DEBUGGING
                if debugging
                    disp(strcat('Contribution to global force for unit ', num2str(i), ' in direction ', num2str(q), ' was: ', char( force_addition )));
                end
            end
        end
    end
end

%% 11) Simplify the global forces as much as possible.

%PROGRESS_BAR
disp('Calling simplify on the global cable forces, in parallel...');
% create the parallel pools. We'll be needing them starting now.
pools = gcp;
% Each simplify has one output, and will take our input
% alongside a set number of simplify steps.
% We'll keep a cell array of all the pools:
% Note that N is the total number of units, including the not-moving unit.
% There are 6 states per unit, and 3 "directions" of (x,z,theta).
global_forces_pools = cell(num_states_per_unit/2, N-1);
% Index along the number of the unit, first.
for i=1:size(global_forces,2)
    % Then, for each of the 3 directions (x,z,theta) per unit:
    for j=1:size(global_forces,1)
        % Create a pool to simplify each of the (three) global forces
        % for this unit
        global_forces_pools{j,i} = parfeval(pools, @simplify, 1, global_forces(j,i), 'Steps', num_simplify_steps);
    end
end

% Then, fetch the outputs from the pools:
disp('Fetching simplified outputs from parallel pool...');
% As above, index first by unit number, then by direction.
for i=1:size(global_forces,2)
    for j=1:size(global_forces,1)
        disp(strcat('     Fetching forces for direction:', num2str(j), ' , for unit number: ', num2str(i)));
        global_forces(j,i) = fetchOutputs(global_forces_pools{j,i});
    end
end


%% 12) Before solving, need to express xi_dot properly.

% We'll have a set of equations that has some dxi(i) in it,
% for some coordinates, but not for others.
% For example, there is no dxi1 term, since that's just equal to xi4.
% To solve, let's express xi_dot, then ask the solver 
% to get xi_dot in terms of xi.
% Iterate over units:
for i=1:N-1
    % The derivative of position coordinates at this index
    % are the velocity coordinates at this index
    unit_index_start = 1 + (i-1)*num_states_per_unit;
    velocity_index_start = unit_index_start + velocity_start_offset;
    unit_index_end = (i)*num_states_per_unit;
    xi_dot(unit_index_start: (velocity_index_start-1)) = xi(velocity_index_start:unit_index_end);
    % Then, an easy way to create the dxi(i) variables
    % is to call fulldiff. This saves us re-naming everything by hand.
    for k=velocity_index_start:unit_index_end
        xi_dot(k) = fulldiff(xi(k), xi(k));
    end
end
% xi_dot should be a column vector.
assume(xi_dot, 'real');
xi_dot = xi_dot';

%% 13) Equate the LHS and RHS of Lagrange's equations and solve!

disp('SOLVING LAGRANGES EQUATIONS...');

% Note that at this point, the tensions are still
% independent variables: we'll plug the actual tensions, as functions
% of the system state, back into the solution once it's computed.

% Create a symbolic array of all the equations that will be solved.
% There will be (number of directions)*(number of moving units) equations.
% For example, with 2 moving units, that's 3*2 = 6 equations.
Lagr_eqns = sym('Lagr_eqns', [(num_states_per_unit/2)*(N-1), 1], 'real');
% Loop through and assign each equation:
for i=1:length(Lagr_eqns)
    % Unfortunately, I've used different indexing for the LHS and RHS of
    % Lagrange's equations. The LHS are in a list, and the RHS are in
    % a (num_states/2) x (N-1) array. For example, 3x1 for 2 tensegrity units
    % (one static, the other moving.)
    % This was because we iterated over the cables for the calculation of the RHS,
    % not the units as in the LHS.
    % However, the list is ordered according to unit, so we can say:
    % (where "direction" is x,z, or theta, which is 1 to 3):
    direction = mod(i-1,num_states_per_unit/2) + 1;
    % This gives, for example: i=2, direction=2, i=8, direction=2, etc.
    % Similarly, calculate the unit number:
    unit_num = ceil(i / (num_states/2));
    % This gives, for example: i=2, unit=1, i=8, unit=3, etc.
    Lagr_eqns(i) = (ddt_L_xi_dot(i) - L_xi(i) == global_forces(direction, unit_num));
    % These look something like:
    % ddt_L_xi_dot(1) - L_xi(1) == global_forces(1,1), ...
    % ddt_L_xi_dot(2) - L_xi(2) == global_forces(2,1), ...
    % ddt_L_xi_dot(3) - L_xi(3) == global_forces(3,1), ...
    % ddt_L_xi_dot(4) - L_xi(4) == global_forces(1,2), ...
    % ddt_L_xi_dot(5) - L_xi(5) == global_forces(2,2), ...
    % ddt_L_xi_dot(6) - L_xi(6) == global_forces(3,2), ...
end

% Finally, solve all these equations simultaneously,
% with the independent variables being the accelerations.
% Make a big vector of the accelerations, noting that
% we get one acceleration result per Lagrange's equation:
accel_vars = sym('accel_vars',size(Lagr_eqns), 'real');
for i=1:length(accel_vars)
    % The acceleration variables are d2xi(1), d2xi(2), d2xi(3), d2xi(6), d2xi(7),...
    % Indexing goes:
    % 4 to 7
    % 5 to 8
    % 6 to 9
    % 7 to 13
    % 8 to 14
    % 9 to 15
    unit_num = ceil(i / (num_states/2));
    accel_index = i + (unit_num-1)*(num_states/2);
    accel_vars(i) = d2xi(accel_index);
end

% SOLVE LAGRANGE'S EQUATIONS! This is the culmination of all the code
% in this script!
d2xi_solved = solve( Lagr_eqns, accel_vars);

%% 14) Create a version of the solved tensions that includes a barrier-function constraint

% One possible way to solve the "negative tensions" problem
% would be to constrain all the cable forces to be nonnegative, manually
% changing any negative tensions to 0 if they are negative.
% Another approach would be to multiply the tension functions by
% some approximation to the step function: if tension is less than 0,
% then this step function makes the tension zero.

% A smooth approximation to the unit-step function is the logistic function:
% https://en.wikipedia.org/wiki/Logistic_function
% This is of the form: f(x) = 1 / (1 + exp( -k*(x-x0))),
% where k is the steepness of the curve,
% and x0 is the midpoint of the curve (which is effectively
% an offset, laterally, from x=0.)
% These constants are declared above as k_logistic, x0_logistic.

% First, do a replace_derivatives step on the tensions, so they don't
% have to be modified later:
tensions_sub = replace_derivatives(tensions, xi, num_states_per_unit, debugging);

% A new symbolic variable for the tensions with barrier:
tensions_sub_barrier = tensions_sub;

% For each tension force, add a barrier function:
for i=1:length(tensions_sub_barrier)
    % The new logistic function looks like
    % f(x) = 1 / (1 + exp( -k*(x-x0)))
    logistic_func = 1 / (1 + ...
        exp(-logistic_k*(tensions_sub_barrier(i) - logistic_x0)));
    % The rectified tension is then the original tension
    % multiplied by this logistic function.
    tensions_sub_barrier(i) = logistic_func*tensions_sub_barrier(i);
end


%% 14) Substitute the solved tensions back into the accelerations and simplify

disp('Finally, substitute tensions back into solved accelerations, and replace with states in the xi vector...');

% This script produces four different forms of the equations of motion,
% and the end-user can simulate or use any of the three depending on 
% what the problem specifies. These four are:
% 1) Separate acceleration and tension functions. This approach
%       requires that the tensions function be called to calculate
%       the cable tensions given a system state and inputs, then these
%       tensions are used in the acceleration calculation function.
%       This approach is useful, for example, when running Model-Predictive Control
%       and trying to use a continous/smooth dynamics function but
%       constrain the tensions to be non-negative.
% 2) Separate acceleration and tension functions, with constraints
%       included for nonnegative tensions. In this case, the
%       tensions function is NOT smooth, since any nonnegative tensions
%       are replaced with 0 tension. This approach is useful for
%       simulating realistic system dynamics while checking the tensions
%       separately for another purpose (not sure what you'd want to do,
%       but here it is just in case!)
% 3) Combined xi_dot function that includes both the acceleration and tension
%       calculations, and does not include constraints on the cable tensions.
%       This function is smooth, but it allows cables to "push".
%       Use this function VERY CAREFULLY, since if the tensions are
%       negative, the simulation becomes incorrect. This approach
%       may be useful if your controller does not operate in regions
%       of the state space where tensions are negative: if so, no need
%       to constrain the cable tensions and make the problem more difficult
%       than it needs to be. ALSO, this is the same sort of dynamics function
%       that Skelton uses in his work, as well as folks working on 6-bar spherical
%       tensegrity systems, since those are fully in tension most of the time.
% 4) Combined xi_dot function that includes constraints on cable tensions.
%       As with the separate tension/accel functions with constraints, this
%       approach has non-smooth dynamics due to the constraint. This approach
%       is probably the easiest and most realistic to use when simulating the
%       dynamics of the system, since you don't have to consider constraints
%       independently, and have only one function to be called instead of two.
%       However, it does "hide" the non-smooth behavior in a sense, in that
%       there is no real way to determine if cables are slack or in tension,
%       other than to use this approach in combination with one of the above.

%%
%PROGRESS_BAR
disp('Preparing acceleration solutions for dynamics approach #1...');
% 1) Create an 'accel' function that does not substitute the tensions that were calculated.
%    Note that d2xi_solved is a column vector, same size as accel_vars from above.
%    The prefix "un_" is used here to designate "un-substituted".

% d2xi_solved_un = sym('d1xi_solved_un', size(accel_vars), 'real');
d2xi_solved_un = sym('d1xi_solved_un', size(accel_vars), 'real');
% d2xi_solved_ty = sym('d1xi_solved_un', size(accel_vars), 'real');

% The solved accelerations, without replacing tensions_un->tensions:
for i=1:(N-1)*3  %length(d2xi_solved_un) = 6
    % The output of 'solve' is a struct, so the getfield
    % function can be used to extract its elements.
    % Note that there are still dxi terms that must be converted back
    % into xi states, via the replace_derivatives function.
    % Note that getfield takes a string as the second argument, not a sym.
    
% i.e., tensions_sub = replace_derivatives(tensions, xi, num_states_per_unit, debugging);
    
 d2xi_solved_un(i) = replace_derivatives(getfield(d2xi_solved, char(d2xi(i))), ...
     xi, num_states_per_unit, debugging);

%   d2xi_solved_un(i) = replace_derivatives(getfield(d2xi_solved_ty, mat2cell(d2xi(i))), ...
%       xi, num_states_per_unit, debugging);

%    d2xi_solved_un(i) = replace_derivatives(d2xi(i), ...
%          xi, num_states_per_unit, debugging);
    
end

%PROGRESS_BAR
disp('Preparing acceleration solutions for dynamics approach #2...');
disp('     tensions_sub_barrier was created above.');

%PROGRESS_BAR
disp('Preparing acceleration solutions for dynamics approach #3...');
% 3) Create a single xi_dot: replace the tensions inside the accelerations
%    This is the same as above, but now, an additional substitution step:
d2xi_solved_sub = sym('d2xi_solved_sub', size(accel_vars), 'real');
% We'll also need to store the solved xi_dot:
xi_dot_soln = sym('xi_dot_soln',[num_states, 1], 'real');
for i=1:length(d2xi_solved_sub)
    % The output of 'solve' is a struct, so the getfield
    % function can be used to extract its elements.
    % Note that there are still dxi terms that must be converted back
    % into xi states, via the replace_derivatives function.
    % Note that getfield takes a string as the second argument, not a sym.
%     d2xi_solved_sub(i) = getfield(d2xi_solved, char(d2xi(i)));
    d2xi_solved_sub(i) = d2xi(i);
    % Then, substitute for the solved tensions:
    d2xi_solved_sub(i) = subs(d2xi_solved_sub(i), tensions_un, tensions);
    % Finally, insert into the solved xi_dot location, and perform
    % the same derivative replacement as above.
    % Indexing is:
    % 1 to 4
    % 2 to 5
    % 3 to 6
    % 4 to 10
    % 5 to 11
    % 6 to 12
    % ...
    unit_num = ceil(i / (num_states/2));
    accel_index = i + unit_num*(num_states/2);
    % Insert into the solution vector:
    xi_dot_soln(accel_index) = replace_derivatives(d2xi_solved_sub(i), ...
        xi, num_states_per_unit, debugging);
    % For the xi_dot_soln vector, the appropriate xi state
    % should also be assigned. For example, the derivative of xi(1)
    % is xi(4), etc.
    % Do that here, since we're indexing into the xi_dot_soln vector anyway.
    % This looks like:
    % xi_dot_soln(1) = xi(4);
    % xi_dot_soln(2) = xi(5);
    % xi_dot_soln(3) = xi(6);
    % xi_dot_soln(7) = xi(10);
    % ...
    % So the indexing in terms of i is:
    % i=1, velocity_index=1, accel_index=4;
    % i=4, velocity_index=7, accel_index=10;
    velocity_index = i + (unit_num-1)*(num_states/2);
    % Insert into solution vector, recalling that 'xi'
    % is still the same symbolic variable as at the top of this script
    % (nothing is assigned to xi.)
    xi_dot_soln(velocity_index) = xi(accel_index);
end

%PROGRESS_BAR
disp('Preparing accelerations solutions for dynamics approach #4...');
% 4) Create a single xi_dot: replace the tensions inside the accelerations
%    This is the same as above, but now, an additional substitution step on
%    the tensions, using the tensions with the barrier function included.
d2xi_solved_sub_barrier = sym('d2xi_solved_sub_barrier', size(accel_vars), 'real');
% We'll also need to store the solved xi_dot:
xi_dot_soln_barrier = sym('xi_dot_soln_barrier',[num_states, 1], 'real');
for i=1:length(d2xi_solved_sub_barrier)
    % The output of 'solve' is a struct, so the getfield
    % function can be used to extract its elements.
    % Note that there are still dxi terms that must be converted back
    % into xi states, via the replace_derivatives function.
    % Note that getfield takes a string as the second argument, not a sym.
    % d2xi_solved_sub_barrier(i) = getfield(d2xi_solved, char(d2xi(i)));
    d2xi_solved_sub_barrier(i) = d2xi(i);
    % Differently than above: do the replace_derivatives step now, before 
    % substitution, since tensions_sub_piecewise already has its derivatives replaced.
    d2xi_solved_sub_barrier(i) = replace_derivatives(d2xi_solved_sub_barrier(i), xi, ...
        num_states_per_unit, debugging);
    % Then, substitute for the solved tensions, using the piecewise 
    % symbolic expression from above:
    d2xi_solved_sub_barrier(i) = subs(d2xi_solved_sub_barrier(i), tensions_un, tensions_sub_barrier);
    % Finally, insert into the solved xi_dot location, and perform
    % the same derivative replacement as above.
    % Indexing is:
    % 1 to 4
    % 2 to 5
    % 3 to 6
    % 4 to 10
    % 5 to 11
    % 6 to 12
    % ...
    % TO-DO: should this be num_states_per_unit?
    unit_num = ceil(i / (num_states/2));
    accel_index = i + unit_num*(num_states/2);
    % Insert into the solution vector:
    xi_dot_soln_barrier(accel_index) = d2xi_solved_sub_barrier(i);
    % For the xi_dot_soln vector, the appropriate xi state
    % should also be assigned. For example, the derivative of xi(1)
    % is xi(4), etc.
    % Do that here, since we're indexing into the xi_dot_soln vector anyway.
    % This looks like:
    % xi_dot_soln(1) = xi(4);
    % xi_dot_soln(2) = xi(5);
    % xi_dot_soln(3) = xi(6);
    % xi_dot_soln(7) = xi(10);
    % ...
    % So the indexing in terms of i is:
    % i=1, velocity_index=1, accel_index=4;
    % i=4, velocity_index=7, accel_index=10;
    velocity_index = i + (unit_num-1)*(num_states/2);
    % Insert into solution vector, recalling that 'xi'
    % is still the same symbolic variable as at the top of this script
    % (nothing is assigned to xi.)
    xi_dot_soln_barrier(velocity_index) = xi(accel_index);
end

%% 15) Write MATLAB functions(s).

disp('Writing MATLAB functions to files...');

% For both the tensions and dlengths_dt, need to replace derivatives
% with states in xi.
dlengths_dt_sub = replace_derivatives(dlengths_dt, xi, num_states_per_unit, debugging);

% Save the geometry struct:
save(two_d_geometry_path, 'two_d_geometry');

% Need to write our lengths, dlengths_dt, tensions, 
% and finally, xi_dot_soln.
% Note that lengths and dlengths_dt are just functions of state,
% but tensions and xi_dot are functions of state and inputs,
disp('     Writing lengths function (not used, for debugging only)...');
matlabFunction(lengths,'file','two_d_spine_lengths','Vars',{[xi]});
disp('     Writing dlengths_dt function (not used, for debugging only)...');
matlabFunction(dlengths_dt_sub,'file','two_d_spine_dlengths_dt','Vars',{[xi]});

%PROGRESS_BAR
disp('     Writing functions for dynamics approach #1:');
disp('     Writing tensions function...');
matlabFunction(tensions_sub,'file','two_d_spine_tensions','Vars',{xi,u});
disp('     Writing accelerations solution, without tensions substituted...');
% CANNOT write into d2xi_solved_un: the 1x6 sym is not writen in as function
matlabFunction(d2xi_solved_un,'file','two_d_spine_accel','Vars',{xi,tensions_un});

%PROGRESS_BAR
disp('     Writing functions for dynamics approach #2:');
disp('     (Note: please use the same accel function as approach 1).');
disp('     Writing the tensions function with barrier included...');
matlabFunction(tensions_sub_barrier,'file','two_d_spine_tensions_barrier','Vars',{xi,u});

%PROGRESS_BAR
disp('     Writing functions for dynamics approach #3:');
disp('     Writing xi_dot function...');
% matlabFunction(xi_dot_soln,'file','two_d_spine_xi_dot','Vars',{xi,u});

%PROGRESS_BAR
disp('     Writing functions for dynamics approach #4:');
disp('     Writing xi_dot, with barrier function...');
% matlabFunction(xi_dot_soln_barrier,'file','two_d_spine_xi_dot_barrier','Vars',{xi,u});

% In order to keep these symbolic variables for later,
% let's write them to a .mat file.%PROGRESS_BAR
disp('     Saving a .mat file with all the symbolic variables...');
% rename the variables so they make more sense.
two_d_spine_tensions = tensions_sub;
two_d_spine_accel = d2xi_solved_un;
two_d_spine_tensions_barrier = tensions_sub_barrier;
two_d_spine_xi_dot = xi_dot_soln;
two_d_spine_xi_dot_barrier = xi_dot_soln_barrier;
save('two_d_dynamics_symbolic_vars.mat', 'two_d_spine_tensions', 'two_d_spine_accel', ...
    'two_d_spine_tensions_barrier', 'two_d_spine_xi_dot', 'two_d_spine_xi_dot_barrier');

%% Script has finished.

%DEBUGGING
disp('Complete!');

