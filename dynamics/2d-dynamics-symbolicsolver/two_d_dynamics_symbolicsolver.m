% two_d_dynamics_symbolicsolver.m

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

%% 1) Prepare the workspace:
clc;
clear all;
close all;

% A flag to turn debugging on and off.
% This is useful to see the symbolic variables that are created.
debugging = 1;

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
%N = 3;
N=2;

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
  
% again, note the ' transpose above.
% Given the a matrix above, we can pick out how many
% point masses there are per unit:
num_pm_unit = size(a,2);
% ... again noting that a is transposed.
  
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
    unit_total_mass = 0.142;
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

% Lastly, define the cables that connect adjacent units.
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
c_vert = 100;
c_saddle = 100;

connections = cell(num_pm_unit, num_pm_unit);
% NOTE that these are assuming that a 'lower' unit is the 'from', 
% and it is connecting to a 'to' that's 'above' it.
% For example, with the inverted-Y spine vertebra, node 4 from a lower
% vertebra connects to node 3 of the vertebra above it, 
% but it's not true that node 3 from a lower vertebra connects to 
% node 4 of the vertebra above it.
% For our specific example, the following nodes are connected,
% with 2 vertical and 2 saddle cables:
connections{2,2} = [k_vert, c_vert];
connections{3,3} = [k_vert, c_vert];
connections{4,2} = [k_saddle, c_saddle];
connections{4,3} = [k_saddle, c_saddle];

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

% NOTE that we constrain these symbolic variables
% to be in the real numbers. Complex-number distance vectors
% don't make sense here.
% The state vector is then
xi = sym('xi', [num_states, 1], 'real');

% We'll be solving for the derviative of the state vector,
% as in the equation xi_dot = f(xi, u), so we need to declare
% that as a symbolic variable also.
% This will be our final result at the end of the script.
xi_dot = sym('xi_dot', [num_states, 1], 'real');

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
% Call it "t" for tension
tensions_unsolved = sym('t', size(tensions));

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
% [ rx_from, rzfrom; rx_to, rz_to; 2, 3]'
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
disp('Substituting system states back into r_dot...');
% Swap out all the 'dxi(whatever)' for the xi+3 coordinate
r_dot = replace_derivatives(r_dot, xi, num_states_per_unit, debugging);

%% 7) Calculate the kinetic and potential energy, and the Lagrangians, for the whole system.

%PROGRESS_BAR
disp('Solving for L = T-V for each unit...');

% TO-DO: check all these calculations. The results
% seem a bit too simple...

% Use the same indexing as in the previous section.
% HOWEVER, we need to calculate the Lagrangian for the first
% unit, too, I think.
% TO-DO: should I be setting Lagrangian(0)=0 since it's not moving,
% or should it be nonzero since the masses have potential energy?
% ANSWER, MAYBE: it doesn't matter, since that Lagrangian is a constant scalar,
% and we never use the Lagrangians directly: they're always as derivatives,
% so no matter what Lagrangian(0) is, diff(Lagrangian(0)) is always 0.
for k=1:N
    % For this unit, calculate the kinetic energy, (1/2)mv^2,
    % as in T = (1/2)* sum( m(i)*r_dot(i)^2), for point masses i.
    % Add in an extra simplify step to make things go faster, later.
    T = sym;
    for p=1:num_pm_unit
        % Add the kinetic energy for this point mass,
        % noting that the vector m has the masses for each point mass.
        % Here, we use r_dot'*r_dot as the 2-norm.
        T = T + (1/2) * m(p)*(r_dot(:,p,k)'*r_dot(:,p,k));
    end
    % Similarly, calculate the potential energy for this unit.
    % The potential energy is gravity times the z component of each unit.
    % Since the point masses may have different mass, do a loop.
    V = sym;
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
        % The first term is the full time derivative of (partial L / partial xi_dot).
        % Note that I'm not storing L_xi_dot as a symbolic variable itself,
        % it gets overwritten here. That's just because it is never used directly,
        % we only calculate it in order to calcualate ddt_L_xi_dot.
        % Note that we must index into 'Lagrangrian' via k+1 and not k,
        % since Lagrangian(1) is the constant value for the not-moving unit.
        L_xi_dot = diff(Lagrangian(k+1), xi(p));
        % Then take the full time derivative:
        ddt_L_xi_dot(count) = fulldiff( Lagrangian(k+1), sym2cell(xi));
        % Replace the dxi* terms in the symbolic variable:
        ddt_L_xi_dot(count) = replace_derivatives(ddt_L_xi_dot(count), xi, num_states_per_unit, debugging);
        
        % The second term is derivatives with respect to the POSITION variables,
        % so subtract away the +3 offset. 
        % My apologies for all this index-shuffling.
        L_xi(count) = diff(Lagrangian(k+1), xi( p - velocity_start_offset));
        
        % Do a quick simplify. Prior work used a parallel pool here,
        % might as well do that again.
        
        
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
% The dimensions of fx(i) are 

% CHECK: do my ddt_L_xi_dot have 2nd-order derivatives in them?
% Jeff's do. That's where he gets the xi_dot terms to solve for.

%% 9) Now, start calculating the cable forces.
% Calculate the lengths of each cable, as well as the rate of change
% of those lengths, so that we can calculate F = kx - c dx/dt.

%PROGRESS_BAR
disp('Calculating cable tensions...')

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
    for k=1:num_cables_per_unit
        for p=1:num_cables_per_unit
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
                dlengths_dt(cable_num) = simplify( ...
                    replace_derivatives(dlengths_dt(cable_num), xi, num_states_per_unit, debugging), num_simplify_steps);
                
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
                tensions(cable_num) = ...
                    connections{k,p}(1) * (lengths(cable_num) - u(cable_num)) ...
                    - connections{k,p}(2) * dlengths_dt(cable_num);
                
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
            F_cable = tensions_unsolved(k) * (tension_points(:,2,k) - tension_points(:,1,k));
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
                    disp(strcat('Contribution to global force for unit ', num2str(i), ' in direction ', num2str(q), ' was ', char( force_addition )));
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
            F_cable = tensions_unsolved(k) * (tension_points(:,1,k) - tension_points(:,2,k));
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
                    disp(strcat('Contribution to global force for unit ', num2str(i), ' in direction ', num2str(q), ' was ', char( force_addition )));
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
pf1 = parfeval(pools, @simplify, 1, global_forces(1,1), 'Steps', num_simplify_steps);
pf2 = parfeval(pools, @simplify, 1, global_forces(2,1), 'Steps', num_simplify_steps);
pf3 = parfeval(pools, @simplify, 1, global_forces(3,1), 'Steps', num_simplify_steps);
%pf4 = parfeval(pools, @simplify, 1, global_forces(1,2), 'Steps', 10);
%pf5 = parfeval(pools, @simplify, 1, global_forces(2,2), 'Steps', 10);
%pf6 = parfeval(pools, @simplify, 1, global_forces(3,2), 'Steps', 10);

disp('Fetching simplified outputs from parallel pool...');
disp('Fx, Unit 1...');
global_forces(1,1) = fetchOutputs(pf1);
disp('Fz, Unit 1...');
global_forces(2,1) = fetchOutputs(pf2);
disp('Ftheta, Unit 1...');
global_forces(3,1) = fetchOutputs(pf3);
%disp('Fx, Unit 2...');
%global_forces(1,2) = fetchOutputs(pf4);
%disp('Fz, Unit 2...');
%global_forces(2,2) = fetchOutputs(pf5);
%disp('Ftheta, Unit 2...');
%global_forces(3,2) = fetchOutputs(pf6);

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

%% 13) Equate the LHS and RHS of Lagrange's equations and solve!

disp('SOLVING LAGRANGES EQUATIONS...');

% The independent variables here are dxi4, dxi5, dxi6
[dxi4, dxi5, dxi6] = solve( ddt_L_xi_dot(1) - L_xi(1) == global_forces(1), ...
                   ddt_L_xi_dot(2) - L_xi(2) == global_forces(2), ...
                   ddt_L_xi_dot(3) - L_xi(3) == global_forces(3), ...
                   xi_dot(4), xi_dot(5), xi_dot(6));
               
% This doesn't work as of 2016-11-19, 7pm.

%% Script has finished.

%DEBUGGING
disp('Complete!');


% SOME BACKUP CODE IF NEEDED:

%     % As before, calculate the indices into the state vector xi.
%     % This is needed to specify the independent variables for differentiation.
%     % The state variables for this unit start at intervals of num_states_per_unit apart,
%     % and end at the next interval of num_states_per unit.
%     % For example, in the 6-state-per-unit spine vertebra, these
%     % intervals are 1-6, 7-12, 13-18, ...
%     unit_index_start = 1 + (k-1)*num_states_per_unit;
%     unit_index_end = (k)*num_states_per_unit;
%         indep_vars = sym2cell(xi(unit_index_start:unit_index_end));
%         %r_dot(:,p,k+1) = fulldiff(r(:,p,k+1), indep_vars);

%     % Pick out the indices of the two connection points
%     % for this specific cable.
%     % The 'from' and 'to' unit can be found the following way:
%     
%     % From: round down the index i divided by the number of cables per unit,
%     % noting that an offset of 1 is required in two places. 
%     % For our Y-spine example, cables 1-4 should have a 'from' of 1,
%     % and cables 5-8 should have a 'from' of 2.
%     from_unit = floor((i-1)/num_cables_per_unit) + 1;
%     % To: we could either round up the division of the index i by that
%     % same thing, or, we could just add one to from_unit, noting that
%     % we're always assuming that cables connect two adjacent units.
%     to_unit = from_unit + 1;

% tension_vec = sym('tension_vec', [num_cables, 2]);
% % Since we're assuming positive tension,
% % that means that each of these stored cable forces
% % will be acting in the + direction against the "upper" of 
% % the two units it connects to.
% % TO-DO: work through what happens when the tension forces
% % are in either the -x or -z direction, as with the saddle cables.

%                 % Now that the tension is expressed as a scalar, project it 
%                 % along the position vector of the cable.
%                 % TO-DO: could we calculate forces as vectors directly and skip this step?
%                 tension_vec(cable_num) = tensions(cable_num) * (to_node - from_node);
%                 % ...note that these vectors are now in the LOCAL coordinate system 
%                 % of the "to" unit.
%                 % Do a simplification step, as always.
%                 tension_vec(cable_num) = simplify(tension_vec(cable_num));




