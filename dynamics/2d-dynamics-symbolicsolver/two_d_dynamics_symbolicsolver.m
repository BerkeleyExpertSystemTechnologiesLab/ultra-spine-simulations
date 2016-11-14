% two_d_dynamics_symbolicsolver
% Berkeley Emergent Space Tensegrities Lab and Andrew P. Sabelhaus
% Copyright 2016
% This script symbolically solves Lagrange's equation for 2-dimensional
%   tensegrity systems which involve repeated "units," each of which are
%   a rigid body that is defined by point masses at nodes, and where all the
%   point masses are rigidly connected. Then, connection points between different
%   "units" are defined, as pairs of nodes that have a cable between them.
%   (TO-DO: what about for the "first" and "final" units?)
%   This is inspired by, and some of the terminology comes from, 
%   the type-II tensegrity spine under development in the BEST Lab.
% Many thanks to Jeff Friesen and Abishek Akella for earlier versions of this script.

% This script generates multiple m-files:
%   duct_accel
%   dlengths_dt
%   lengths

% Prepare the workspace:
clc;
clear all;
close all;

% A flag to turn debugging on and off.
% This is useful to see the symbolic variables that are created.
debugging = 1;

% Throughout this script, I output some messages that show the progress
% of the script during its calculations.
% Those are labelled 'PROGRESS_BAR'.

% Some other parameters to look at:
% the 'equal_masses' flag in the following section,
% ...

%% Geometric and physical parameters of the 2D tensegrity system

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
% Then, for N-1 of these structures, symbolic calculations
% of the velocities ...

% This script assumes that the first unit is fixed at the origin,
% and that the first node of this first unit is at exactly (0,0,0).
% There are N units, including the fixed one:
N = 3;

% Here, we define the positions of each of the point masses
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


%% Create the symbolic variables for the solver to use

%PROGRESS_BAR
disp('Creating symbolic variables...');

% The size of the system state varies with number of units.
% For each unit, there are six state variables: [x; z; theta];
% and their derivatives [dx; dz; dtheta];,
% which are the location of the origin of the coordinate system,
% and the rotation about that center of the coordinate system,
% for each unit.
% (TO-DO): what happens when the origin is NOT the center of mass???????????????????????????????????
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

%% TO-DO: do these dynamics need to change now that our spine
% is not symmetric in 3D?

% Answer: I think not. If we were calculating the rigid body dynamics,
% not using point masses, then yes we'd have to calculate the center of mass
% and express the rotational velocity/acceleration around that COM. But here,
% we're just using point masses. The rotation here is one coordinate frame
% versus another, which can be completely arbitrary, as along as it's consistent.

%% Next, constrain the first unit. 

%PROGRESS_BAR
disp('Adding constraints on the first unit...');

% The positions of its nodes are equal to a in both the
% local and global coordinate frames, since this unit does not move.
r(:,:,1) = a;

% The velocities of each point mass here are exactly zero.
% Make a matrix of zeros that fits exactly.
r_dot(:,:,1) = zeros(size(r_dot(:,:,1)));

%% Then, express the coordinates of each point mass in terms of the system states.

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
    disp(strcat('Assigning point mass locations for unit number: ', num2str(k+1)));
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
        % TO-DO: what's the direction of multiplication for rotation matrices????
        r(:,p,k+1) = R(:,:,k)*a(:,p) + xi(x_index:z_index);
    end
end

%% Similarly, express the velocities of each point mass in terms of xi.

%PROGRESS_BAR
disp('Assigning point mass velocities in terms of system states...');

% For each of the moving units...
for k=1:N-1
    %PROGRESS_BAR
    disp(strcat('Calculating symbolic derivatives of point masses for unit number: ', num2str(k+1)));
    % As before, calculate the indices into the state vector xi.
    % This is needed to specify the independent variables for differentiation.
    % The state variables for this unit start at intervals of num_states_per_unit apart,
    % and end at the next interval of num_states_per unit.
    % For example, in the 6-state-per-unit spine vertebra, these
    % intervals are 1-6, 7-12, 13-18, ...
    unit_index_start = 1 + (k-1)*num_states_per_unit;
    unit_index_end = (k)*num_states_per_unit;
    % For each of the point masses in this unit:
    for p=1:num_pm_unit
        % The velocity is the full derivative of position (with respect to time.)
        % Credit goes to Tim Jorris for the fulldiff function.
        % Note, however, that the independent variables must be passed in as a cell array.
        indep_vars = sym2cell(xi(unit_index_start:unit_index_end));
        r_dot(:,p,k+1) = fulldiff(r(:,p,k+1), indep_vars);
    end
end

% Now, to make the simplification easier on the solver later, substitute
% all the derivatives that fulldiff calculated back into the symbolic variables xi.
% Fulldiff writes these variables as, for example, dxi1, dxi1, etc., which are really
% xi4, xi5, etc.

%PROGRESS_BAR
disp('Substituting system states back into r_dot...');

% Though we know that there are only 3 variables here, let's still use 
% the num_states_per_unit variable to calculate which states are positions and
% which are velocities.
% This variable should be 3:
velocity_start_offset = num_states_per_unit/2;
% Iterate through all the units:
for k=1:N-1
    % Calculate the start and end indices for the positions at this index
    % This will be, for example, 1, 7, 13, ...
    unit_index_start = 1 + (k-1)*num_states_per_unit;
    % This will be, for example, 3, 10, 15, ...
    unit_position_end = unit_index_start + velocity_start_offset - 1;
    % For each of the position variables in this unit:
    for p=unit_index_start:unit_position_end
        % The old variable for this specific position variable starts with a 'd', as output
        % by fulldiff:
        oldvalue = strcat('d', char(xi(p)));
        % The new variable is the corresponding velocity state
        newvalue = char(xi(p+velocity_start_offset));
        
        %DEBUGGING
        if debugging
            disp(strcat('     p is: ', num2str(p)));
            disp(strcat('     oldvalue is: ', oldvalue));
            disp(strcat('     newvalue is: ', newvalue));
        end
        
        % Add the field/value pair to the struct
        %r_dot_sub_vars.eval(field) = value;
        
        % Perform the substitution for this field/value pair
        r_dot = subs(r_dot, oldvalue, newvalue);
    end
end

%% Calculate the kinetic and potential energy, and the Lagrangian, for the whole system.

%% Script has finished.

%DEBUGGING
disp('Complete!');


