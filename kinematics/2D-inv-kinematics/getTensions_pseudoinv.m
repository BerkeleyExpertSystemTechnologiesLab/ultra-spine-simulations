%% Find Cable Forces and Rest Lengths in 2D Spine
% This function minimizes the force density values of the cables in a 2D
% spine consisting of two stacked tetrahedra. The bottom tetra is fixed,
% and the location of the top tetra is defined by the input xi. We assume
% that the spine is sitting on a surface, so that there are vertical
% reaction forces at each of the two bottom nodes in contact. The function
% returns the cable tensions and the corresponding rest lengths.
%
% Authors: Drew Sabelhaus, Mallory Daly, and Ellande Tang
% Created: 12/8/16
% Modified: 2/7/18

% THIS FUNCTION USES THE PSEUDOINVERSE FORMULATION of the inverse
% kinematics. The optimization is now inequality constrained. We'll see if
% that works better.

function [tensions, restLengths, A, p, qOpt_lax] = getTensions_pseudoinv(xi, spineParameters, minCableTension)
% Outputs:
%   There are four solvers at work here, and this function is designed to
%   return the A, p, and qOpt for all four. They are:
%
%   1) Mal/Ellande, Equality Constraint (original version)
%   2) Skelton, Equality Constraint (uses algorithmic formulation of Aq=p
%       and is thus more theoretically grounded)
%   3) Mal/Ellande, Inequality Constraint. Easier for quadprog to solve.
%   4) Skeltn, Inequality Constraint. This is the Friesen ICRA 2014 paper,
%       and is what we want to use from now on (if it works...)

%% Spine Parameters

% all parameters hard-coded.

% Number of bars, cables, and nodes
r = 6; % bars
s = 4; % cables
n = 8; % nodes
% Parameter that we need for the reaction forces: 

% Geometric parameters
ll = spineParameters.l; % m, length of long bars
h = spineParameters.h; % m, height from top to bottom of tetra
ls = h/2; % m, length of short bar
w = sqrt(ll^2-(h/2)^2); % m, width from center of tetra

% Mass and force parameters
g = spineParameters.g; % m/s^2, acceleration due to gravity
M = spineParameters.total_m; % kg/tetra
springConstant = spineParameters.k_vert; % structure has vertical and horizontal k, but they're the same, so ignore for now
m_node = spineParameters.m; % kg/node FOR ALL NODES!
% Thus, change it to one node:
m_node = m_node(1); % To-do: check this!!!
M = spineParameters.total_m; % kg/tetra
springConstant = spineParameters.k_vert;

%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."

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
 
 % Checked on 2018-02-08 and confirmed that this matches Jeff Friesen's
 % paper's formulation of the connectivity matrix. In (s+r) x n, 10x8.

% For later, let's store a reduced C matrix that only contains the cables.
% s=4 cables
Cs = C(1:s,:);

%% Nodal Positions
% Coordinate system such that nodes 3 and 4 are on the x axis and nodes 1
% and 2 are centered on the z axis.

% Nodal positions of bottom tetrahedra
%           1    2    3    4
x_bot = [   0   -w    w    0]';
z_bot = [   0 -h/2 -h/2  h/2]';
% Nodes are:
%   1) center node
%   2) bottom left
%   3) bottom right
%   4) top

% Nodal positions and rotation of top tetra
% NOT THE NODES! This is a different dimension that x, z_bot
x_rigid_top = xi(1);
z_rigid_top = xi(2);
theta = xi(3);

% So... let's calculate the nodal positions of the top vertebra.

% Rotation matrix for given angle of theta
% For our coordinate system, with Z upwards, X to the right, then a rotation 
% about +Y (into the page) of +theta turns the vertebra clockwise. This is
% because we're looking at the vertebra "from behind."


% THIS ROTATION MATRIX was wrong, with the (-) in the (2,1) element.
% Now moved to the (1,2) element.
rot = [ cos(theta), -sin(theta);
       sin(theta), cos(theta)];
   
% Need to multiply each node by the rotation matrix and add the (x,z)
% offset from the xi state vector. Copy out (x,z) offset to a 2x4 matrix to
% make this easier.
xz_offset = repmat([x_rigid_top; z_rigid_top], 1, 4);

% Translate and rotate position of fixed nodes to get coordinates of free
% nodes
xz_top = rot*[x_bot'; z_bot'] + xz_offset;

% Extract x and z of free nodes
% NOTE: these vectors are now the same size as x, z_bot.
x_top = xz_top(1,:)';
z_top = xz_top(2,:)';

% Combined nodal positions
% As per Jeff's formulation, each of these should be vectors in R^n, where
% n is the number of nodes. So, we're looking at each as an 8-vector.
x = [x_bot; x_top];
z = [z_bot; z_top];

% DEBUGGING
% Plot nodal positions
% ...looks good. 
% figure
% plot(x_bot,z_bot,'k.','MarkerSize',10)
% hold on
% plot(x_top,z_top,'r.','MarkerSize',10)
% return;

%% Lengths of Bars and Cables

% Rows 1-4 are cables
% Rows 5-10 are bars
l = [norm([x(2),z(2)]-[x(6),z(6)]); %  1
     norm([x(3),z(3)]-[x(7),z(7)]); %  2
     norm([x(4),z(4)]-[x(6),z(6)]); %  3
     norm([x(4),z(4)]-[x(7),z(7)]); %  4
     ll;                            %  5
     ll;                            %  6
     ls;                            %  7
     ll;                            %  8
     ll;                            %  9
     ls];                           % 10

% Cable diagonal length matrix
l_cables = l(1:s);
L_cables = diag(l_cables);

%% Solve for Reaction Forces
% Assume spine is sitting on a surface. Then there are vertical reaction
% forces at nodes 2 and 3.
  
% Mallory and Ellande version: solve by hand.
% Solve AR*[R2; R3] = bR, where AR will always be invertible
AR = [1 1; 0 (x(3)-x(2))];
bR = [2*M*g; M*g*(x(1)-x(2))+M*g*(x(5)-x(2))];
R = AR\bR;
R2 = R(1);
R3 = R(2);

% Better version (to publish later): solve algorithmically,
% specify the reaction forces in a certain way and then solve for them.
% Lots of difficult indexing here so we'll save for later...

% First: how about we specify which nodes are fixed. We assume that these
% experience external reaction forces that must be balanced out.

% Here, nodes 2 and 3 contact the ground. Assume built-in.
%        
%fixed = [0; % 1
%         1; % 2
%         1; % 3
%         0; % 4
%         0; % 5
%         0; % 6
%         0; % 7
%         0]; % 8

% ...leave this for future research. I'm gonna have to do it some day, like
% for the Laika model, but for now I'm just going to let Mallory and
% Ellande's calculations stand. (I did the math by hand and I agree with
% them, although it would be interesting to see if adding X-direction reactions
% will change anything. Maybe makes indeterminate? Built in vs. rollers
% might not matter, either, since we're keeping the nodes in the same
% position.

% Maybe future research into tensegrity systems, is to
% include reaction forces into the optimization. See something like this
% example problem,
% http://www.unm.edu/~bgreen/ME360/Statics%20-%20Truss%20Problem.pdf, which
% might imply that we augment the cable densities 'q' with some reaction
% forces 'R'.

% To envision why we need these external reaction forces: imagine a
% tensegrity triangle, floating in space, with the same gravitational force
% applied to all nodes. Gravity will cancel out, overall. Then, imagine if
% one edge of the triangle hits thr ground. The top node will now be
% "lower" than if it was free floating, if tensions are equivalent, so in
% other words, the tensions must be DIFFERENT in order to keep the
% structure in the same position.



%% Equilibrium Force Equations - Both Mal\Ellande and Skelton


% First: the Mallory/Ellande approach.
% Solve sum forces, moments = 0, get a system of 4 equations (x, z, for
% vertebra 1, 2.)

% Create vector of distance differences
dx = C*x;
dz = C*z;

% Some extra math to try out what we write in linear algebra for the paper,
% "L" paper is "Lc" here
% Lc = [ C*x; C*z]
% Af = 

% Define A*q = p, where q is a vector of the cable force densities and p is
% a vector of the external forces. Note that A is not a full rank matrix
% (not invertible).
A = [ -dx(1) -dx(2) -dx(3) -dx(4);  % horizontal forces, bottom tetra
       dx(1)  dx(2)  dx(3)  dx(4);  % horizontal forces, top tetra
      -dz(1) -dz(2) -dz(3) -dz(4);  % vertical forces, bottom tetra
       dz(1)  dz(2)  dz(3)  dz(4)]; % vertical forces, top tetra

% In the manual derivation formulation, each row of A is the force balance
% in one dimension, for all cable tensions. Example: A(1,:) is the
% x-coordinate balance for all 4 of the cables, size( A(1,:), 2) is 4.

% Moments about each tetra.
% Takes node numbers, which are indexed into the x and z matrices, and calculates the force by using those positions.
% Inputs are:
%   a = node number, for center of mass node
%   b = node number, for the end of the moment arm (where the cable is applying a force)
%   c = node number, the "other end" of the cable. This is required to calculate a force in Newtons, since 
%           the q vector is force DENSITY, in N/m. 
qfun = @(a,b,c) (x(b)-x(a))*(z(c)-z(b)) - (z(b)-z(a))*(x(c)-x(b));
% The bottom row appended to A represents the contribution of each of the 4 cables, contributing to the moment around 
% the Y-axis of the top vertebra.
% Here, note that the "center of mass" node is 5 for the top vertebra, 1 for the center of mass of the bottom node.
A = [A;
     qfun(1,2,6) qfun(1,3,7) qfun(1,4,6) qfun(1,4,7);
     qfun(5,6,2) qfun(5,7,3) qfun(5,6,4) qfun(5,7,4)];

%disp('A, augmented, from manual derivation:');
%A


% disp('A, initially, from manual derivation:');
% A
% size(A)
% disp('Rank of A with rigid body derivation:');
% rank(A)


% Alternatively: Skelton/Friesen formulation, using C ' diag(C x) etc.

% A= [C_A' *diag(C_A*tetraNodes(:,1));
%        C_A' *diag(C_A*tetraNodes(:,2));
%        C_A' *diag(C_A*tetraNodes(:,3))];

% Jeff's tetraNodes seem to be [x, z]. C_A is "full C, including cables and
% bars."

% For reference:
% tetraNodesPreTransform = [L/2   0     -h/2  1; %A
%                          -L/2   0     -h/2  1; %B                 
%                          0    -L/2    h/2  1; %C
%                          0     L/2    h/2  1];%D

% 

% disp('Debugging: lets see the Skelton full A matrix, with size:');
% 
A_skelton = [ C' * diag(C * x);
              C' * diag(C * z)];
% size(A_skelton)
% 
% disp('Rank of A with graph formulation:');
% rank(A_skelton)

% ...this is size s+r = 10, so the size of q will be 10.
          
% Seems different. Different dimensions, at least!
% Let's see what it looks like when we only pick out the cables, and ignore
% the bars. (After all, we don't care about the compressive forces in the
% bars right now.)
% We've already cut out the C to include cables only (Cs),
%A_skelton_c = [ Cs' * diag(Cs * x);
%                Cs' * diag(Cs * z)]
%size(A_skelton_c)

% IMPORTANT NOTE: On 2018-03-30, talked with Jeff. Apparently the 
% reduction in rows of A must happen AFTER the pseudoinverse, below.
% So, DO NOT use the A_skelton_c matrix.

% DEBUGGING:
% Maybe we can't just cut the configuration matrix off just here. Jeff does
% it after calculating A, and going by the columns. Let's do it that way
% and compare. s=4 here. 
%A_skelton_c_after = A_skelton(:,1:s)
% Result: this is the same as A_skelton_c.

%% Skelton's external force vector

% Let's do the forces applied to each node, so we can use the skelton
% configuration. (Drew thinks this gets rid of the need for moment
% balance???)
% We'll store this as the p_skelton vector. Copied from my notes below:

% In the manual derivation (Mal/Ellande) formulation, each row of A (below) is the force balance
% in one dimension, for all cable tensions. Example: A(1,:) is the
% x-coordinate balance for all 4 of the cables, size( A(1,:), 2) is 4.

% In the Skelton_c formulation, one row of A is one cable's contribution to
% the force balance, for one node. Example: A(1,:) is the forces at node 0,
% where no cables touch, thus has no contribution to the force balance.

% The p that's used here is p = [p_x; p_z], where p_x is the vector of
% external loads applied to each node (in order) in the x-direction.
% So, to compare it to Mallory/Ellande's p, we'd want to... take this p and
% make sure all the x-forces balance out for one rigid body. For example,
% that would be forces at the coordinates for x(1:4) for the first vertebra, which would then correspond
% to p(1:4). And, p(5:8) == x(5:8), p(9:12) = z(1:4), etc. (but here, the x
% and z are distances, and I'm just talking indices.) 
% Maybe the comparison is p_mallory(1) == p_skelton(1:4), which sums
% x-forces for the bottom rigid body, for example.

% This leads to a p_skelton of size 2 * n, since p_x is of size n and p_z
% is also of size n.
% Importantly, this is INDEPENDENT of the number of bars AS WELL AS the
% number of cables!

% We'll need to go gravity here for everything, but then also include the
% vertical reaction forces for nodes 2 and 3 (which sit on the ground.)

p_skelton = zeros(16,1);
% v1, x, 4 nodes:
p_skelton(1:4) = zeros(4,1);
% v2, x, 4 nodes:
p_skelton(5:8) = zeros(4,1);
% v1, z, 4 nodes:
p_skelton(9:12) = - m_node * g; % not the total mass of one vert!
% v2, z, 4 nodes:
p_skelton(13:16) = -m_node * g;
% BUT we should adjust by adding the reaction forces to node 2 and 3 of the
% lower vertebra. Otherwise, the problem becomes infeasible, the tensegrity
% always falls / cannot have zero force.
% Nodes 2 and 3 in Z are n + 2, n+ 3, since 8 nodes means 8 coordinates in
% the X direction first. That's 10 and 11. (not 14, 15 since those would be
% the top vertebra instead.

% On 2018-03-30: it seems that Mallory and Ellande have defined R2 and R3
% to be downward?? They subtract R2 and R3, not add.
p_skelton(10) = p_skelton(10) + R2;
p_skelton(11) = p_skelton(11) + R3;

% DEBUGGING: Let's see if the Skelton formulation can solve for a robot
% floating in space (no gravity, no external forces.)
%p_skelton = zeros(16,1);

% disp('p_skelton, size 2*n x 1, is:')
% p_skelton
% 
% disp('Comparison: p (from Mallory/Ellande) is:');
% p

% YES! This works. p_skelton(1:4) == p_mallory(1), to numerical precision. Except, the signs seem
% to be flipped for the x and z for the 2nd vertebra, elements p_mallory(3)
% and 4 versus p_skelton(9:12), 13:16.
% Might need to flip gravity...?

% disp('Do the external forces sum to zero? Sum is:')
% sum(p_skelton)

% %% A quick check on the number of solutions for both sets of equations.
% 
% disp('Checking the number of solutions for both formulations of Aq=p:');
% 
% % Solutions exist to nonsquare Aq=p if A pin(A) p = p.
% 
% disp('A pinv(A) p , for Mal/Ellande:');
% 
% A * pinv(A)
% A * pinv(A) * p
% 
% A * pinv(A) * p - p

%% Mallory / Ellande's External Force Vector

% Note R2 and R3 cannot be negative, but this constraint is not imposed
% here. However, a condition under which R2 or R3 becomes negative creates
% an infeasible problem for the cable tensions, so the issue solves itself.
% It's possible that a more generalized problem would need to solve for
% reaction forces using a solver in order to impose this constraint.

% p = [ 0; 0; M*g-R2-R3; M*g; 0; (R2-R3)*w];

% The p vector is the right-hand side of the force balance.
% The zeros come from \sum F = 0 in the X direction, for vertebra 1 and vertebra 2.
% 3rd, 4th rows are \sum F = 0 in the Z direction, for vertebra 1 and vertebra 2.
% 5th row is moment balance, = 0.
%disp('p, manual derivation:');

% ON 2018-03-18: ...shouldn't this be -mg, and then +R2, +R3? Change it and
% see what happens:
%p = [ 0; 0; M*g-R2-R3; M*g; 0];
% ...seems to be a much more reasonable answer. Rrrrgh need to re-run the
% MPC now at least.
% On 2018-04-05: did sum of moments around both vertebrae, now p \in R
% (3*d), or: 2 x-balance, 2 z-balance, 2-moment balance.
% Also, the moment balance for the bottom vertebra has an external
% component, which is the reaction forces at nodes 2 and 3, so these are -(x(2)*R2) + (x(3)*R3).
% Or, in the terms we'll use in the paper, (x(1) - x(2))*R2 + (x(3) - x(1))*R3

% Written by hand:
%p = [ 0; 0; -M*g + R2 + R3; -M*g; -(x(2)*R2) + x(3)*R3; 0];

% BETTER VERSION USING LINEAR ALGEBRA:
% We should be able to collapse the p_skelton into this p, by using a
% kronecker product to collapse each set of 4 rows into one.
% That's I_4 (kron) 1_s^\top, e.g., four rows of blocks of row vectors of ones.

p_collapsed = kron( eye(4), ones(1,4) ) * p_skelton;

% Then add back the moment balance constraints.
p = [ p_collapsed; -(x(2)*R2) + x(3)*R3; 0];

% ...just to confirm, we did the direction properly, right???
%p = [ 0; 0; M*g - R2 - R3; M*g; -(x(2)*R2) + x(3)*R3; 0];
% ...result: this seemed not to make a difference. That sounds right, the
% vertebra is moving "way too fast" anyway.

%% Solve Problem for Minimized Cable Tension - Mal/Ellande, Equality Constraint

% disp('Solving Mal/Ellande, Equality Constraint:');
% 
% % Solve with QUADPROG
% % The objective function, described by 0.5 * x' * H * x + f' * x
% % This minimizes the total force density in the cables, quadratically.
% H = 2*eye(4);
% f = zeros(4,1);
% % The force and moment balances are constraints for the optimization
% Aeq = A;
% beq = p;
% % Inequalities ensure a minimum cable tension
% Aineq = -L_cables;
% bineq = -minCableTension*ones(s,1);
% % opts = optimoptions(@quadprog,'Display','notify-detailed');
% [qOpt, ~, exitFlag] = quadprog(H,f,Aineq,bineq,Aeq,beq);
% 
% disp('qOpt, Mal/Ellande, Equality:');
% qOpt
% 
% if exitFlag == 1
%     tensions = L_cables*qOpt; % N
%     restLengths = l_cables - tensions/springConstant;
%     if any(restLengths <= 0)
%         display('WARNING: One or more rest lengths are negative. Position is not feasible with current spring constant.')
%     end
% else
%     display(['Quadprog exit flag: ' num2str(exitFlag)])
%     tensions = Inf*ones(s,1);
%     restLengths = -Inf*ones(s,1);
% end

% DEBUGGING:
% Out of curiosity, let's see what we get when the moments are ignored:
%Aeq_ig = A(1:4, :);
%beq_ig = p(1:4);
%[qOpt_ig, ~, ~] = quadprog(H, f, Aineq, bineq, Aeq_ig, beq_ig);

%qOpt_ig
% seems to be the same for now? Is there a good reason why we can ignore
% moments in the structure? YES, because this is just a point mass system.
% That does not, however, mean we can ignore moments for the *whole
% structure* when calculating the reaction forces!! E.g., the way that Mal
% and Ellande did this, cutting off the last bit of "p" does not seem to
% make a difference, BUT we did need to do sum moments when calculating R2,
% R3 that went into p(3:4) as force balance.

%% Solve Problem for Minimized Cable Tension - Skelton/Friesen, Equality Constraint

% % BUT for the skelton formulation, we really do need to
% % include +R2 and +R3 in the p vector.
% 
% % Here, Drew does not think we can ignore the bars. Later, when we do the
% % pseudoinverse, maybe we could - although I still have to test Jeff's
% % theory on the matter, and maybe it's easier to just include all of them
% % and then put no minimization weight on the bars and no constraints (e.g.
% % set H = [eye; zeros] instead of eye, and c = [cs; 0] or something.)
% 
% disp('Using the Skelton formulation, equality:');
% 
% %DEBUGGING
% %A_skelton
% %p_skelton
% 
% % From above, we had H = 2*eye(4); f = zeros(4,1); since the objective
% % function is only q'q for now (NOT RELAXING YET).
% % Now, we're doing the full 10-vector for q, since need to include the bars
% % in the equality constraint.
% % BUT! Need to place zero weight on the bar forces. so,
% H_sk = zeros(s+r, s+r);
% H_sk(1:s, 1:s) = 2 * eye(s); % s+r = 10;
% f_sk = zeros(s+r, 1); % zero vector size 10
% 
% % The force and moment balances are constraints for the optimization
% % A q = p,
% Aeq_sk = A_skelton;  % size 16 x 10
% beq_sk = p_skelton;  % size 16 x 1
% % ...these work with q size 10 x 1
% 
% % Inequalities ensure a minimum cable tension. With
% %l_cables = l(1:s);
% %L_cables = diag(l_cables);
% % ...the L_cables matrix is diagonal with size s x s. We need to insert it into
% % a zeros matrix with s+r columns now, to align with q of size s+r x 1.
% Aineq_sk = zeros(s, s+r);
% Aineq_sk(1:s, 1:s) = -L_cables;
% % We can then keep the original b_inequality, since Aineq_sk * q is of size
% % s x 1. E.g., only need the cable force constraints, and padding with
% % zeros such that the dimensions fit for q.
% bineq_sk = -minCableTension*ones(s,1);
% 
% disp('Skelton Formulation Solution, Equality Constraint:');
% 
% % let's do the optimization using the skelton formulation. Should be:
% 
% % DEBUGGING: is there a solution without constraining the cables?
% % E.g. is there any solution to A_skelton * q = p?.
% [qOpt_sk, ~, exitFlag] = quadprog(H_sk, f_sk, [], [], Aeq_sk, beq_sk);
% 
% %[qOpt_sk, ~, exitFlag] = quadprog(H_sk, f_sk, Aineq_sk, bineq_sk, Aeq_sk, beq_sk);
% 
% qOpt_sk

%% Solve Problem for Minimized Cable Tension - Mal/Ellande, Inequality Constraint
% As of 2018-03-30, this is the version returned as 'restlengths' and is
% used in the 2D MPC.

disp('Relaxing the Mal/Ellande formulation:');

% Let's just try to relax the Mal/Ellande formulation and see if the answer
% changes.

%disp('Solving the inequality relaxed version:');

% As = A (it's already strings only, column-wise, and isn't 3n rows but it
% should still match the other matrices.
As = A;
% Need the pseudoinverse
Apinv = pinv(As);
% The quantity V is I - A+ A
% What's the side of A+ * A?
% Answer: seems to be 4. Is that num cables???
ApA = Apinv * As; % is, necessarily, square.
ApA_dim = size(ApA, 1); % ...so either dimension can be taken.
V = eye(ApA_dim) - ApA;

% Now we can formulate the relaxed problem:
H_lax = V' * V;
f_lax = V' * Apinv * p;
A_ineq_lax = -V;
% Calculate the minumum force densities from the min cable tension and
% length. l_cables is an 8 dimensional vector. (capital L is diag version.)
% let's just do element-wise division for ease.
min_force_densities = minCableTension*ones(s,1) ./ l_cables;
% the b term here is Apinv * p - c. 
% MIGHT BE NEGATIVE OF THIS. actually, prob not. Direction of inequality
% flipped between Jeff's paper and the use of quadprog.
b_ineq_lax = Apinv * p - min_force_densities;


% Finally, let's see if we can solve. Dropping the equality terms.
[wOpt_lax, ~, ~] = quadprog(H_lax, f_lax, A_ineq_lax, b_ineq_lax);

% ...actually, this is the w vector, NOT q. Still need to do q=A_g*F + V*w;
qOpt_lax = Apinv*p + V*wOpt_lax
% ...this is EXACTLY THE SAME AS THE EQUALITY CONSTRAINT SOLUTION!!!

% Change the rest lengths to the relaxed values and see what changes in the
% MPC simulation:
tensions = L_cables*qOpt_lax; % N
restLengths = l_cables - tensions/springConstant;

% ...this seems to work now. We should check and confirm that both
% solutions stabilize the vertebrae.


%% Solve Problem for Minimized Cable Tension - Skelton/Friesen, Inequality Constraint

% Jeff's code for DuCTT:
% % Solve with YALMIP
% yalmip('clear')
% q = sdpvar(s,1);
% obj = q'*q;
% constr = [A*q == p, L_cables*q >= minCableTension*ones(s,1)];
% options = sdpsettings('solver','quadprog','verbose',0);
% sol = optimize(constr,obj,options);

% Jeff's code is:
%     C_A=[C; B];
%     A= [C_A' *diag(C_A*tetraNodes(:,1));
%         C_A' *diag(C_A*tetraNodes(:,2));
%         C_A' *diag(C_A*tetraNodes(:,3))];
%     A_g = pinv( A);
%     A_g_A=A_g*A;
%     V=(eye(length(A_g_A))-A_g_A);
%     [Q,R,E] = qr(V);
%     [m , n] = size(R);
%     j=1;
%     i=1;
%     while i<=m
%         if norm(R(i,:))>10^-12
%             R_new(j,:)=R(i,:);
%             j=j+1;
%         else
%             i=m;
%         end
%         i=i+1;
%     end
%     V=Q(:,1:j-1);
%     
%     % Run the actual optimization for the inverse kinematics
%     $ quadprog is (H, f, A, b), 0.5*x'Hx + f'x, Ax <= b
%     w = quadprog( V(1:(N-1)*8,:)' * K_scale * V(1:(N-1)*8,:), ... % H, which is w'V'Vw 
%        V(1:(N-1)*8,:)'*K_scale*A_g(1:(N-1)*8,:)*F, ... % f, which is
%        w'*V'*A+p (Jeff has factor of 2 in paper, but since quadprog has
%        0.5 * H, works out.) This is not "+" but "A+", pinv.
%        -V(1:(N-1)*8,:), ... % A, -Vw
%        A_g(1:(N-1)*8,:)*F-pretension, ... % b, A+p - c , c=pretension.
%        [],[],[],[],[],options); % misc stuff.
%     q=A_g*F + V*w;

% disp('Relaxing the Skelton formulation:');
% 
% % 2018-03-30: ISSUE: the inequality constraint that constrains the cables
% % plus bars. We'd have to flip some signs to make sure that the densitites
% % for the bars were compressive! So, alternatively, let's do Jeff's code
% % and remove the force densities of the bars.
% 
% % Let's try relaxing the skelton formulation to see if we get a feasible
% % answer that way. Maybe it's just that quadprog has difficulty
% % initializing? (Or is it really that no solutions exist, and we've
% % formulated the problem incorrectly?)
% % See the comments below for what each term means.
% 
% % As of 2018-03-30, let's use the full A_skelton matrix, and then 
% % hack out the bars part AFTER the pseudoinverse.
% A_sk = A_skelton; % renamed.
% A_sk_pinv = pinv(A_sk);
% % We want to remove the last r columns of A_sk and the last r rows of
% % A_sk_pinv, so that A_sk is 16 x 4 and A_sk_pinv is 4 x 16
% A_sk = A_sk(:, 1:s);
% A_sk_pinv = A_sk_pinv(1:s, :);
% 
% % continuing with Jeff's derivation...
% ApA_sk = A_sk_pinv * A_sk;
% ApA_sk_dim = size(ApA_sk, 1); % ...so either dimension can be taken.
% V_sk = eye(ApA_sk_dim) - ApA_sk;
% 
% % Now we can formulate the relaxed problem:
% H_sk_lax = V_sk' * V_sk;
% f_sk_lax = V_sk' * A_sk_pinv * p_skelton;
% 
% % The inequality constraint is really V * w >= -A_sk_pinv * p + c
% % Or just flip the signs, since quadprog does A*x <= b, so this becomes
% % - V*w <= A_sk_pinv * p - c
% A_sk_ineq_lax = -V_sk;
% % Calculate the minumum force densities from the min cable tension and
% % length. l_cables is an s-dimensional vector. (capital L is diag version.)
% % let's just do element-wise division for ease.
% min_force_densities = minCableTension*ones(s,1) ./ l_cables;
% b_sk_ineq_lax = A_sk_pinv * p_skelton - min_force_densities;
% 
% % Finally, let's see if we can solve. Dropping the equality terms.
% [wOpt_sk_lax, ~, ~] = quadprog(H_sk_lax, f_sk_lax, A_sk_ineq_lax, b_sk_ineq_lax);
% 
% % ...actually, this is the w vector, NOT q. Still need to do q=A_g*F + V*w;
% qOpt_sk_lax = A_sk_pinv * p_skelton + V_sk * wOpt_sk_lax;
% % ...this is EXACTLY THE SAME AS THE EQUALITY CONSTRAINT SOLUTION!!!
% 
% disp('Optimal q, relaxed Skelton formulation:');
% qOpt_sk_lax

%% Check Distance Vectors Symbollically

% x1 = sym('x1','real');
% x2 = sym('x2','real');
% x3 = sym('x3','real');
% x4 = sym('x4','real');
% x5 = sym('x5','real');
% x6 = sym('x6','real');
% x7 = sym('x7','real');
% x8 = sym('x8','real');
% x = [x1 x2 x3 x4 x5 x6 x7 x8]';
% 
% z1 = sym('z1','real');
% z2 = sym('z2','real');
% z3 = sym('z3','real');
% z4 = sym('z4','real');
% z5 = sym('z5','real');
% z6 = sym('z6','real');
% z7 = sym('z7','real');
% z8 = sym('z8','real');
% z = [z1 z2 z3 z4 z5 z6 z7 z8]';
% 
% C*x
% C*z
% 
% q1 = sym('q1','real');
% q2 = sym('q2','real');
% q3 = sym('q3','real');
% q4 = sym('q4','real');
% qs = [q1 q2 q3 q4]'

end

% end of function.