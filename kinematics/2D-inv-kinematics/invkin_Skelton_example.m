%% Invkin Skelton Example
% Authors: Drew Sabelhaus
% Modified: 3/31/18

% This script tests out the Skelton formulation of Aq=p, for one vertebra
% without cables, just to get a better understanding of how quadprog works
% with the A matrix.

% function used to return [tensions, restLengths, A, p, A_skelton, qOpt, qOpt_lax]

%% prep the workspace
clear all;
close all;
clc;

%% some constants
dynamics_path = '../../dynamics/2d-dynamics-symbolicsolver';
addpath(dynamics_path);
load('two_d_geometry.mat'); % loads a struct named "two_d_geometry"
spineParameters = two_d_geometry;

%% Spine Parameters

% all parameters hard-coded.
% Let's add one node, hanging from the center, to experiment with positive
% tensions. at z=-h/3 x=0. 
% Let's connect three cables, from now-node 5 to nodes 1, 2, 3.

% Number of bars, cables, and nodes
r = 3; % bars
s = 3; % cables
n = 5; % nodes
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
m_node = m_node(1);
M = spineParameters.total_m; % kg/tetra
springConstant = spineParameters.k_vert;

%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."

% Full connectivity matrix
% Row 1-3 are the cables
% Rows 2-4 are bars (used to be 5-10)
% Columns 1-4 are bottom tetra nodes. Column 5 is the new hanging mass.
%    1  2  3  4  5
C = [1  0  0  0  -1;  %  1
     0  1  0  0  -1;  %  2
     0  0  1  0  -1;  %  3
     1 -1  0  0  0;   %  4
     1  0 -1  0  0;   %  5
     1  0  0 -1  0];  %  6
 
% Checked on 2018-02-08 and confirmed that this matches Jeff Friesen's
% paper's formulation of the connectivity matrix. In (s+r) x n, 10x8.

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

% Combined nodal positions, with the hanging mass now at x=0, z=-h/3
x = [x_bot; 0];
z = [z_bot; -h/3];

% % DEBUGGING
% % Plot nodal positions
% % ...looks good. 
% figure
% plot(x, z,'k.','MarkerSize',10)
% hold on
% return;

%% Lengths of Bars and Cables

% UNUSED RIGHT NOW - change for future, conversion from  min force
% densitites to min cable forces.

% Row 1 is the cable, from node 1 to 5
% Rows 2-4 are bars (was 5-10)
l = [norm([x(1),z(1)]-[x(5),z(5)])  % 1
     ll;                            %  2
     ll;                            %  3
     ls];                           %  4

%% Solve for Reaction Forces
% Assume spine is sitting on a surface. Then there are vertical reaction
% forces at nodes 2 and 3.
  
% Mallory and Ellande version: solve by hand.
% Solve AR*[R2; R3] = bR, where AR will always be invertible
% Now, include the force in z from the hanging mass, for which:
% Sum forces in z needs another g*m_node in bR
% Sum moments needs another g*m_node*(dist from node 2) for node 5
AR = [1 1; ...
      0 (x(3)-x(2))];
  
bR = [M*g + m_node*g; ...
      M*g*(x(1)-x(2)) + m_node*g*(x(5)-x(2))];
  
R = AR\bR;
R2 = R(1);
R3 = R(2);



%% Equilibrium Force Equations - Skelton

% 
disp('Debugging: lets see the Skelton full A matrix, with size:');

A = [ C' * diag(C * x);
      C' * diag(C * z)]
size(A)

% ...this is size s+r = 3, so the size of q will be 3.

%% Skelton's external force vector

% Let's do the forces applied to each node, so we can use the skelton
% configuration.

% We don't need to sum moments here, since nodes are points.
% p is an 10-vector: elements (1:5) are external forces in x, and (6:10) are
% external forces in z. So:

p = zeros(2*n, 1);
% Add the mass to the z-direction for all nodes
p(6:10) = -m_node *g;
% Add the z reaction forces at nodes 2 and 3 (elements 7,8)
p(7) = p(7) + R2;
p(8) = p(8) + R3;

%% Solve Problem for Minimized Cable Tension - Skelton/Friesen, Equality Constraint

% BUT for the skelton formulation, we really do need to
% include +R2 and +R3 in the p vector.

% Here, Drew does not think we can ignore the bars. Later, when we do the
% pseudoinverse, maybe we could - although I still have to test Jeff's
% theory on the matter, and maybe it's easier to just include all of them
% and then put no minimization weight on the bars and no constraints (e.g.
% set H = [eye; zeros] instead of eye, and c = [cs; 0] or something.)

% disp('Using the Skelton formulation, equality:');
% 
% H = zeros(s+r, s+r); % since s=0, this is 3x3, since q \in 3.
% H(1:r, 1:r) = 2 * eye(r); % s+r = 3;
% f = zeros(s+r, 1); % zero vector size 3
% 
% % The force and moment balances are constraints for the optimization
% % A q = p,
% Aeq = A;  % size 16 x 10
% beq = p;  % size 16 x 1
% % ...these work with q size 10 x 1
% 
% % Inequalities ensure a minimum cable tension. With
% %l_cables = l(1:s);
% %L_cables = diag(l_cables);
% % ...the L_cables matrix is diagonal with size s x s. We need to insert it into
% % a zeros matrix with s+r columns now, to align with q of size s+r x 1.
% %Aineq_sk = zeros(s, s+r);
% %Aineq_sk(1:s, 1:s) = -L_cables;
% % We can then keep the original b_inequality, since Aineq_sk * q is of size
% % s x 1. E.g., only need the cable force constraints, and padding with
% % zeros such that the dimensions fit for q.
% %bineq_sk = -minCableTension*ones(s,1);
% 
% disp('Skelton Formulation Solution, Equality Constraint:');
% 
% % let's do the optimization using the skelton formulation. Should be:
% 
% % DEBUGGING: is there a solution without constraining the cables?
% % E.g. is there any solution to A_skelton * q = p?.
% [qOpt_sk, ~, exitFlag] = quadprog(H, f, [], [], Aeq, beq);
% 
% %[qOpt_sk, ~, exitFlag] = quadprog(H_sk, f_sk, Aineq_sk, bineq_sk, Aeq_sk, beq_sk);
% 
% qOpt_sk


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

disp('Relaxing the Skelton formulation:');

% Let's try relaxing the skelton formulation to see if we get a feasible
% answer that way. Maybe it's just that quadprog has difficulty
% initializing? (Or is it really that no solutions exist, and we've
% formulated the problem incorrectly?)
% See the comments below for what each term means.

% ALSO: IMPORANT! If we remove the bars' rows/columns here, then ApA
% becomes = I, V=0, and our problem degenerates. So at least in this case,
% we can't select elements out of A in order to do ApA.
% Conclusion: we will need to do the optimization for ALL tensions, BUT, we
% can zero out the constraints for the bar elements by pre-multiplying the
% "q" with a diagonal matrix of 1s and 0s, see below.

% The pseudoinverse:
A_pinv = pinv(A);

% continuing with Jeff's derivation...
ApA = A_pinv * A;
ApA_dim = size(ApA, 1); % ...so either dimension can be taken.
V = eye(ApA_dim) - ApA;

% Now we can formulate the relaxed problem: (factoring out a 2),
H_lax = V' * V;
f_lax = V' * A_pinv * p;

% Jeff's inequality constraint comes from minimum cable tensions greater
% than 0, which we won't enforce here, since the rods are in compression.

% OBSERVATION: this solves whereas the equality constraint does not!
% Intuition here is that we could probably have used YALMIP and not gone
% through this whole process... but it would be important to discuss in the
% paper that we could relax it and have things be OK.

% Let's put a constraint for positive force density q >= c for tension
% elements, and no constraint for the compression elements.
% That's the same as -q <= -c for q(1).
% We can multiply q by a diagonal matrix to change its sign, where S is the
% full dimension of q, and its upper block is minus the identity matrix of
% size s, and lower block is zeros.
% Let's make it the identity then change the sign for the s-block.

% S multiplies q, so should be size s+r (same as q).
S = zeros(s+r, s+r);
% multiple the first s elements of q by -1
S(1:s, 1:s) = -eye(s);

% c also needs to be [ones(s)*c; 0], so no "other side" for the compression
% elements.
% The homogenous solution seemed to be -1.2 something something for the bars, so try -1
c_tension = 5;
c = [ones(s,1)* c_tension; zeros(r, 1)]; %c is size s+r x 1, or length(q) x 1

% So we have S q <= -c,
% S (A_pinv * p + V w) <= -c, (was A_pinv * p + V w <= -c,)
% S*V*w <= -S * A_pinv * p - c, (was V * w <= -A_pinv * p - c)

A_ineq_lax = S * V;
b_ineq_lax = -S * A_pinv * p - c;

% This seems to work!!! 

% Just to check and confirm that positive tensions would NOT solve, let's
% do V * w => -A_pinv * p, or -V * w <= A_pinv * p
%A_ineq_lax = -V;
%b_ineq_lax = A_pinv * p;

% Correct, this is not feasible.

% Calculate the minumum force densities from the min cable tension and
% length. l_cables is an s-dimensional vector. (capital L is diag version.)
% let's just do element-wise division for ease.
%min_force_densities = minCableTension*ones(s,1) ./ l_cables;
%b_sk_ineq_lax = A_pinv * p_skelton - min_force_densities;

% Finally, let's see if we can solve. Dropping the equality terms.
[wOpt_sk_lax, ~, ~] = quadprog(H_lax, f_lax, A_ineq_lax, b_ineq_lax);

% ...actually, this is the w vector, NOT q. Still need to do q=A_g*F + V*w;
qOpt_sk_lax = A_pinv * p + V * wOpt_sk_lax;

disp('Optimal q, relaxed Skelton formulation:');
qOpt_sk_lax

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

% end of function.