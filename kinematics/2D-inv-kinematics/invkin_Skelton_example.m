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

% Number of bars, cables, and nodes
r = 3; % bars
s = 0; % cables
n = 4; % nodes
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
m_node = spineParameters.m; % kg/node
M = spineParameters.total_m; % kg/tetra
springConstant = spineParameters.k_vert;

%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."

% Full connectivity matrix
% Rows 1-3 are bars (used to be 5-10)
% Columns 1-4 are bottom tetra nodes
%    1  2  3  4   
C = [1 -1  0  0;  %  1
     1  0 -1  0;  %  2
     1  0  0 -1];  %  3
 
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

% Combined nodal positions
% As per Jeff's formulation, each of these should be vectors in R^n, where
% n is the number of nodes. So, we're looking at each as an 8-vector.
x = [x_bot];
z = [z_bot];

% DEBUGGING
% Plot nodal positions
% ...looks good. 
% figure
% plot(x_bot,z_bot,'k.','MarkerSize',10)
% hold on
% return;

%% Lengths of Bars and Cables

% Rows 1-3 are bars (was 5-10)
l = [ll;                            %  1
     ll;                            %  2
     ls];                           %  3

%% Solve for Reaction Forces
% Assume spine is sitting on a surface. Then there are vertical reaction
% forces at nodes 2 and 3.
  
% Mallory and Ellande version: solve by hand.
% Solve AR*[R2; R3] = bR, where AR will always be invertible
AR = [1 1; 0 (x(3)-x(2))];
bR = [M*g; M*g*(x(1)-x(2))];
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
% p is an 8-vector: elements (1:4) are external forces in x, and (5:8) are
% external forces in z. So:

p = zeros(2*n, 1);
% Add the mass to the z-direction for all nodes
p(5:8) = -m_node *g;
% Add the z reaction forces at nodes 2 and 3 (elements 6,7)
p(6) = p(6) + R2;
p(7) = p(7) + R3;


%% Solve Problem for Minimized Cable Tension - Skelton/Friesen, Equality Constraint

% BUT for the skelton formulation, we really do need to
% include +R2 and +R3 in the p vector.

% Here, Drew does not think we can ignore the bars. Later, when we do the
% pseudoinverse, maybe we could - although I still have to test Jeff's
% theory on the matter, and maybe it's easier to just include all of them
% and then put no minimization weight on the bars and no constraints (e.g.
% set H = [eye; zeros] instead of eye, and c = [cs; 0] or something.)

disp('Using the Skelton formulation, equality:');

H = zeros(s+r, s+r); % since s=0, this is 3x3, since q \in 3.
H(1:r, 1:r) = 2 * eye(r); % s+r = 3;
f = zeros(s+r, 1); % zero vector size 3

% The force and moment balances are constraints for the optimization
% A q = p,
Aeq = A;  % size 16 x 10
beq = p;  % size 16 x 1
% ...these work with q size 10 x 1

% Inequalities ensure a minimum cable tension. With
%l_cables = l(1:s);
%L_cables = diag(l_cables);
% ...the L_cables matrix is diagonal with size s x s. We need to insert it into
% a zeros matrix with s+r columns now, to align with q of size s+r x 1.
%Aineq_sk = zeros(s, s+r);
%Aineq_sk(1:s, 1:s) = -L_cables;
% We can then keep the original b_inequality, since Aineq_sk * q is of size
% s x 1. E.g., only need the cable force constraints, and padding with
% zeros such that the dimensions fit for q.
%bineq_sk = -minCableTension*ones(s,1);

disp('Skelton Formulation Solution, Equality Constraint:');

% let's do the optimization using the skelton formulation. Should be:

% DEBUGGING: is there a solution without constraining the cables?
% E.g. is there any solution to A_skelton * q = p?.
[qOpt_sk, ~, exitFlag] = quadprog(H, f, [], [], Aeq, beq);

%[qOpt_sk, ~, exitFlag] = quadprog(H_sk, f_sk, Aineq_sk, bineq_sk, Aeq_sk, beq_sk);

qOpt_sk


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

% Let's put a constraint that q <= -c, bars in compression, force density greater
% than c.
% That's A_pinv * p + V w <= -c,
% The inequality constraint is V * w <= -A_pinv * p - c
% The homogenous solution seemed to be -1.2 something something, so try -1
c = -1;
A_ineq_lax = V;
b_ineq_lax = -A_pinv * p + c;

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
% ...this is EXACTLY THE SAME AS THE EQUALITY CONSTRAINT SOLUTION!!!

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