%% Invkin Skelton Example
% Authors: Drew Sabelhaus
% Modified: 3/31/18

% This script tests out the Skelton formulation of Aq=p, hard-coding
% various numbers of vertebrae and cables.

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
% Let's do the two vertebrae, with 4 nodes and 3 bars each.

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
m_node = m_node(1);
M = spineParameters.total_m; % kg/tetra
springConstant = spineParameters.k_vert;

%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."

% Full connectivity matrix
% Row 1-4 are the cables, connecting vertebrae at (top, bottom): 2->2,
% 3->3, 4->2, 4->3. That's 2->6, 3->7, 4->6, 4->7.
% Rows 5-10 are bars 
% Columns 1-4 are bottom tetra nodes, 5-8 are top.
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
     0  0  0  0  1  0  0 -1]; %  10
 
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

% Let's just translate the top one up by 0.1 for now.
% Also do a little bit of x movement.
x_top = x_bot - 0.005;
z_top = z_bot + 0.09;

% Combined nodal positions, both tetra.
x = [x_bot; x_top];
z = [z_bot; z_top];

% % DEBUGGING
% Plot nodal positions
% ...looks good. 
figure
plot(x, z,'k.','MarkerSize',10)
hold on
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
% Now, include the force in z from the top vertebra, for which:
% Sum forces in z needs another g*M in bR
% Sum moments needs another g*M*(dist from node 5)
AR = [1 1; ...
      0 (x(3)-x(2))];
  
bR = [M*g + M*g; ...
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

% ...this is size s+r = 10, so the size of q will be 10 (though we only
% care about the first four, the cables.)

%% Skelton's external force vector

% Let's do the forces applied to each node, so we can use the skelton
% configuration.

% We don't need to sum moments here, since nodes are points.
% p is an 2*n-vector: elements (1:n) are external forces in x, and (n+1:2n) are
% external forces in z. That's p \in 16, first 8 elements, second 8
% elements.

p = zeros(2*n, 1);
% Add the mass to the z-direction for all nodes
p(9:16) = -m_node *g;
% Add the z reaction forces at nodes 2 and 3 (elements n+2, n+3 = 10, 11)
p(10) = p(10) + R2;
p(11) = p(11) + R3;


%% Check the existence of solutions

disp('Checking existence conditions for solutions:');

A * pinv(A)
A * pinv(A) * p
A * pinv(A) * p - p

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
c_tension = 0;
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