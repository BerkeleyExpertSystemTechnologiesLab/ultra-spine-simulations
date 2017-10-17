function [ q, A, p, tensions] = InvKin( C , x, y, z, forcesZ, momentsX, momentsY, momentsZ, coms, fixed, minCableTension  )
% Inverse Kinematics Solver V1.0
%   By Ellande Tang, Mallory Daly, Drew Sabelhaus
% Last edited: Oct 17, 2017

% Inputs:
%   C = connectivity matrix. Elements are 0, -1, or 1.
%       - Size is (s + r) by n, where s is the number of cables, r is num
%           rods, and n is number of nodes.
%       - See H.-J. Schek's "The Force Density Method for Form Finding and Computation of General Networks."
%       - Note that C needs to be for the whole structure, not just one body.
%   x, y, z = the positions of each of the nodes. 
%       - Each vector is size n by 1, same dimension as columns(C).
%   forcesZ = vector of forces applied to each body.
%       - Size is m by 1, where m is the number of bodies.
%       - In the future, we'll calculate M automatically somehow, since
%           it's related to n, but that's not supported yet.
%       - This is the force in only the Z direction, so it's a scalar.
%   momentsX = same as forces, m by 1.
%   momentsY = same form as momentsX
%   momentsZ = unused? Is this supported yet?
%   coms = centers of mass. This is a vector of indexes into the C matrix.
%       - Size is m by 1
%       - This script assumes that the bodies are SYMMETRIC, so that one
%           of the nodes is the center of mass. Future work: calculate the
%           COM automatically.
%       - Example: for the tetrahedral spine, with 5 nodes, we specify node
%           one as the center node, so that's the center of mass. With two
%           vertebrae, COM nodes are 1 and 6, so coms = [1 6]'
%   fixed = indices of which nodes are "fixed" to the ground.
%       - Size is k by 1, where k is the number of fixed nodes.
%       - This script works by first solving for the reaction forces at
%           these nodes, then solves for the tensions in the rest of the
%           structure. So, for some applications, this is a bit arbitrary,
%           but needs to happen so that the structure doesn't "fall through
%           the floor" when gravity is applied.
%   minCableTension = minimum force in the cables.
%       - If set to 0, some cables will be likely be untensioned. For most
%           engineering designs, we don't want this.
%       - So, minCableTension is a parameter that we control, and you can
%           change this based on the design that's made.

% Outputs:
%   q = vector of force densities
%       - Size is s by 1 (number of cables, column vector.)
%       - Note that these are NOT the tensions, but the force densities.
%           Multiply by cable length (distance between nodes) to get
%           forces.
%   A = matrix used for Sum(Forces)=0 calculation, A q = p.
%       - See the Schek paper or other research papers (Friesen 2014) for
%           examples. A contains all the information about the cable
%           locations and lengths, calculated from the inputs.
%   p = vector of applied forces to each body.
%       - p is calculated by combining the reaction forces (calculated as
%           part of this script) and the external forces passed in to this
%           script (e.g., forcesZ.)
%       - If this inverse kinematics problem has a solution, there should
%           be cable forces q that exactly balance out p, e.g., A q = p.
%   tensions = the actual cable tensions. This is q * L_cables.

% NOTES:
%   External Z moments not currently supported
%   External X,Y forces not currently supported
%   Can accomodate 2D cases, for fixed nodes, use only 2 elements. 
%       (That assumes an X,Z coordinate system with Z being vertical)

% TO-DO:
%   Check against the pseudo-inverse methods from the Friesen paper.
%       (Is that why we get infeasible solutions sometimes?)
%       Answer: YES, this needs to change to solve using the pseudoinverse
%       of A. Also, need to check that we're not constraining rods to be
%       themselves or something like that? Friesen only uses part of A, I
%       think...
%   Include additional constraints on the cable tensions. Example would
%       be that cables along the same "strip" of material have the same
%       tension. E.g, constrain some q to be the same.

%% Quick multi-output splitter function
% ex: [x, y, z] = dealer([1 2 3]);
toCSV = @(in) in{:};
dealer = @(in) toCSV(num2cell(in));

%% Setup And counting

% Debugging statement:
disp('Starting the inverse kinematics solver for tensegrity systems...');

% Create vector of distance differences
dx = C*x;
dy = C*y;
dz = C*z;

% number of nodes
n =  size(C,2); %nodes

bodies = length(coms); %number of separate bodies

% identify non-com nodes and connected cables of each body

fixedBody = zeros(length(fixed),0);

for a = 1:bodies
    % bars connected to a COM
    temp1 = C(~((C(:,coms(a)))==0),:);
    % removing COM index from list
    temp1(:,coms(a)) = 0;
    % nodes of each body sans COM node
    body{a} = find(sum(abs(temp1),1))';
    % Cables atached to each body
    cables{a} = find(sum(abs(C(((C(:,coms(a))) == 0),body{a})),2));
    % Identify which body is fixed
    for b = 1:length(fixed)
        if any(body{a} == fixed(b))
            fixedBody(b) = a;
        end
    end
end

% indices of all cables calculated from connection matrix
allCables = find(sum(abs(C(:,coms)),2) == 0);

s = length(allCables); %number of cables


%% Lengths of Bars and Cables
% Necessary later to calculate cable lengths
l = sqrt(dx.^2+dy.^2+dz.^2);

% Cable diagonal length matrix
l_cables = l(1:s);
L_cables = diag(l_cables);

%% Solve for Reaction Forces
% Assume spine is sitting on a surface. Then there are vertical reaction
% forces at designated nodes

% TO-DO: this is hard-coded for the tetrahedral spine, in either 2D or 3D
% versions. Instead, needs to be more general, to accommodate arbitrary
% fixed bodies. For example, have reaction forces at the center of the
% vertebrae.

[Rforces,RmomentsX,RmomentsY,RmomentsZ] = deal(zeros(bodies,1));
if length(fixed) == 3
    AR = [1 1 1;...
        % X moment
        0 (y(fixed(2))-y(fixed(1))) (y(fixed(3))-y(fixed(1)));...
        % Y moments
        0 (x(fixed(2))-x(fixed(1))) (x(fixed(3))-x(fixed(1)))];
    bR = [-sum(forcesZ);...
        -(sum(momentsX)+sum(forcesZ.*(y(coms)-y(fixed(1)))));...
        -(sum(momentsY)-sum(forcesZ.*(x(coms)-x(fixed(1)))))];
else
    AR = [1 1;...
        % Y moments
        0 -(x(fixed(2))-x(fixed(1)))];
    bR = [-sum(forcesZ);...
        -(sum(momentsY)-sum(forcesZ.*(x(coms)-x(fixed(1)))))];
end
% bR(all((AR == 0)'),:) = [];
% AR(all((AR == 0)'),:) = [];

R = AR\bR;
% bR = [sum of forces; sum of momentsX; sum of momentsY];

% Adds forces to "forces on structure" matrices
for a = 1:length(fixed)
    targetBody = fixedBody(a);
    Rforces(targetBody) = Rforces(targetBody)+R(a);
    RmomentsX(targetBody) = RmomentsX(targetBody) + (y(fixed(a))-y(coms(targetBody)))*R(a);
    RmomentsY(targetBody) = RmomentsY(targetBody) - (x(fixed(a))-x(coms(targetBody)))*R(a);
end


%% Equilibrium Force Equations
% Creates matrices for solution, re-works the old linear algebra to
% directly solve for q. Generally speaking, the original C'QC = f is
% reformulated as
% [[dx dx ... dx] .* C]'*q = f

for a = 1:bodies
    Cprime(a,:) = sum(C(allCables,body{a}),2)'; %extracts only cable and per body component of connectivity matrix
end

% Moment shorthand function, format of inputs is COM node, connected node on body, non-body cable node
Mom =  @(a,b,c) cross([x(b)-x(a),y(b)-y(a),z(b)-z(a)]',[x(c)-x(b),y(c)-y(b),z(c)-z(b)]');

for a = 1:bodies
    %x,y,z equilibrium criteria
    A1(a,:) = Cprime(a,:).*dx(allCables)';
    A2(a,:) = Cprime(a,:).*dy(allCables)';
    A3(a,:) = Cprime(a,:).*dz(allCables)';
    p1(a,1) = 0;
    p2(a,1) = 0;
    p3(a,1) = forcesZ(a) + Rforces(a);
    p4(a,1) = -(momentsX(a) + RmomentsX(a));
    p5(a,1) = -(momentsY(a) + RmomentsY(a));
    p6(a,1) = -(momentsZ(a) + RmomentsZ(a));
    
    % Creates terms associated with moment equilibrium
    % X moment
    A4(a,:) = zeros(1,s);
    % Y moment
    A5(a,:) = zeros(1,s);
    % Z moment
    A6(a,:) = zeros(1,s);
    for b = 1:s % Cable index; goes through every cable
        if any(cables{a} == allCables(b)) % Is the current cable in the current body's list of attached cables?
            connectedNodes = find(C(b,:)); % Pair of nodes associated with current cable
            if any(body{a} == connectedNodes(1)) % Find which of the two nodes is the one connected to the body
                [A4(a,b),A5(a,b),A6(a,b)] = dealer(Mom(coms(a),connectedNodes(1),connectedNodes(2)));
            else
                [A4(a,b),A5(a,b),A6(a,b)] = dealer(Mom(coms(a),connectedNodes(2),connectedNodes(1)));
            end
        end
    end
end
A = [A1;A2;A3;A4;A5;A6];
p = [p1;p2;p3;p4;p5;p6];

%% Solve Problem for Minimized Cable Tension

% Solve with QUADPROG
% min (1/2 x'Hx + f'x)
% Aineq*x < bineq
% Aeq*x = beq
H = 2*eye(s);
f = zeros(s,1);

% Solve only selected body
% body4Check = 7;
% Aeq = A(body4Check:bodies:body4Check+5*bodies,:);
% beq = p(body4Check:bodies:body4Check+5*bodies,:);

Aeq = A;
beq = p;
Aineq = -L_cables;
bineq = -minCableTension*ones(s,1);
% opts = optimoptions(@quadprog,'Display','notify-detailed');
opts = optimoptions(@quadprog,'MaxIterations',2000);
[q, ~, exitFlag] = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],[],opts);

if exitFlag == 1
    tensions = L_cables*q; % N
    %     restLengths = l_cables - tensions/springConstant;
    %     if any(restLengths <= 0)
    %         display('WARNING: One or more rest lengths are negative. Position is not feasible with current spring constant.')
    %     end
else
    display(['Quadprog exit flag: ' num2str(exitFlag)])
    disp('These tensions will not keep the bodies stationary...');
    tensions = L_cables*q; % N
    %tensions = Inf*ones(s,1);
    %     restLengths = -Inf*ones(s,1);
end


%% Graphs Structure based on Inputs

for a = 1:size(C,1)
    b = find(C(a,:));
    ln(a) = line([x(b(1)),x(b(2))],[y(b(1)),y(b(2))],[z(b(1)),z(b(2))]);
    if any(allCables == a)
        set(ln(a),'color',[0 0 1]);
    else
        set(ln(a),'color',[0 0 0]);
    end
end

line(x(coms),y(coms),z(coms),'linestyle','none','marker','o')

set(gca,...
    'DataAspectRatioMode','manual',...
    'DataAspectRatio',[1 1 1],'view',[45 45],...
    'plotBoxAspectRatioMode','manual',...
    'plotBoxAspectRatio',[1 1 1]);

rotate3d on

end

