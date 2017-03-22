function [ q, A, p] = InvKin( C , x, y, z, forcesZ, momentsX, momentsY, momentsZ, coms, fixed, minCableTension  )
% Inverse Kinematics Solver V1.0
%   By Ellande Tang, credit to Mallory Daly
% Last edited: March 7th, 2017

%% Quick multi-output splitter function
% ex: [x, y, z] = dealer([1 2 3]); 
toCSV = @(in) in{:};
dealer = @(in) toCSV(num2cell(in));

%% Setup And counting

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

[Rforces,RmomentsX,RmomentsY,RmomentsZ] = deal(zeros(bodies,1));

% Solve AR*[R1; R2] = bR, where AR will always be invertible

AR = [1 1 1;...
    % X moment
    0 (y(fixed(2))-y(fixed(1))) (y(fixed(3))-y(fixed(1)));...
    % Y moments
    0 (x(fixed(2))-x(fixed(1))) (x(fixed(3))-x(fixed(1)))];
bR = [-sum(forcesZ);...
    -sum(forcesZ.*(y(coms)-y(fixed(1))));...
    -sum(forcesZ.*(x(coms)-x(fixed(1))))];
R = AR\bR;
% bR = [sum of forces; sum of momentsX; sum of momentsY];

% Adds forces to "forces on structure" matrices
for a = 1:length(fixed)
    targetBody = fixedBody(a);
    Rforces(targetBody) = Rforces(targetBody)+R(a);
    RmomentsX(targetBody) = RmomentsX(targetBody) - (y(fixed(a))-y(coms(targetBody)))*R(a);
    RmomentsY(targetBody) = RmomentsY(targetBody) + (x(fixed(a))-x(coms(targetBody)))*R(a);
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
    p4(a,1) = momentsX(a) + RmomentsX(a);
    p5(a,1) = momentsY(a) + RmomentsY(a);
    p6(a,1) = momentsZ(a) + RmomentsZ(a);
    
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
body4Check = 7;
Aeq = A(body4Check:bodies:body4Check+5*bodies,:);
beq = p(body4Check:bodies:body4Check+5*bodies,:);
% 
% Aeq = A;
% beq = p;
Aineq = -L_cables;
bineq = -minCableTension*ones(s,1);
% opts = optimoptions(@quadprog,'Display','notify-detailed');
opts = optimoptions(@quadprog,'MaxIterations',200);
[q, ~, exitFlag] = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],[],opts);

if exitFlag == 1
    tensions = L_cables*q; % N
    %     restLengths = l_cables - tensions/springConstant;
    %     if any(restLengths <= 0)
    %         display('WARNING: One or more rest lengths are negative. Position is not feasible with current spring constant.')
    %     end
else
    display(['Quadprog exit flag: ' num2str(exitFlag)])
    tensions = Inf*ones(s,1);
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

