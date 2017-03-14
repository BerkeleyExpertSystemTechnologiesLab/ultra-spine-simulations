clear
close all

load('datav4.mat')

% Geometric parameters
minCableTension = 0;


% Mass and force parameters
g = gravity; % m/s^2, acceleration due to gravity
m = m_rod;

forces = [-m*g -m*g -m*g -m*g -m*g -m*g -m_payload*g]';
momentsX = [0 0 0 0 0 0 0];
momentsY = [0 0 0 0 0 0 0];
momentsZ = [0 0 0 0 0 0 0];
coms = [15 16 17 18 19 20 21]';

% Fixed nodes should both be part of the same body and not the COM
% NOTE: future revisions will improve upon this limitation, current design
% essentially chooses a vertebra to fix in space rather than the more
% physically accurate model of designating particular nodes as resting on a
% surface.
fixed = [1 5 9]';
bodies = length(coms); %number of separate bodies

%% Connectivity Matrix
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks."

% Create vector of distance differences
dx = C*x;
dy = C*y;
dz = C*z;

% number of nodes
n =  size(C,2); %nodes

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
    
    % Identify which body is fixed (limited to one fixed body at the
    % moment)
    for b = 1:length(fixed)
        if any(body{a} == fixed(b))
            fixedBody(b) = a;
        end
    end
    
end

% indices of all cables calculated from connection matrix
allCables = find(sum(abs(C(:,coms)),2) == 0);


s = length(allCables); %number of cables


[ q, A, p ] = InvKin( C , x, y, z, forcesZ, momentsX, momentsY, momentsZ, coms, fixed, minCableTension  );