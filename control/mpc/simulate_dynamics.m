% Abishek Akella, Jeff Friesen

function systemStates = simulate_dynamics(systemStates, restLengths, inputs, dt, links, noise)
% Simulates dynamics from the system states for each link and the rest
% lengths of the cables. Requires those as inputs in addition to simulation
% time dt, and number of links (minus the base fixed link). Input noise (1
% or 0) determines whether or not Gaussian noise is added to all states.
%
% Needs to be manually modified when more links are being added - functions
% getTension and duct_accel require state inputs from all links
%
% This script uses Runge-Kutta. The 'K' constants are as per that method.

stateX{1} = 0; % First link is fixed, not important
Te{1} = 0;
initial_states = [];
for k = 2:(links+1)
    stateX{k} = systemStates(k-1, :);
    initial_states = [initial_states; stateX{k}'];
end

% Runge-Kutta step 1
% For the current state, get all the forces on the rigid bodies (via
% tensions)
for k = 2:(links+1)
    Te{k} = getTensions(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12), ...
        stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12), ...
        stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

% Calculate the accelerations applied to the rigid bodies (this is where
% the dynamics are called)
K = duct_accel(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

% Reformat the vectors to incorporate the velocities and accelerations into
% one
K1 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1;
    kk = 6*(k-2)+1;
    K1(kn:(kn+5)) = stateX{k}(7:12);
    K1((kn+6):(kn+11)) = K(kk:(kk+5));
end 


% Runge-Kutta step 2
% For the current state, get all the forces on the rigid bodies (via
% tensions)
for k = 2:(links+1)
    Te{k} = getTensions(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12), ...
        stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12), ...
        stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

% Apply the accelerations to the current state
for k = 2:(links+1)
    n = 12*(k-2)+1;
    stateX{k} = stateX{k} + (K1(n:(11+n))*dt/2);
end

% Calculate the accelerations applied to the rigid bodies (this is where
% the dynamics are called)
K = duct_accel(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

K2 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1;
    kk = 6*(k-2)+1;
    K2(kn:(kn+5)) = stateX{k}(7:12);
    K2((kn+6):(kn+11)) = K(kk:(kk+5));
end 
%
for k = 2:(links+1)
    Te{k} = getTensions(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12), ...
        stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12), ...
        stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

for k = 2:(links+1)
    n = 12*(k-2)+1;
    stateX{k} = stateX{k} + (K2(n:(11+n))*dt/2);
end

K = duct_accel(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

K3 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1;
    kk = 6*(k-2)+1;
    K3(kn:(kn+5)) = stateX{k}(7:12);
    K3((kn+6):(kn+11)) = K(kk:(kk+5));
end 

for k = 2:(links+1)
    Te{k} = getTensions(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12), ...
        stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12), ...
        stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

for k = 2:(links+1)
    n = 12*(k-2)+1;
    stateX{k} = stateX{k} + (K3(n:(11+n))*dt);
end

K = duct_accel(stateX{2}(1),stateX{2}(2),stateX{2}(3),stateX{2}(4),stateX{2}(5),stateX{2}(6),stateX{2}(7),stateX{2}(8),stateX{2}(9),stateX{2}(10),stateX{2}(11),stateX{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    stateX{3}(1),stateX{3}(2),stateX{3}(3),stateX{3}(4),stateX{3}(5),stateX{3}(6),stateX{3}(7),stateX{3}(8),stateX{3}(9),stateX{3}(10),stateX{3}(11),stateX{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    stateX{4}(1),stateX{4}(2),stateX{4}(3),stateX{4}(4),stateX{4}(5),stateX{4}(6),stateX{4}(7),stateX{4}(8),stateX{4}(9),stateX{4}(10),stateX{4}(11),stateX{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

K4 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1;
    kk = 6*(k-2)+1;
    K4(kn:(kn+5)) = stateX{k}(7:12);
    K4((kn+6):(kn+11)) = K(kk:(kk+5));
end 

combinedState = initial_states + dt/6*(K1+2*K2+2*K3+K4)';

if (noise)
    for k = 1:links
        systemStates(k, 1) = combinedState(12*(k-1)+1) + 0.001*randn(1); % x
        systemStates(k, 2) = combinedState(12*(k-1)+2) + 0.001*randn(1); % y
        systemStates(k, 3) = combinedState(12*(k-1)+3) + 0.001*randn(1); % z
        systemStates(k, 4) = combinedState(12*(k-1)+4) + 0.001*randn(1); % T
        systemStates(k, 5) = combinedState(12*(k-1)+5) + 0.001*randn(1); % G
        systemStates(k, 6) = combinedState(12*(k-1)+6) + 0.001*randn(1); % P
        systemStates(k, 7) = combinedState(12*(k-1)+7) + 0.0005*randn(1); % dx
        systemStates(k, 8) = combinedState(12*(k-1)+8) + 0.0005*randn(1); % dy
        systemStates(k, 9) = combinedState(12*(k-1)+9) + 0.0005*randn(1); % dz
        systemStates(k, 10) = combinedState(12*(k-1)+10) + 0.0005*randn(1); % dT
        systemStates(k, 11) = combinedState(12*(k-1)+11) + 0.0005*randn(1); % dG
        systemStates(k, 12) = combinedState(12*(k-1)+12) + 0.0005*randn(1); % dP
    end
else
    for k = 1:links
        systemStates(k, 1) = combinedState(12*(k-1)+1);
        systemStates(k, 2) = combinedState(12*(k-1)+2);
        systemStates(k, 3) = combinedState(12*(k-1)+3);
        systemStates(k, 4) = combinedState(12*(k-1)+4);
        systemStates(k, 5) = combinedState(12*(k-1)+5);
        systemStates(k, 6) = combinedState(12*(k-1)+6);
        systemStates(k, 7) = combinedState(12*(k-1)+7);
        systemStates(k, 8) = combinedState(12*(k-1)+8);
        systemStates(k, 9) = combinedState(12*(k-1)+9);
        systemStates(k, 10) = combinedState(12*(k-1)+10);
        systemStates(k, 11) = combinedState(12*(k-1)+11);
        systemStates(k, 12) = combinedState(12*(k-1)+12);
    end
end


end