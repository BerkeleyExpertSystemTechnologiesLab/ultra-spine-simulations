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

% This script uses the Runge-Kutta (RK4) method for temporal discretization
% in order to forward simulate the system's dynamics. The method calculates
% the rate of change of states using the continuous-time dynamics at
% the start, end, and midpoint of each time step and uses the weighted sum to
% approximate the system states at the next time instant. Incorporating
% these additional terms minimizes truncation error relative to the
% first-order solution (Euler's method).
% 
% K1 = f(t, systemStates{t})
% K2 = f(t + dt/2, systemStates{t} + K1*dt/2)
% K3 = f(t + dt/2, systemStates{t} + K2*dt/2)
% K4 = f(t + dt, systemStates{t} + K3*dt)
%
% systemStates{t+1} = systemStates{t} + (dt/6)*(K1 + 2*K2 + 2*K3 + K4)

% Define the magnitude of the noise, for both the position and velocity:
noise_mag_pos = 0.0005; % Was 0.001
noise_mag_vel = 0.0002; % Was 0.0005

tempState{1} = 0; % First link is fixed (Placeholder for state info)
Te{1} = 0;
initial_states = [];
for k = 2:(links+1)
    tempState{k} = systemStates(k-1, :);
    initial_states = [initial_states; tempState{k}'];
end

% Runge-Kutta step 1
% Calculating the rate of change for each state in the system. These
% correspond to linear/angular velocities for the first six states and
% their corresponding accelerations for the second six states.

for k = 2:(links+1)
    Te{k} = getTensions(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12), ...
        tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12), ...
        tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

% Duct_accel solves for the six accelerations for X, Y, Z and orientation
% variables (Tau, Gamma, Phi). Since the velocities are already encoded in
% the state information (final six states of 12-state vector), K1 is
% subsequently built from the last six states of the state vector and these
% six acceleration components.

K = duct_accel(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

K1 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1; % Indexing
    kk = 6*(k-2)+1;
    K1(kn:(kn+5)) = tempState{k}(7:12);
    K1((kn+6):(kn+11)) = K(kk:(kk+5));
end 


% Runge-Kutta step 2
% Again, solving for rate of change for each state. This time the
% calculation takes place at the midpoint of the interval and uses the
% first-order estimate of the states at this point (systemStates + K1*dt/2).

for k = 2:(links+1)
    Te{k} = getTensions(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12), ...
        tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12), ...
        tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

for k = 2:(links+1)
    n = 12*(k-2)+1; % Indexing
    tempState{k} = tempState{k} + (K1(n:(11+n))*dt/2);
end

% Duct_accel solves for the six accelerations for X, Y, Z and orientation
% variables (Tau, Gamma, Phi). Since the velocities are already encoded in
% the state information (final six states of 12-state vector), K2 is
% subsequently built from the last six states of the state vector and these
% six acceleration components.

K = duct_accel(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

K2 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1;
    kk = 6*(k-2)+1;
    K2(kn:(kn+5)) = tempState{k}(7:12);
    K2((kn+6):(kn+11)) = K(kk:(kk+5));
end 

% Runge-Kutta step 3
% Again, solving for rate of change for each state. Once again the
% calculation takes place at the midpoint of the interval and uses the
% first-order estimate of the states at this point (systemStates + K2*dt/2).

for k = 2:(links+1)
    Te{k} = getTensions(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12), ...
        tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12), ...
        tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

for k = 2:(links+1)
    n = 12*(k-2)+1; % Indexing
    tempState{k} = tempState{k} + (K2(n:(11+n))*dt/2);
end

% Duct_accel solves for the six accelerations for X, Y, Z and orientation
% variables (Tau, Gamma, Phi). Since the velocities are already encoded in
% the state information (final six states of 12-state vector), K3 is
% subsequently built from the last six states of the state vector and these
% six acceleration components.

K = duct_accel(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

K3 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1; % Indexing
    kk = 6*(k-2)+1; % Indexing
    K3(kn:(kn+5)) = tempState{k}(7:12);
    K3((kn+6):(kn+11)) = K(kk:(kk+5));
end 

% Runge-Kutta step 4
% Again, solving for rate of change for each state. This time the
% calculation takes place at the endpoint of the interval and uses the
% first-order estimate of the states at this point (systemStates + K4*dt).

for k = 2:(links+1)
    Te{k} = getTensions(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12), ...
        tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12), ...
        tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12), ...
        inputs(k-1, 1),inputs(k-1, 2),inputs(k-1, 3),inputs(k-1, 4),inputs(k-1, 5),inputs(k-1, 6),inputs(k-1, 7),inputs(k-1, 8), ...
        restLengths(1), restLengths(2), restLengths(3), restLengths(4), restLengths(5), restLengths(6), restLengths(7), restLengths(8), k);
end

for k = 2:(links+1)
    n = 12*(k-2)+1; % Indexing
    tempState{k} = tempState{k} + (K3(n:(11+n))*dt);
end

% Duct_accel solves for the six accelerations for X, Y, Z and orientation
% variables (Tau, Gamma, Phi). Since the velocities are already encoded in
% the state information (final six states of 12-state vector), K3 is
% subsequently built from the last six states of the state vector and these
% six acceleration components.

K = duct_accel(tempState{2}(1),tempState{2}(2),tempState{2}(3),tempState{2}(4),tempState{2}(5),tempState{2}(6),tempState{2}(7),tempState{2}(8),tempState{2}(9),tempState{2}(10),tempState{2}(11),tempState{2}(12),Te{2}(1),Te{2}(2),Te{2}(3),Te{2}(4),Te{2}(5),Te{2}(6),Te{2}(7),Te{2}(8), ...
    tempState{3}(1),tempState{3}(2),tempState{3}(3),tempState{3}(4),tempState{3}(5),tempState{3}(6),tempState{3}(7),tempState{3}(8),tempState{3}(9),tempState{3}(10),tempState{3}(11),tempState{3}(12),Te{3}(1),Te{3}(2),Te{3}(3),Te{3}(4),Te{3}(5),Te{3}(6),Te{3}(7),Te{3}(8), ...
    tempState{4}(1),tempState{4}(2),tempState{4}(3),tempState{4}(4),tempState{4}(5),tempState{4}(6),tempState{4}(7),tempState{4}(8),tempState{4}(9),tempState{4}(10),tempState{4}(11),tempState{4}(12),Te{4}(1),Te{4}(2),Te{4}(3),Te{4}(4),Te{4}(5),Te{4}(6),Te{4}(7),Te{4}(8));

K4 = zeros(1, links*12);
for k = 2:(links+1)
    kn = 12*(k-2)+1; % Indexing
    kk = 6*(k-2)+1; % Indexing
    K4(kn:(kn+5)) = tempState{k}(7:12);
    K4((kn+6):(kn+11)) = K(kk:(kk+5));
end 

% Weighting approximations to determine system states at {t+1}.
combinedState = initial_states + dt/6*(K1+2*K2+2*K3+K4)';

% Add Gaussian noise, if desired.
if (noise)
    for k = 1:links
        % Used to be: systemStates(k, 1) = combinedState(12*(k-1)+1) + 0.001*randn(1); % x 
        systemStates(k, 1) = combinedState(12*(k-1)+1) + noise_mag_pos*randn(1); % x
        systemStates(k, 2) = combinedState(12*(k-1)+2) + noise_mag_pos*randn(1); % y
        systemStates(k, 3) = combinedState(12*(k-1)+3) + noise_mag_pos*randn(1); % z
        systemStates(k, 4) = combinedState(12*(k-1)+4) + noise_mag_pos*randn(1); % T
        systemStates(k, 5) = combinedState(12*(k-1)+5) + noise_mag_pos*randn(1); % G
        systemStates(k, 6) = combinedState(12*(k-1)+6) + noise_mag_pos*randn(1); % P
        % Used to be: systemStates(k, 7) = combinedState(12*(k-1)+7) + 0.0005*randn(1); % dx
        systemStates(k, 7) = combinedState(12*(k-1)+7) + noise_mag_vel*randn(1); % dx
        systemStates(k, 8) = combinedState(12*(k-1)+8) + noise_mag_vel*randn(1); % dy
        systemStates(k, 9) = combinedState(12*(k-1)+9) + noise_mag_vel*randn(1); % dz
        systemStates(k, 10) = combinedState(12*(k-1)+10) + noise_mag_vel*randn(1); % dT
        systemStates(k, 11) = combinedState(12*(k-1)+11) + noise_mag_vel*randn(1); % dG
        systemStates(k, 12) = combinedState(12*(k-1)+12) + noise_mag_vel*randn(1); % dP
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