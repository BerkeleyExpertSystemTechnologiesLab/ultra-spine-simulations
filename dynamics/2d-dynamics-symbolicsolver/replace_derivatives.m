function [ result ] = replace_derivatives( input, xi, n, debugging )
%replace_derivatives.m
%   Copyright 2016 Andrew P. Sabelhaus
%   Berkeley Emergent Space Tensegrities Lab, UC Berkeley

%   This function replaces the symbolic derivatives
%   that fulldiff calculates, with the corresponding variables
%   inside a state vector xi.
%   It is used when calling fulldiff on some symbolic expression,
%   where dxi1, dxi2, ... are actually other states in xi,
%   such as dxi1 = xi4, dxi2 = xi5, etc.

% This is designed for dynamics equations with repeated rigid bodies,
% where the state vector consists of positions and their derivatives,
% in order.

% Inputs:
%   input, the symbolic expression that should be modified
%   xi, the state vector. A column vector \in R^(something x 1).
%   n, the number of states per rigid body (e.g., how to
%       'divide up' the xi state vector.)
%   debugging, a flag to turn debugging on or off.
% Outputs:
%   result, which is input with its 'd' symbolic variables replaced.

% First, let's calculate how many rigid bodies there are.
N = size(xi,1)/n;

% Though we know that there are only 6 variables per unit here, let's still use 
% the n variable to calculate which states are positions and
% which are velocities, so this could be used for 2D and 3D dynamics in the future.
% This variable should be 3: (when n=6)
velocity_start_offset = n/2;

% Create the 'result' variable.
result = input;

% Iterate through all the units:
for k=1:N
    % Calculate the start and end indices for the positions at this index
    % This will be, for example, 1, 7, 13, ...
    unit_index_start = 1 + (k-1)*n;
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
        
        % Perform the substitution for this field/value pair
        result = subs(result, oldvalue, newvalue);
    end
end

end

