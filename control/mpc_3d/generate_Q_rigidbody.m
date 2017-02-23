% generate_Q_rigidbody.m
% Copyright 2016 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function generates a weighting matrix for use in MPC/LQR. (e.g., N-many rigid bodies.)

function Q = generate_Q_rigidbody( weights, weighting_ratio, bodies_do_not_track, N )
% Inputs:
%   weights = a column vector of either size 4 or size 12.
%       If 4 weights: per tetrahedra, these would be one each for xyz, tgp, dot xyz, dot tgp.
%       Here, xyz are the longitudinal coordinates, and Theta Gamma Phi are the euler angles/rotations.
%   weighting_ratio = ratio of weighting between rigid bodies, from first to last, 
%       as an adjustment to weights.
%   bodies_do_not_track = cell array of which bodies to zero out. For ex., {} tracks all, {1,2} tracks the third body only, etc.
%   N = number of rigid bodies.
% Outputs:
%   Q = a 12*N by 12*N matrix with the specified weights along the diagonal.

% Check which set of weights were passed in.
% If there are only 4 weights, then assign them block-diagonally in sections of three
% to each of the rigid body states (this is like weighting xyz, angles, dot xyz, dot angles.

% Sanity checks
assert( size(bodies_do_not_track, 2) <= N, 'Error! which_bodies is greater than N. Cannot weight bodies that are nonexistant!');

if( size(weights,1) == 4)
    % Generate block weighting matrix.
    % Inputs are four weights to be assigned to blocks 1-3, 4-6, 7-9, and 10-12.
    Q = zeros( 12*N, 12*N);

    % Create small 3x3 matrices to copy around.
    Q1_3 = weights(1) * eye(3);
    Q4_6 = weights(2) * eye(3);
    Q7_9 = weights(3) * eye(3);
    Q10_12 = weights(4) * eye(3);

    % For each rigid body:
    for i=1:N
        % Assign the 3x3 matrices according to position.
        % Block 1 to 3:
        Q( 12*(i-1) + 1: 12*(i-1) + 3, 12*(i-1) + 1: 12*(i-1) + 3) = Q1_3;
        % Block 4 to 6:
        Q( 12*(i-1) + 4: 12*(i-1) + 6, 12*(i-1) + 4: 12*(i-1) + 6) = Q4_6;
        % Block 7 to 9:
        Q( 12*(i-1) + 7: 12*(i-1) + 9, 12*(i-1) + 7: 12*(i-1) + 9) = Q7_9;
        % Block 10 to 12:
        Q( 12*(i-1) + 10: 12*(i-1) + 12, 12*(i-1) + 10: 12*(i-1) + 12) = Q10_12;
    end
    
elseif( size(weights,1) == 12)
    % Generate with individual weights for all 12 of the states.
    Q = zeros( 12*N, 12*N);

    % Create a 12x12 matrix with each element of 'weights' on the diagonal.
    Q_one_body = diag(weights);

    % For each rigid body:
    for i=1:N
        % Assign the 12x12 block for this position in the larger Q.
        Q( 12*(i-1) + 1: 12*(i-1) + 12, 12*(i-1) + 1: 12*(i-1) + 12 ) = Q_one_body;
    end
    
else
    error('Weights must be column vectors of either size 4 or size 12.');
end

% Apply the weighting ratio to the different rigid bodies.
% This is used, for example, to weight the top vertebra of the ULTRA Spine
% heavier than the lower vertebrae.
% Since, with MPC, these weights get very large as the horizon gets past 10 or so,
% this weighting ratio should be very small.
assert( weighting_ratio <= 4, 'Weighting ratio is greater than 4, optimization will likely fail, choose a smaller ratio.');

% Create a list of how these blocks should be re-weighted
% Blocks (rigid bodies) are weighted linearly from first (no weighting) to last (weighting equal to weighting_ratio.)
weighting_ratios_by_block = linspace(1, weighting_ratio, N);

% Multiply each block by the ratio.
for i=1:N
    % Pick out the current block
    current_block = Q( 12*(i-1) + 1: 12*(i-1) + 12, 12*(i-1) + 1: 12*(i-1) + 12 );
    % Re-weight it
    %current_block = ((i / N) * weighting_ratio) * current_block; 
    current_block = weighting_ratios_by_block(i) * current_block;
    % Stick this block back in place
    Q( 12*(i-1) + 1: 12*(i-1) + 12, 12*(i-1) + 1: 12*(i-1) + 12 ) = current_block;
end

% Finally, set to zero all the blocks/bodies which are not tracked.
for i=1:size(bodies_do_not_track, 2)
    % Zero out this specific block
    % this block is:
    k = bodies_do_not_track{i};
    Q( 12*(k-1) + 1: 12*(k-1) + 12, 12*(k-1) + 1: 12*(k-1) + 12 ) = zeros(12);
end

% End function.

