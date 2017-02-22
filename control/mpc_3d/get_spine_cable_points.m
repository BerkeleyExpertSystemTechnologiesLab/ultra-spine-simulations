% get_spine_cable_points.m
% Copyright 2015 Abishek Akella, Andrew P. Sabelhaus
% This function calculates the positions of all the cables of the spine. This script does NOT use the anchor points.

function String_pts = get_spine_cable_points(Tetra, anchor)

% Inputs:
% Tetra: cell array with the coordinates in 3D of the four outer nodes of each spine vertebra (see ultra_spine_mpc for more info)
% anchor: not used here.

% First, get the number of tetrahedra in this spine configuration. (We use the terms vertebra and tetrahedra interchangeably here.)
num_vertebra = size(Tetra, 2);
% But, we are only going to plot cables for the moving tetrahedra, so:
num_moving_tetra = num_vertebra - 1;

% An array to store all the coordinates of the cables
String_pts = [];
% For each of the 3 moving tetrahedra: (TO-DO: CHECK IF THIS WORKS FOR MORE THAN 3 MOVING TETRAS)
for k = 1:num_moving_tetra
   String_pts = [String_pts; ...
                % Each of the nodes of a tetra is a row vector: ex., 4th node is (4,:).
                Tetra{k}(1,:); Tetra{k+1}(1,:); ...
                Tetra{k}(3,:); Tetra{k+1}(3,:); ...
                Tetra{k}(3,:); Tetra{k+1}(2,:); ...
                Tetra{k}(2,:); Tetra{k+1}(2,:); ...
                Tetra{k}(4,:); Tetra{k+1}(4,:); ...
                Tetra{k}(4,:); Tetra{k+1}(1,:)];
end