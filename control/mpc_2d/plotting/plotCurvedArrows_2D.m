% plotCurvedArrows_2D.m
% adapts the arrows in the visualization of one spine vertebra into 2D,
% by adding arrows for the Z, X axes and \theta rotation.
% Andrew P. Sabelhaus (adapted from Abishek Akella), march 2018.

function [] = plotCurvedArrows_2D()
    hold on
    %CreateCurvedArrow3([0.2 0 -0.03], [0.2 0.03 0], [0.2, 0, 0]);
    %CreateCurvedArrow3([0.2 0.03 0], [0.2, 0, 0.03], [0.2, 0, 0]);
    %CreateCurvedArrow3([0.2 0 0.03], [0.2 -0.03 0], [0.2, 0, 0]);
    
    %CreateCurvedArrow3([0 0.2 -0.03], [0.03 0.2 0], [0, 0.2, 0]);
    %CreateCurvedArrow3([0.03 0.2 0], [0 0.2, 0.03], [0, 0.2, 0]);
    %CreateCurvedArrow3([0 0.2 0.03], [-0.03 0.2 0], [0, 0.2, 0]);
    
    %CreateCurvedArrow3([0, -0.03, .2], [0.03, 0, .2], [0, 0, .2])
    %CreateCurvedArrow3([0.03, 0, .2], [0, 0.03, .2], [0, 0, .2])
    %CreateCurvedArrow3([0, 0.03, .2], [-0.03, 0, .2], [0, 0, .2])
    
    CreateStraightArrow([0 0 0], [0.3 0 0], 'r');
    CreateStraightArrow([0 0 0], [0 0.3 0], 'g');
    CreateStraightArrow([0 0 0], [0 0 0.3], 'k');
    hold off
    %daspect([1 1 1]);
end

function [h] = CreateCurvedArrow3(from, to, center, count)
       
    % Inputs
    if (nargin < 4), count = 15; end
    if (nargin < 3), center = [0 0 0]; end
    center = center(:); from = from(:); to = to(:);

    % Start, stop and normal vectors    
    start = from - center; rstart = norm(start);
    stop = to - center; rstop = norm(stop);
    angle = atan2(norm(cross(start,stop)), dot(start,stop));
    normal = cross(start, stop); normal = normal / norm(normal);

    % Compute intermediate points by rotating 'start' vector
    % toward 'end' vector around 'normal' axis
    % See: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    phiAngles = linspace(0, angle, count);
    r = linspace(rstart, rstop, count) / rstart;
    intermediates = zeros(3, count);
    a = center(1); b = center(2); c = center(3);
    u = normal(1); v = normal(2); w = normal(3); 
    x = from(1); y = from(2); z = from(3);
    for ki = 1:count,
        phi = phiAngles(ki);
        cosp = cos(phi); sinp = sin(phi);
        T = [(u^2+(v^2+w^2)*cosp)  (u*v*(1-cosp)-w*sinp)  (u*w*(1-cosp)+v*sinp) ((a*(v^2+w^2)-u*(b*v+c*w))*(1-cosp)+(b*w-c*v)*sinp); ...
             (u*v*(1-cosp)+w*sinp) (v^2+(u^2+w^2)*cosp)   (v*w*(1-cosp)-u*sinp) ((b*(u^2+w^2)-v*(a*u+c*w))*(1-cosp)+(c*u-a*w)*sinp); ...   
             (u*w*(1-cosp)-v*sinp) (v*w*(1-cosp)+u*sinp)  (w^2+(u^2+v^2)*cosp)  ((c*(u^2+v^2)-w*(a*u+b*v))*(1-cosp)+(a*v-b*u)*sinp); ...
                      0                    0                      0                                1                               ];
        intermediate = T * [x;y;z;r(ki)];
        intermediates(:,ki) = intermediate(1:3);
    end

    % Draw the curved line
    % Can be improved of course with hggroup etc...
    X = intermediates(1,:);
    Y = intermediates(2,:);
    Z = intermediates(3,:);    
    tf = ishold;
    if (~tf), hold on; end
    h = line(X,Y,Z);       
    set(h, 'LineWidth', 3);
    if (~tf), hold off; end
end

function [h] = CreateStraightArrow(from, to, lineSpec)
    h = quiver3(from(1), from(2), from(3), to(1)-from(1), to(2)-from(2), to(3)-from(3), lineSpec);
    set(h, 'LineWidth', 3);
    h.MaxHeadSize = 0.5;
end