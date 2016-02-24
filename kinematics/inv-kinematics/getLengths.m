% This function calculates the distance between points in 3D, passed in as n-vectors
% @author Jeff Friesen
% @param X,Y,Z n-vectors of points in R-3
% @return L a single (n-1) vector of these normed lengths

function [ L ] = getLengths( X,Y,Z )
diffX=diff(X);
diffY=diff(Y);
diffZ=diff(Z);
L=zeros(length(diffX),1);
for i=1:length(diffX)
    L(i)=norm([diffX(i) diffY(i) diffZ(i)],2);
end
end

