function [AdjIndex] = AdjIndexGen(K)
%ADJINDEXGEN Computes optimum adjacency grid cell index for intersection
%search.
%
% K - Current number of segments in vessel tree
%
% AdjIndex - Optimum adjacency grid cell index
%
% This function uses an experimentally derived function to compute the
% optimum number of grid spaces to consider for checking for segment
% intersections, after forming a new connection.

%Minimum value is 2
AdjIndex = round(2+sqrt(0.003*K),0);

end

