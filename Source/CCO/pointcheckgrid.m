function [Seg,dcrit] = pointcheckgrid(Seg,RandPoint,K,N,dcrit,Gridx,Gridy)
% POINTCHECKGRID calculates dcrit for terminal node candidate.
%
% Seg - Segment connection matrix
% RandPoint - 200x2 array of random points
% K - Current number of segments in vessel tree
% N - Index of RandPoint to determine terminal node candidate
% dcrit - critical distance of terminal node candidate
% Gridx - GridID of terminal node candidate in x direction
% Gridy - GridID of terminal node candidate in y direction
%
% Seg - Segment connection matrix with new terminal node
% dcrit - new (if any) critical distance
%
% This function returns Dcrit calculated if it is smaller than previous dcrits based on
% distance from end and centre points

    %Find segments that are within 1 grid space of terminal node
    GridArray = find((abs(Seg(:,12)-Gridx) <= 1 | abs(Seg(:,14)-Gridx) <= 1) & (abs(Seg(:,13)-Gridy) <= 1 | abs(Seg(:,15)-Gridy) <= 1));

    % Reset Dmin
    dmin = 1000000;
    
    %Loop over all segments that are within 1 grid space of terminal node
    for i = 1:length(GridArray)

        I = GridArray(i);

        % Calculate distance between ends and centre and various other points

        t(1) = abs(sqrt((RandPoint(N,1) - Seg(I,3))^2 + (RandPoint(N,2) - Seg(I,4))^2));
        t(2) = abs(sqrt((RandPoint(N,1) - Seg(I,1))^2 + (RandPoint(N,2) - Seg(I,2))^2));

        distx = (Seg(I,1)+Seg(I,3))/2;
        disty = (Seg(I,2)+Seg(I,4))/2;
        t(3) = abs(sqrt((RandPoint(N,1) - distx)^2 + (RandPoint(N,2) - disty)^2));

        distx = Seg(I,1) - (Seg(I,1)-Seg(I,3))/4;
        disty = Seg(I,2) - (Seg(I,2)-Seg(I,4))/4;
        t(4) = abs(sqrt((RandPoint(N,1) - distx)^2 + (RandPoint(N,2) - disty)^2));

        distx = Seg(I,1) - (Seg(I,1)+Seg(I,3))*3/4;
        disty = Seg(I,2) - (Seg(I,2)+Seg(I,4))*3/4;
        t(5) = abs(sqrt((RandPoint(N,1) - distx)^2 + (RandPoint(N,2) - disty)^2));

        % Find the smallest distance from all segments
        if min(t) < dmin
            dmin = min(t);
        end
    end

    % If this distance is greater than the previous largest distance update the
    % value
    if dmin > dcrit
        dcrit = dmin;
        Seg(K+2,3) = RandPoint(N,1);
        Seg(K+2,4) = RandPoint(N,2);
    end

end

