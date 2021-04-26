function [SegArray] = SegDistanceArray(Seg,K)
% SEGDISTANCEARRAY finds the closest 25 segments of a terminal node
%
% Seg - Segment connection matrix
% K - Current number of segments in vessel tree
%
% SegArray - 25x1 array of segment IDs that are the 25 closest neighbours
% of the terminal node
%
% This function loops through all segments and calculates end point and mid
% point distances. For each segment, it stores the minimum distance and
% finally sorts it into ascending order before pruning down to the first
% 25.
    
    %Loop through all segments
    for I = 1 : K
            
            %Calculate various distances
            t(1) = abs(sqrt((Seg(K+2,3) - Seg(I,3))^2 + (Seg(K+2,4) - Seg(I,4))^2));
            t(2) = abs(sqrt((Seg(K+2,3) - Seg(I,1))^2 + (Seg(K+2,4) - Seg(I,2))^2));

            distx = (Seg(I,1)+Seg(I,3))/2;
            disty = (Seg(I,2)+Seg(I,4))/2;
            t(3) = abs(sqrt((Seg(K+2,3) - distx)^2 + (Seg(K+2,4) - disty)^2));

            distx = Seg(I,1) - (Seg(I,1)-Seg(I,3))/4;
            disty = Seg(I,2) - (Seg(I,2)-Seg(I,4))/4;
            t(4) = abs(sqrt((Seg(K+2,3) - distx)^2 + (Seg(K+2,4) - disty)^2));

            distx = Seg(I,1) - (Seg(I,1)+Seg(I,3))*3/4;
            disty = Seg(I,2) - (Seg(I,2)+Seg(I,4))*3/4;
            t(5) = abs(sqrt((Seg(K+2,3) - distx)^2 + (Seg(K+2,4) - disty)^2));
            
            %Stores minimum distance (and seg ID) into matrix d
            d(I,2) = min(t);
            d(I,1) = I;
    end
    
    %Sorts segment distances into ascending order (closest first)
    d = sortrows(d,2);
    
    %Finds the closest 25 segments
    if length(d(:,1)) > 25
        SegArray = d(1:25,1); 
    else
        SegArray = d(:,1);
    end
end

