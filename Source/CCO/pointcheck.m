function [Seg,dcrit] = pointcheck(Seg,RandPoint,K,N,dcrit)
% Returns Dcrit calculated if it is smaller than previous dcrits based on
% distance from end and centre points

% Reset Dmin

dmin = 1000000;
for I = 1:K
    
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

