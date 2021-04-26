function [Int] = LineIntersectDau(Int,line1,line2)
% Int is set to 1 if line intersects, if line does not intersect int stays
% the same

% Checks it intersects with the line based upon their line equations
%line1 = [Seg(Z,1), Seg(Z,2); Seg(Z,3), Seg(Z,4)]
%line2 = [RandBiPoint(1), RandBiPoint(2); X1, Y1]


slope = @(line) (line(2,2) - line(1,2))/(line(2,1) - line(1,1));
m1 = slope(line1);
m2 = slope(line2);

intercept = @(line,m) line(1,2) - m*line(1,1);
b1 = intercept(line1,m1);
b2 = intercept(line2,m2);
xintersect = (b2-b1)/(m1-m2);
yintersect = m1*xintersect + b1;

isPointInside = @(xint,myline) ...
    (xint >= myline(1,1) && xint <= myline(2,1)) || ...
    (xint >= myline(2,1) && xint <= myline(1,1));
inside = isPointInside(xintersect,line1) && ...
         isPointInside(xintersect,line2);
if inside == 1 && line1(1) ~= line2(1) && line1(3) ~= line2(3)
    Int = 1;
end

end
