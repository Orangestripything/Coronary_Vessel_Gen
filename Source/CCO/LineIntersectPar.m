function [Int] = LineIntersectPar(Seg,Int,RandBiPoint,X1,Y1)
% Int is set to 1 if line intersects, if line does not intersect int stays
% the same
% Check it intersects using determinate
% Different functions are used for the parent and daughter as for some
% reason they do not work on each other

% Checks it intersects with the line
x=[Seg(1) Seg(3) RandBiPoint(1) X1];
y=[Seg(2) Seg(4) RandBiPoint(2) Y1];
dt1=det([1,1,1;x(1),x(2),x(3);y(1),y(2),y(3)])*det([1,1,1;x(1),x(2),x(4);y(1),y(2),y(4)]);
dt2=det([1,1,1;x(1),x(3),x(4);y(1),y(3),y(4)])*det([1,1,1;x(2),x(3),x(4);y(2),y(3),y(4)]);

if(dt1<=0 & dt2<=0)
    Int = 1;         %If lines intesect
end
end

