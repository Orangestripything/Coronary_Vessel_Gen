function [r] = CircBound(Aperf,P,MaxPoint)
%Provides the radius of the circle currently
r = sqrt(Aperf*P/(pi()*MaxPoint));
end

