function [Seg,x_vector,y_vector] = GridGen(K,Seg,r)
% GRIDGEN Generates an optimised grid mesh and assigns grid cell IDs to each
%segment.
%
% K - Current number of segments in vessel tree
% Seg - Segment connection matrix
% r - radius of current perfusion area
%
% Seg - Segment connection matrix with updated grid IDs
% x_vector - grid cell boundaries in x coordinate
% y_vector - grid cell boundaries in y coordinate
%
% The function first finds the optimum grid cell size using an
% experimentally derived function, to compute x_vector and y_vector.
% After scaling the grid cell boundaries according to r, the function loops
% through each segment and assigns grid cell IDs for its proximal and
% distal points.
    
    %User specified constants
    S = 200;    %(250 - 300)
    G = 8;      %(6-8)
    
    %Determine optimum grid cell size
    gridcellsize = round(sqrt(K((K-S)>0)/G),0);
    
    %Initialise grid cell boundary arrays
    x_vector = 0;
    y_vector = 0;
    
    %Ensures S > K
    if ~isempty(gridcellsize)
        
        %Calculate grid cell boundaries
        x_vector = [1:gridcellsize];
        y_vector = [1:gridcellsize];
        
        %Scale grid cell boundaries
        scale = 2/x_vector(end);
        x_vector = x_vector * scale;
        y_vector = y_vector * scale;

        %Loop through all segments
        for i = 1 : length(Seg(:,1))
            %Proximal x and y 
            Seg(i,12) = find(x_vector >= Seg(i,1)/r, 1, 'first');
            Seg(i,13) = find(y_vector >= Seg(i,2)/r, 1, 'first');
            %Distal x and y
            Seg(i,14) = find(x_vector >= Seg(i,3)/r, 1, 'first');
            Seg(i,15) = find(y_vector >= Seg(i,4)/r, 1, 'first'); 
        end
    end
    
end

