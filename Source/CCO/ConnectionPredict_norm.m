function [BifPoints] = ConnectionPredict_norm(trainedModelx, trainedModely,Seg,SegArray,K)
% CONNECTIONPREDICT_NORM Predicts optimum bifurcation point for a list of
% segment connections, using feature normalisation.
%
% trainedModelx - ML Regression model to predict bifurcation x coordinate
% trainedModely - ML Regression model to predict bifurcation y coordinate
% Seg - Segment connection matrix
% SegArray - Array of Segment IDs to be considered
% K - Current number of segments in vessel tree
%
% BifPoints - 25x2 array of bifurcation point predictions.
%
% The function uses features Px, D1x, GPx, Py, D1y, GPy, Radius, l3, OrigSeg to
% predict the optimum bifurcation point of a segment connection.

    %Filter to the closest 25 segments of the terminal node
    LclSeg = Seg(SegArray,:);   %Filtered segment matrix
    
    %Convert x and y coordinates to be relative to terminal node point
    LocalSeg(:,1:2) = Seg(SegArray,[1 3])-Seg(K+2,3);   %x coordinate
    LocalSeg(:,4:5) = Seg(SegArray,[2 4])-Seg(K+2,4);   %y coordinate
    
    %Check to see if GP segment is not the root segment
    try
        %Convert x and y coordinates to be relative to proximal point of
        %segment of concern
        LocalSeg(:,3) = Seg(Seg(SegArray,5),1) - Seg(SegArray,1);   %x coordinate
        LocalSeg(:,6) = Seg(Seg(SegArray,5),2) - Seg(SegArray,2);   %y coordinate
    catch
        %Use zero value if it is the root segment
        LocalSeg(:,3) = 0;  
        LocalSeg(:,6) = 0;
    end
    
    %Compute range of area of consideration
    rangex = max(LocalSeg(:,1:2),[],2) - min([min(LocalSeg(:,1:2),[],2),zeros(20,1)],[],2); %max x coordinate range
    rangey = max(LocalSeg(:,3:4),[],2) - min([min(LocalSeg(:,3:4),[],2),zeros(20,1)],[],2); %max y coordinate range
    
    %Normalise connection points
    LocalSeg(:,1:3) = LocalSeg(:,1:3)./rangex;  %x coordinate
    LocalSeg(:,4:6) = LocalSeg(:,4:6)./rangey;  %y coordinate
    
    %Obtained remaining features
    LocalSeg(:,7) = Seg(SegArray,8);    %Radius
    LocalSeg(:,8) = ((LclSeg(:,1)-LclSeg(:,3)).^2 + (LclSeg(:,2)-LclSeg(:,4)).^2).^0.5; %Length of segment
    LocalSeg(:,9) = SegArray<21;    %Presence of original seed segment (1 = True)
    
    %Convert to table
    T = array2table(LocalSeg,'VariableNames',{'Px','D1x','GPx','Py','D1y','GPy','Radius','l3','OrigSeg'});
    
    %Predict x and y bifurcation points using pretrained model
    BifPoints(:,1) = trainedModelx.predictFcn(T).*rangex + Seg(K+2,3);  %Bifurcation x coordinate
    BifPoints(:,2) = trainedModely.predictFcn(T).*rangey + Seg(K+2,4);  %Bifurcation y coordinate
end

