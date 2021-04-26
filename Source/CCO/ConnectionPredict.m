function [BifPoints] = ConnectionPredict(trainedModelx, trainedModely,Seg,SegArray,K)
% CONNECTIONPREDICT Predicts optimum bifurcation point for a list of
% segment connections.
%
% trainedModelx - ML Regression model to predict bifurcation x coordinate
% trainedModely - ML Regression model to predict bifurcation y coordinate
% Seg - Segment connection matrix
% SegArray - Array of Segment IDs to be considered
% K - Current number of segments in vessel tree
%
% BifPoints - 25x2 array of bifurcation point predictions.
%
% The function uses features Px, D1x, Py, D1y, Radius, l3, OrigSeg to
% predict the optimum bifurcation point of a segment connection.
    
    %Filter to the closest 25 segments of the terminal node
    LclSeg = Seg(SegArray,:);   %Filtered segment matrix
    
    %Convert x and y coordinates to be relative to terminal node point
    LocalSeg(:,1:2) = LclSeg(:,[1 3]) - Seg(K+2,3); %x coordinate
    LocalSeg(:,3:4) = LclSeg(:,[2 4]) - Seg(K+2,4); %y coordinate
    
    %Obtained remaining features
    LocalSeg(:,5) = Seg(SegArray,8);    %Radius
    LocalSeg(:,6) = ((LclSeg(:,1)-LclSeg(:,3)).^2 + (LclSeg(:,2)-LclSeg(:,4)).^2).^0.5; %Length of segment
    LocalSeg(:,7) = SegArray<21;    %Presence of original seed segment (1 = True)
    
    %Convert to table
    T = array2table(LocalSeg,'VariableNames',{'Px','D1x','Py','D1y','Radius','l3','OrigSeg'});
    
    %Predict x and y bifurcation points using pretrained model
    BifPoints(:,1) = trainedModelx.predictFcn(T) + Seg(K+2,3);  %Bifurcation x coordinate
    BifPoints(:,2) = trainedModely.predictFcn(T) + Seg(K+2,4);  %Bifurcation y coordinate
end

