% First load the data to be processed
% This While loop creates the strahler order for each segment

for i = 1:10
    mydir  = pwd;
    idcs   = strfind(mydir,'/');
    newdir = strcat(mydir(1:idcs(end-2)),'06_Data/RawData/Unmodified/4000seg/Seg_ML/');

    Seg = readmatrix(strcat(newdir,'4000_seg_ML_',num2str(i),'.csv'));

    % Initialise the loop
    Seg(:,12) = 0;
    Increase = 1;
    N = 0;
    MaxSegment = 8001;

    % While the Strahler order is still increasing
    while Increase == 1
        Increase = 0;


        k = Seg(:,12); % k is the Current strahler orders of all segments
        Location = find(k==N); % Find all segments of the current order

        % Loop through all segments of this order
        for K = 1:size(Location)
            if Seg(Location(K),5) ~= 0

                %find the parent and daughter segment ID's
                Parent = Seg(Location(K),5);
                Daughter1 = Seg(Parent,6);
                Daughter2 = Seg(Parent,7);

                % If both Daughters orders are the same increase the parents by 1
                if Seg(Daughter1,12) == Seg(Daughter2,12)
                    Seg(Parent,12) = N +1;
                    Increase = 1;
                else
                    % Otherwise set the parent to be the higher of the daughters
                    if Seg(Daughter1,12) >= Seg(Daughter2,12)
                        Seg(Parent,12) = Seg(Daughter1,12);
                    else
                        Seg(Parent,12) = Seg(Daughter2,12);
                    end
                end
            end
        end
        N= N+1;
    end
    % The root is always the highest order and is not set within this due to
    % the parent of it being 0 causing an error
    Seg(1,12) = N-1;


    % The diamater based method is repeated 5 times which is more than enough
    % for convergences stated in Kassab's method
    for F = 1:5
        % maxBif is the highest strahler order taken from the bifurcation level
        % previously
        MaxBif = max(Seg(:,12));

        % This finds the mean and standard deviation of each order
        for N = 0:MaxBif
            k = Seg(:,12);
            Location = find(k==N);
            Radius = 0;
            for K = 1:size(Location)
                Radius = Radius + Seg(Location(K),9);
                deviationarray(K) = Seg(Location(K),9);
            end
            Radius = Radius/K;
            StrahlerPlot(N+1,1)= N;
            StrahlerPlot(N+1,2) = Radius;
            StrahlerPlot(N+1,3) = std(deviationarray);
        end

        % Set the diameter boundaries of each order
        for N = 0:MaxBif-1
            StrahlerPlot(N+1,4) = (StrahlerPlot(N+2,2) + StrahlerPlot(N+2,3) + StrahlerPlot(N+1,2) - StrahlerPlot(N+1,3))/2;
        end
        % Set the Highest strahler order boundary to be just below the root to
        % restrict it to be exclusively the root segment
        StrahlerPlot(MaxBif,4)=StrahlerPlot(MaxBif,4)*0.99;

        % Loop through each segment and set its new diamter based strahler
        % order
        for K = 1:MaxSegment
            for I = 1:MaxBif
                if Seg(K,9) > StrahlerPlot(I,4);
                    Seg(K,12)= I;
                end
            end
        end

    end

    % Find the mean and standard deviations of each order for the plot this
    % could be removed by repeating the loop one more time without completing
    % the end section 
    MaxBif = max(Seg(:,12));

    for N = 0:MaxBif
        k = Seg(:,12);
        Location = find(k==N);
        Radius = 0;
        for K = 1:size(Location)
            Radius = Radius + Seg(Location(K),9);
            deviationarray(K) = Seg(Location(K),9);
            Size = size(Location);
            StrahlerPlot(N+1,4) = Size(1);
        end
        Radius = Radius/K;
        StrahlerPlot(N+1,1)= N;
        StrahlerPlot(N+1,2) = Radius;
        StrahlerPlot(N+1,3) = std(deviationarray);
    end

    % This finds the legnths mean and standard deviation at each order
    for N = 0:MaxBif
        k = Seg(:,12);
        Location = find(k==N);
        Length = 0;
        for K = 1:size(Location)
            Length = Length + sqrt((Seg(Location(K),1)-Seg(Location(K),3))^2 + (Seg(Location(K),2)-Seg(Location(K),4))^2);
            deviationarray(K) = sqrt((Seg(Location(K),1)-Seg(Location(K),3))^2 + (Seg(Location(K),2)-Seg(Location(K),4))^2);
            Size = size(Location);
            StrahlerPlot(N+1,4) = Size(1);
        end
        Length = Length/K;
        StrahlerPlot(N+1,1)= N;
        StrahlerPlot(N+1,5) = Length;
        StrahlerPlot(N+1,6) = std(deviationarray);
    end

    % Increase the order to begin at the arteriole level
    % Depends on minimum segment diameter
    StrahlerPlot(:,1) = StrahlerPlot(:,1) + 5


    % Scale the radius to be diamater in millimeters
    StrahlerPlot(:,2) = StrahlerPlot(:,2)*2000;

    % Plot the results
    plot(StrahlerPlot(:,1),StrahlerPlot(:,2))
    csvwrite(strcat('StrahlerOrder_sample',num2str(i),'.csv'),StrahlerPlot);
    hold on

    % Load and plot comparison data if needed
    % load KassabLADMorph
    % plot(KasssabLADMorph(:,1),KasssabLADMorph(:,2),'b')
    % plot(KasssabLADMorph(:,1),KasssabLADMorph(:,3),'b')
    % plot(KasssabLADMorph(:,1),KasssabLADMorph(:,4),'b')

    title('Comparison with Morphometric data')
    xlabel('Strahler Order')
    ylabel('Average Segment diameter (mm)')
    clear
end
% legend('Model','Comparison Data')