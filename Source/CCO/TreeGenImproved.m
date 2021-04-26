%SCRIPT: TreeGenImproved generates a two-dimensional coronary artieral
%vessel tree using a modified constrained constructive optimisation.

tic

% Set Starting Conditions
MaxPoints = 500;                    % Total Number of terminal segments
MaxSegment = (MaxPoints*2) + 1;     % Total Number of segments
MaxPoints = MaxPoints+1;            % For Caluclating r
Aperf = pi()*0.05^2;                % Area To be filled
load TestDataPatient                % Loads Patient Data
Seg = TestDataPatient;              % Assigns data to Seg structure
Points = 200;                       % Random Points created for each iteration
BiPoints = 200;                     % Number of bifurcation Points Tested

%Material Parameters
u = 3.6*10^-3;                      % Blood Viscosity
p = ((1.33*10^4)-(7.98*10^3));      % Blood Pressure
Qtot = 8.33*10^-6;                  % Total Flow Rate to LAD
Q = Qtot/ (MaxPoints+1);            % Set the flow for each terminal end
K = 21;                             % Start after the initial 21 Segments

cnt = 1;                            %cnt increases with every successful connection

%Load machine learning regression models used for bifurcation prediction
load('trainedModelx_simple.mat');   %Regression model for bifurcation x coordinate
load('trainedModely_simple.mat');   %Regression model for bifurcation y coordinate

%Initialise arrays to store computation time for PSA and BSA
PntSel = zeros(MaxPoints,1);        %time spent finding a terminal node
BifSel = zeros(MaxPoints,1);        %time spent connecting a terminal node

%Initialise array to store number of failed attempts
atmpts = zeros(MaxPoints,1);     %number of failed attempts at finding a successful terminal node

%Set initial grid resolution to 3x3
vectorsize = 3;

%Calculate grid cell boundaries
x_vector = [1:vectorsize];
y_vector = [1:vectorsize];

%Scale grid cell boundaries
scale = 2/x_vector(end);
x_vector = x_vector * scale;
y_vector = y_vector * scale;

%Set value of K to start search localisation & ML regression
GridStart = 201;    %Value of K to start search localisation
MLStart = 201;      %Value of K to start ML regression

%Calculate and assign grid cell ID for initial seed segments
for i = 1:K
    point = [Seg(i,1) Seg(i,2) Seg(i,3) Seg(i,4)];
    Seg(i,12) = find(x_vector >= point(1), 1, 'first');
    Seg(i,13) = find(y_vector >= point(2), 1, 'first');
    Seg(i,14) = find(x_vector >= point(3), 1, 'first');
    Seg(i,15) = find(y_vector >= point(4), 1, 'first');
end

%Loop until desired number of terminal segments have been generated
while K <= MaxSegment-2             % -2 as it is K + 2

    %Find optimum AdjIndex
    AdjIndex = AdjIndexGen(K);  %Adjacency Gridcell Index

    %Increase number of failed attempts by 1
    atmpts(cnt,1) = atmpts(cnt,1) + 1;

    %Start PSA algorithm timer
    timer_point = tic;

    P = (K-1)/2 + 2;                   % P is number of points currently
    r = CircBound(Aperf,P,MaxPoints);  % Calculate r of current Circle
    Seg(:,1:4) = Seg(:,1:4)*r;         % Scale Segments

    dthresh = sqrt(pi()*(r^2)/P);      % Set Distance threshold
    dcrit = 0;                         % Reset Critical Distance
    N = 1;                             % Start Counter for Points
    RandPoint = r*2*rand(Points,2);    % Set random points

    %Update grid resolution every 200 segment iterations
    if mod((K-1)/2,200) == 0
        [Seg,x_vector,y_vector] = GridGen(K,Seg,r);
    end

    % Continue to loop until a satisfactory point is found
    while N <= Points

        C = sqrt((RandPoint(N,1)-r)^2 + (RandPoint(N,2)-r)^2);% Check it lies within circle

        if C < r
            %Use search localised dcrit calculation if K > GridStart
            if K > GridStart
                
                %Find grid IDs of terminal node
                Gridx = find(x_vector >= RandPoint(N,1)/r, 1, 'first');
                Gridy = find(y_vector >= RandPoint(N,2)/r, 1, 'first');
                
                %Calculate dcrit using search localisation
                [Seg,dcrit] = pointcheckgrid(Seg,RandPoint,K,N,dcrit,Gridx,Gridy);
            else
                %Calculate dcrit
                [Seg,dcrit] = pointcheck(Seg,RandPoint,K,N,dcrit);
            end
        end
        N = N+1;
    end

    %End PSA algorithm timer and add this to PntSel
    PntSel(cnt,1) = PntSel(cnt,1) + toc(timer_point);

    % The new Point has been found and is now going to be connected

    % Loop through every possible segment connection
    % M is the connection currently being tested

    % Set a high Vmin to begin with
    % if no point is found Vmin will still be 1000 at the end
    
    %Find segment IDs of 25 closest segments
    SegArray = SegDistanceArray(Seg,K);

    Vmin = 1000;

    %Srart BSA algorithm timer
    bif_point = tic;
    
    %Predict optimum bifurcation points for all segments in SegArray
    if K > MLStart
        [PredBifPoints] = ConnectionPredict(trainedModelx_simple, trainedModely_simple,Seg,SegArray,K);
    end
    
    %Loop through all segments in SegArray to form potential connections
    for SegM = 1:length(SegArray)
        
        %Set M to segment ID
        M = SegArray(SegM);

        Int = 0;     % Set it to not initially intersect

        % If it is connecting to the original Segments it may only
        % bifurcate within a restricted range

        if M <= 21
            minx = min([Seg(M,1);Seg(M,3)]);
            maxx = max([Seg(M,1);Seg(M,3)]);

            miny = min([Seg(M,2);Seg(M,4)]);
            maxy = max([Seg(M,2);Seg(M,4)]);

            % Murrays law is given a large leeway for these due to the
            % restriction in bifurcation points
            r1murrayscale = 1.5;
            r2murrayscale = 0.5;
        else

            % Regular Murrays Law leniency
             r1murrayscale = 1.2;
             r2murrayscale = 0.8;
            % Finds the area for the new bifurcation points
            minx = min([Seg(M,1);Seg(M,3);Seg(K+2,3)]);
            maxx = max([Seg(M,1);Seg(M,3);Seg(K+2,3)]);

            miny = min([Seg(M,2);Seg(M,4);Seg(K+2,4)]);
            maxy = max([Seg(M,2);Seg(M,4);Seg(K+2,4)]);
        end

        % Generate Random Bifurcation Points within allocated area
        RandBiPoint = rand(BiPoints,2);
        RandBiPoint(:,1) = minx + (maxx - minx)*RandBiPoint(:,1); 
        RandBiPoint(:,2) = miny + (maxy - miny)*RandBiPoint(:,2);
        
        %Generate Random Bifurcation Points within narrowed area based on
        %ML prediction
        if K > MLStart
            
            %Number of points to iterative through is reduced due to ML
            %prediction
            BiPoints = 20;
            
            %reduce area of consideration to +/-10% of allocated area
            lclminx = PredBifPoints(SegM,1) + ((maxx - minx)/10);   %min x coordinate
            lclmaxx = PredBifPoints(SegM,1) - ((maxx - minx)/10);   %max x coordinate

            lclminy = PredBifPoints(SegM,2) + ((maxy - miny)/10);   %min y coordinate
            lclmaxy = PredBifPoints(SegM,2) - ((maxy - miny)/10);   %max y coordinate
            
            %Generate 20 random points within narrowed area
            RandBiPoint = rand(BiPoints,2);                         %20x2 matrix of random numbers
            
            %Scale RandBiPoint x range and y range
            RandBiPoint(:,1) = lclminx + (lclmaxx - lclminx)*RandBiPoint(:,1);
            RandBiPoint(:,2) = lclminy + (lclmaxy - lclminy)*RandBiPoint(:,2);
        end

        % Loop through bifurcation points finding the optimal within
        % constraints
        for  F = 1:BiPoints

            % Check it lies within circle
            C = sqrt((RandBiPoint(F,1)-r)^2 + (RandBiPoint(F,2)-r)^2);

            if C < r

                % Calculate l and r to find total volume of connection
                l1 = sqrt((RandBiPoint(F,1)-Seg(M,1))^2 + (RandBiPoint(F,2)-Seg(M,2))^2);
                l2 = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);
                l3 = sqrt((RandBiPoint(F,1)-Seg(K+2,3))^2 + (RandBiPoint(F,2)-Seg(K+2,4))^2);

                % Calculates using Poiseuille
                r1 = ((Seg(M,8)+1)*Q*l1*8*u/(p*pi()))^(1/4);
                r2 = (Seg(M,8)*Q*l2*8*u/(p*pi()))^(1/4);
                r3 = (Q*l3*8*u/(p*pi()))^(1/4);

                Vnew = pi()*((l1*r1^2)+(l2*r2^2)+(l3*r3^2));

                % If this is a smaller volume than previously
                % found to be optimal

                if Vnew <= Vmin

                    % Find the angle of the bifurcation
                    a = sqrt((Seg(M,3)-Seg(K+2,3))^2 + (Seg(M,4)-Seg(K+2,4))^2);
                    b = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);
                    c = sqrt((RandBiPoint(F,1)-Seg(K+2,3))^2 + (RandBiPoint(F,2)-Seg(K+2,4))^2);
                    Angle = acosd((b^2 + c^2 - a^2)/(2*b*c));

                    if Angle <= 80

                        % Find the angle between the bifuraction and parent
                        a = sqrt((Seg(M,1)-Seg(K+2,3))^2 + (Seg(M,2)-Seg(K+2,4))^2);
                        b = sqrt((RandBiPoint(F,1)-Seg(M,1))^2 + (RandBiPoint(F,2)-Seg(M,2))^2);
                        c = sqrt((RandBiPoint(F,1)-Seg(K+2,3))^2 + (RandBiPoint(F,2)-Seg(K+2,4))^2);
                        Angle = acosd((b^2 + c^2 - a^2)/(2*b*c));

                        if Angle >= 130

                            % Find the angle between the bifuraction and parent
                            a = sqrt((Seg(M,1)-Seg(M,3))^2 + (Seg(M,2)-Seg(M,4))^2);
                            b = sqrt((RandBiPoint(F,1)-Seg(M,1))^2 + (RandBiPoint(F,2)-Seg(M,2))^2);
                            c = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);
                            Angle = acosd((b^2 + c^2 - a^2)/(2*b*c));

                            if Angle >= 130

                                % Check the new bifurcation falls under Murray's law
                                % and do not grow in radius
                                r1murray = (r2^3 + r3^3)^(1/3)*r1murrayscale;
                                r2murray = (r2^3 + r3^3)^(1/3)*r2murrayscale;
                                rdiff = abs(r1 - r1murray);

                                if r1 <= r1murray && r1>= r2murray &&  r1 > r2 && r1 > r3

                                    % Checks the parent Segment and its bifurcation obey Murray's law
                                    ParSeg = Seg(M,5);

                                    % As long as it is not segment 1
                                    if ParSeg ~= 0

                                        % Finds Segment Data
                                        DauSeg1 = Seg(ParSeg,6);
                                        DauSeg2 = Seg(ParSeg,7);
                                        QDauSeg1 = Seg(DauSeg1,8);
                                        QDauSeg2 = Seg(DauSeg2,8);
                                        l1 = sqrt((Seg(ParSeg,1)-Seg(ParSeg,3))^2 + (Seg(ParSeg,2)-Seg(ParSeg,4))^2);
                                        l2 = sqrt((Seg(DauSeg1,1)-Seg(DauSeg1,3))^2 + (Seg(DauSeg1,2)-Seg(DauSeg1,4))^2);
                                        l3 = sqrt((Seg(DauSeg2,1)-Seg(DauSeg2,3))^2 + (Seg(DauSeg2,2)-Seg(DauSeg2,4))^2);

                                        % If it is the new segment being tested the flow increases

                                        if DauSeg1 == M
                                            QDauSeg1 = QDauSeg1+1;
                                            l2 = sqrt((Seg(DauSeg1,1)-RandBiPoint(F,1))^2 + (Seg(DauSeg1,2)-RandBiPoint(F,2))^2);
                                        end
                                        if DauSeg2 == M
                                            QDauSeg2 = QDauSeg2+1;
                                            l3 = sqrt((Seg(DauSeg2,1)-RandBiPoint(F,1))^2 + (Seg(DauSeg2,2)-RandBiPoint(F,2))^2);
                                        end

                                        % Calculate r and Murrays law
                                        r1 = ((Seg(ParSeg,8)+1)*Q*l1*8*u/(p*pi()))^(1/4);
                                        r2 = (QDauSeg1*Q*l2*8*u/(p*pi()))^(1/4);
                                        r3 = (QDauSeg2*Q*l3*8*u/(p*pi()))^(1/4);

                                        r1murray = (r2^3 + r3^3)^(1/3)*r1murrayscale;
                                        r2murray = (r2^3 + r3^3)^(1/3)*r2murrayscale;
                                        rdiff = abs(r1 - r1murray);
                                    end

                                    if  r1 <= r1murray && r1>= r2murray
                                        if r1 > r2 && r1> r3

                                            % Check the other Segment and its daughters for Murray's
                                            DauSeg1 = Seg(M,6);
                                            DauSeg2 = Seg(M,7);

                                            % Provded it is not terminal
                                            if DauSeg1 ~=0
                                                l2 = sqrt((Seg(DauSeg1,1)-Seg(DauSeg1,3))^2 + (Seg(DauSeg1,2)-Seg(DauSeg1,4))^2);
                                                l3 = sqrt((Seg(DauSeg2,1)-Seg(DauSeg2,3))^2 + (Seg(DauSeg2,2)-Seg(DauSeg2,4))^2);
                                                l1 = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);

                                                r1 = ((Seg(M,8)+1)*Q*l1*8*u/(p*pi()))^(1/4);
                                                r2 = (Seg(DauSeg1,8)*Q*l2*8*u/(p*pi()))^(1/4);
                                                r3 = (Seg(DauSeg2,8)*Q*l3*8*u/(p*pi()))^(1/4);

                                                r1murray = (r2^3 + r3^3)^(1/3)*r1murrayscale;
                                                r2murray = (r2^3 + r3^3)^(1/3)*r2murrayscale;
                                                rdiff = abs(r1 - r1murray);
                                            end

                                            if  r1 <= r1murray && r1>= r2murray
                                                if r1 >r2 
                                                    
                                                    %Use search localised
                                                    %intersection check if
                                                    %K > GridStart
                                                    if K > GridStart
                                                        %Find gridcell IDs
                                                        %of bifurcation
                                                        %point
                                                        Gridx = find(x_vector >= RandBiPoint(F,1)/r, 1, 'first');
                                                        Gridy = find(y_vector >= RandBiPoint(F,2)/r, 1, 'first');
                                                        
                                                        %Find segments
                                                        %within +/-
                                                        %AdjIndex
                                                        %gridspaces of
                                                        %bifurcation point
                                                        GridArray = find((abs(Seg(:,12)-Gridx) <= AdjIndex | abs(Seg(:,14)-Gridx) <= AdjIndex) & (abs(Seg(:,13)-Gridy) <= AdjIndex | abs(Seg(:,15)-Gridy) <= AdjIndex));
                                                        IntK = length(GridArray);
                                                    else
                                                        IntK = K;
                                                    end

                                                    % Intersection Check
                                                    for SegZ = 1:IntK
                                                        
                                                        %Use localised
                                                        %area if K >
                                                        %GridStart
                                                        if K > GridStart
                                                           Z = GridArray(SegZ);
                                                        else
                                                           Z = SegZ;
                                                        end

                                                        % Skip the intersection check for the point it connects to
                                                        % Checks if it intersects with any of the previous segments

                                                        if M ~= Z && Int ==0
                                                            % This checks the new daughter segment
                                                            line1 = [Seg(Z,1), Seg(Z,2); Seg(Z,3), Seg(Z,4)];
                                                            line2 = [RandBiPoint(F,1), RandBiPoint(F,2); Seg(K+2,3), Seg(K+2,4)];
                                                            Int = LineIntersectDau(Int,line1,line2);

                                                            % This checks the old daughter segments
                                                            line2(2,1)= Seg(M,3);
                                                            line2(2,2) = Seg(M,4);
                                                            Int = LineIntersectDau(Int,line1,line2);

                                                            % This tests the parent segment
                                                            if Seg(Z,1) ~= Seg(M,1)
                                                                if Seg(Z,3) ~= Seg(M,1)
                                                                    Int = LineIntersectPar(Seg(Z,:),Int,RandBiPoint(F,:),Seg(M,1),Seg(M,2));
                                                                end
                                                            end
                                                        end
                                                    end

                                                    % If it has passed the
                                                    % checks save it
                                                    if Int == 0
                                                        Vmin = Vnew;
                                                        NewBiffx = RandBiPoint(F,1);
                                                        NewBiffy = RandBiPoint(F,2);
                                                        L = M;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    Int =0;

    %End BSA algorithm timer and add this to BifSel
    BifSel(cnt,1) = BifSel(cnt,1) + toc(bif_point);

    % This saves all the new data necessary
    if Vmin ~= 1000 && Int ==0

        %display percentage completion in Command Window
        fprintf('%s\n\r',strcat(num2str(100*K/(MaxSegment-2)),'%'));

        %Enriched Segment data
        SegEnriched(cnt,1) = Seg(L,1);      %Parent x coordinate
        SegEnriched(cnt,2) = Seg(L,2);      %Parent y coordinate
        SegEnriched(cnt,3) = Seg(L,3);      %Sibling x coordinate
        SegEnriched(cnt,4) = Seg(L,4);      %Sibling y coordinate
        SegEnriched(cnt,5) = Seg(K+2,3);    %Terminal node x coordinate
        SegEnriched(cnt,6) = Seg(K+2,4);    %Terminal node y coordinate
        SegEnriched(cnt,7) = Seg(L,8);      %Flow rate of segment
        SegEnriched(cnt,8) = r;             %radius of perfusion area
        SegEnriched(cnt,9) = P;             %total current segment count
        D = find(Seg(:,5)==Seg(L,5),2);     %Daughter IDs of parent segment
        D = D(D~=L);                        %Sibling segment ID
        GP = Seg(L,5);                      %Grandparent segment ID
        SegEnriched(cnt,10) = Seg(D,3);     %Parent 2 x coordinate
        SegEnriched(cnt,11) = Seg(D,4);     %Parent 2 y coordinate
        
        %check if Grandparent segment is not root segment
        if GP ~= 0          
         
            SegEnriched(cnt,12) = Seg(GP,1);    %Grandparent x coordinate
            SegEnriched(cnt,13) = Seg(GP,2);    %Grandparent y coordinate
            
            %Lengths of parent segment, parent 2 segment and grandparent
            %segment
            l1 = sqrt((Seg(GP,1)-Seg(GP,3))^2 + (Seg(GP,2)-Seg(GP,4))^2);
            l2 = sqrt((Seg(D,1)-Seg(D,3))^2 + (Seg(D,2)-Seg(D,4))^2);
            l3 = sqrt((Seg(L,1)-Seg(L,3))^2 + (Seg(L,2)-Seg(L,4))^2);

            %Calculate bifurcation angle beta
            a = sqrt((Seg(GP,1)-Seg(L,3))^2 + (Seg(GP,2)-Seg(L,4))^2);
            b = l1;
            c = l3;
            beta = acosd((b^2 + c^2 - a^2)/(2*b*c));

            %Calculate angle between grandparent and parent segment           
            a = sqrt((Seg(D,3)-Seg(L,3))^2 + (Seg(D,4)-Seg(L,4))^2);
            b = l2;
            c = l3;
            alpha = acosd((b^2 + c^2 - a^2)/(2*b*c));
            
            %Store above variables in SegEnriched
            SegEnriched(cnt,15) = l1;
            SegEnriched(cnt,16) = l2;
            SegEnriched(cnt,17) = l3;
            SegEnriched(cnt,18) = beta;
            SegEnriched(cnt,19) = alpha;
        end
        
        %Find if segment is a seed segment
        if L < 21
            SegEnriched(cnt,14) = 1;    %1 = True
        else
            SegEnriched(cnt,14) = 0;    %0 = False
        end
        
        %Bifurcation Point
        SegEnriched(cnt,20) = NewBiffx; %Bifurcation x coordinate
        SegEnriched(cnt,21) = NewBiffy; %Bifurcation y coordinate

        %Increase segment iteratino counter by 1
        cnt = cnt + 1;

        % Sets the points for the new segments
        Seg(K+1,1) = NewBiffx;
        Seg(K+1,2) = NewBiffy;
        Seg(K+1,3) = Seg(L,3);
        Seg(K+1,4) = Seg(L,4);

        % Changes the parent segments end point
        Seg(L,3) = NewBiffx;
        Seg(L,4) = NewBiffy;

        % Changes the starting point of the segment 3,4 already saved
        Seg(K+2,1) = NewBiffx;
        Seg(K+2,2) = NewBiffy;

        % Assigns the parent to it
        Seg(K+1,5) = L;
        Seg(K+2,5) = L;

        % Assigns parents to old daughters
        Daughter1 = Seg(L,6);
        Daughter2 = Seg(L,7);

        % K + 1 always connecting segment so becomes new parent
        if Daughter1 >= 1
            Seg(Daughter1,5) = K+1;
            Seg(Daughter2,5) = K+1;
        end

        % Assign the daughter to the new connecting segment
        Seg(K+1,6) = Seg(L,6);
        Seg(K+1,7) = Seg(L,7);

        % Assign the new daughters to the parent
        Seg(L,6) = K+1;
        Seg(L,7) = K+2;

        % Assigns flow to the new point and the connection point
        Seg(K+2,8) = 1;
        Seg(K+1,8) = Seg(L,8);

        %Assign grid cell IDs to newly formed segments
        gridcellarray = [L,K+1,K+2];
        for grd = 1:length(gridcellarray)
            g = gridcellarray(grd);
            point = [Seg(g,1) Seg(g,2) Seg(g,3) Seg(g,4)]/r;
            Seg(g,12) = find(x_vector >= point(1), 1, 'first');
            Seg(g,13) = find(y_vector >= point(2), 1, 'first');
            Seg(g,14) = find(x_vector >= point(3), 1, 'first');
            Seg(g,15) = find(y_vector >= point(4), 1, 'first');

        end

        V = 1;

        % Recursively Increases flow up the tree to the root
         while Seg(L,5) >=1
             Seg(L,8) = Seg(L,8) + 1;
             L = Seg(L,5);
         end

        % Root's flow always increased by 1 and is not affected by previous
        Seg(1,8) = Seg(1,8) + 1;

        % Calculates the radius of the new segments
        l = sqrt((Seg(K,1)-Seg(K,3))^2 + (Seg(K,2)-Seg(K,4))^2);
        dp = p/25;
        Seg(K,9)= ((Seg(K,8)*Q*l*8*u/(p*pi()))^(1/4));

        l = sqrt((Seg(K+1,1)-Seg(K+1,3))^2 + (Seg(K+1,2)-Seg(K+1,4))^2);
        Seg(K+1,9)= ((Seg(K+1,8)*Q*l*8*u/(p*pi()))^(1/4));
        K = K + 2; 
    end

    Seg(:,1:4) = Seg(:,1:4)/r;       % Rescaled back to neutral between (0-2)
    NewBiffx= 0;                     % Reset Bifx = 0
    NewBiffy = 0;                    % Reset Bify = 0
end

%End Timer for complete generation process
CompTime = toc;

% For some reason it is one extra at the end
Seg(1,8) = Seg(1,8) - 1;

%Scale Tree to Final size
P = (K-1)/2 + 2;                   % P is number of points created
r = CircBound(Aperf,P,MaxPoints);  % Calculate r of current Circle
Seg(:,1:4) = Seg(:,1:4)*r;
rsize = r;

% Goes through each segment and determines its bifurcation level

for I = 1:K
    V = 0;
    T = I;
    while Seg(T,5) >= 1
        T = Seg(T,5);
        V = V+1;
    end
    Seg(I,11) = V;
end


MaxBif = max(Seg(:,11));

% Calculates r for every segment based on the length and the flow rate
% through it
    VTot = 0;
    for E = 1:MaxSegment
        l = sqrt((Seg(E,1)-Seg(E,3))^2 + (Seg(E,2)-Seg(E,4))^2);
        Seg(E,9) = (Seg(E,8)*Q*l*8*u/((p/MaxBif)*pi()))^(1/4);
        VTot = VTot+ (pi()*(l*Seg(E,9)^2));
    end

% Scales the radius for the plot so that its line width is not x10^-5
Seg(:,10) = Seg(:,9)*10000;

% Plots visual representation of vessels 
figure(1)
title('Generated Coronary Vessel Tree')
for E = 1:MaxSegment
    x1 = [Seg(E,1);Seg(E,3)];
    x2 = [Seg(E,2);Seg(E,4)];
    r = Seg(E,10);
    plot(x1,x2,'black','Linewidth',r)
    hold on
end
% Scales the plot so that it is centred (can be removed)
xlim([0 0.1044])
ylim([0 0.1044])
h = circle(rsize,rsize,rsize);             % Plot circle boundary

% Plots the average radius against the number of bifuractions down the tree

    figure(2)
    title('Average Segment Diameter Along Vessel Tree')
    xlabel('Bifurcation level')
    ylabel('Average Segment diameter (mm)')

    % Calculates the average radius for each level
    for N = 0:MaxBif
        k = Seg(:,11);
        Location = find(k==N);
        Radius = 0;
        for K = 1:size(Location)
            Radius = Radius + Seg(Location(K),9);
        end
        Radius = Radius/K;
        RadiusPlot(N+1,1)= N;
        RadiusPlot(N+1,2) = Radius;
    end

    % Plots the graph and converts to diameter and milimeters
    RadiusPlot(1,1) = 0;
    RadiusPlot(1,2) = Seg(1,9);
    RadiusPlot(:,2) = RadiusPlot(:,2)*2000;
    plot(RadiusPlot(:,1),RadiusPlot(:,2));
    
%Plot graph of Computation Time for both algorithms
figure(3)
plot(BifSel);
hold on
plot(PntSel);
xlabel('Segment Iteration');
ylabel('Time (s)');
yyaxis right
plot(atmpts,'-k');
ylabel('Failed Attempts');
ylim([-30,40]);
xlim([0,MaxPoints])
legend('BifSelectionTime','PointSelectionTime','Failed attempts','FontSize',14);
title(strcat('Computation Time for',{' '},num2str(MaxPoints),' Segment Generation'),'fontweight','normal','FontSize',18);
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

VTot;
