function IMS = StitchFibers(IMS,settings)

% 1. Skeleton Cleaning

% 2. Label segments and Make EndLib

% 3. Vectorize Segments - once we are here, we do not go back

% 4. Break high curvatures

% 4.5. Remove Non-conforming widths (measure widths)

% 5. Rebuild EndLib

% 6. Match Search, Ranking/Pairing, Percolation

% 7. Vectorize Fibers

% 8. Break high curvatures again

% 9. Return Fibers for Length/Width/Orientational Analysis

hwait = waitbar(0,'Stitching Fibers...');

%% Image Characteristics
% load(IMS)
w = IMS.nmWid;                              % width of image in nm
% S = IMS.skelTrim;                         % The skeleton segments

%% Hard Coded:
SearchLat = 0.02*w;                         % Now fraction of image width (was 90);
SearchLong = .04*w;                         % was 200;
MinLength = settings.maxBranchSizenm;       % was 30;         % Any segment less than 'maxBranchSizenm' nm long will be cleaned out
ODiffTol = 40;                              % Probably going to phase this out for a curvature estimate


%% Skeleton Cleaning - BETA - should incorporate with main filter
waitbar(0.1,hwait,'Cleaning Skeleton...');
skelClose = bwmorph(imclose(IMS.skelTrim,strel('disk',1)),'skeleton',Inf);
S = cleanSkel(skelClose,settings.maxBranchSize);
Sbranch = bwmorph(S,'branchpoints');
BigBranch = imdilate(Sbranch,ones(3));
S_segs = S&~BigBranch;                      % Find branch points, dilate them, remove the dilated points to separate everything

[m,n] = size(S_segs);
pixdim = IMS.nmPix;                         % size of a pixel in nm
pixarea = pixdim^2;
MinSegLen = ceil(MinLength/pixdim);
SFib = bwareaopen(S_segs,MinSegLen,8);      % Let's not worry about those shorties (<30nm) for a bit

% It has come to my attention that some skeletal segments have holes which
% cause problems. Here we remove those segments.

SFib = RemoveHoles(SFib);


%% Building up Utility Variables and Lookup functions
waitbar(0.2,hwait,'Building Utility Variables...');
RP = regionprops(SFib,'Area','Orientation');    % Grab region props for all the segments
SegLabels = bwlabel(SFib,8);                       % Create an image where their labels correspond to the order regionprops found them in
EndsImg = bwmorph(SFib,'endpoints');            % Also create a binary matrix of the segment endpoints
NumSegs = length(RP);
NumEnds = 2*NumSegs;

End2Sub = cell(NumSegs,2);                      % Cell array that converts segment number into subscript indices of its endpoints in image
Sub2End = cell(m,n);                            % Reverse lookup that converts subscript indices into [segNum, endpoint 1 or 2]
EndLib = struct();                              % Structure Array that will store lots of info on endpoints, will be numSegs x 2

for i = 1:NumSegs
    IsoSeg = SegLabels == i;                       % IsoSeg = an isolated segment i.e. only white pixels where this segment is
    SegEnds = IsoSeg.*EndsImg;                  % SegEnds = just the endpoints of this segment
    [Endi Endj] = find(SegEnds);                % Indices of the end of Segend i
    End2Sub{i,1} = [Endi(1),Endj(1)];
    End2Sub{i,2} = [Endi(2),Endj(2)];
    RP(i).Label = i;
end

for s = 1:NumSegs
    for ep = 1:2
        Coords = End2Sub{s,ep};
        EndLib(s,ep).EPCoord = End2Sub{s,ep};
        EndIndex = (ep-1)*NumSegs+s;
        EndLib(s,ep).EndIndex = EndIndex;       % Sometimes linear indices are useful for finding endpoints
        EndLib(s,ep).Label = s;
        Sub2End{Coords(1),Coords(2)} = [s ep];
    end
end

IMS.SegLabels = SegLabels;
IMS.EndLib = EndLib;


%% Get Active Contour Fits of all segments
waitbar(0.3,hwait,'Fitting Fiber Segments...');

IMS = fitAllSegments(IMS,settings);
% Now we will have IMS.fibSegs, which has the vectorized versions of each
% segment


%% Make Search Kernels

waitbar(0.4,hwait,'Running Projection Search...')
% First build the search kernel, which will look out from each endpoint of
% a segment for other endpoints
% It is 90 x 200 nm, a horizontal bar

Kh = ceil(SearchLat/pixdim);                    % How far to search laterally (from the sides of the segment)
Kw = ceil(SearchLong/pixdim);                   % How far to search longitudinally (out in the direction of the segment)

Kernel = ones(Kh,Kw);                           % The search kernel
KStart = zeros(Kh,Kw);                          % The point in the search kernel that goes at the segment endpoint
KStart(floor(Kh/2),3) = 1;
KStart(floor(Kh/2+1),3) = 1;                    % Had to add this in due to pixel loss in rotation

% KStart:
% | 0  0  0  0  0  0  0  0  0  0 |              The one is where the endpoint goes
% | 0  0  0  0  0  0  0  0  0  0 |
% | 0  0  1  0  0  0  0  0  0  0 | ---------->  This is the search
% | 0  0  0  0  0  0  0  0  0  0 |              direction
% | 0  0  0  0  0  0  0  0  0  0 |

% save('FLdebug')


%% Run The Algo

disp('Running Search...')

for j = 1:NumSegs 
    for ep = 1:2
        
        %% Choosing a Search Angle
        % If it's endpoint 1, Search in the opposite direction of xy(2)-xy(1)
        % If it's endpoint 2, search in the direction of xy(end)-xy(end-1)
        % Keeping in mind that x = j and y = i
    
        if ep==1
            SearchVec = IMS.fibSegs(j).xy(:,1) - IMS.fibSegs(j).xy(:,3);
            SearchVec = SearchVec./norm(SearchVec);
        else
            SearchVec = IMS.fibSegs(j).xy(:,end) - IMS.fibSegs(j).xy(:,end-2);
            SearchVec = SearchVec./norm(SearchVec);
        end
        
        SearchAngle = ImVec2RealAng(SearchVec);             % i.e. a search vector of [1, -1] gives a real angle of +45 deg.
        IMS.EndLib(j,ep).SearchVec = SearchVec;
        IMS.EndLib(j,ep).SearchAngle = SearchAngle;
        
        %% Project Search Boxes to find Match Endpoints
        
        RotKern = imrotate(Kernel,SearchAngle);             % Rotate the kernels to that angle
        RotKStart = imrotate(KStart,SearchAngle);           % Rotate the kernel start point to that angle
        
        EPMatch = KernSearch(EndsImg,RotKern,RotKStart,End2Sub{j,ep});
        IMS.EndLib(j,ep).MatchSubs = EPMatch;               % This is the matrix of subscript indices of endpoints that match with this endpoint

        if not(isempty(EPMatch))
            IMS.EndLib(j,ep).MatchSegs = GetMatchSegs(IMS.EndLib(j,ep),SegLabels);
            IMS.EndLib(j,ep).MatchEnds = GetMatchEnds(EPMatch,Sub2End);
        else
            IMS.EndLib(j,ep).MatchSegs = [];
            IMS.EndLib(j,ep).MatchEnds = [];
        end
    end
end


%% Ranking Matches

waitbar(0.7,hwait,'Scoring Matches...')
disp('Ranking Results...')

MatchMatrix1 = sparse(zeros(NumEnds));                       % EndIndex = (en-1)*NumEnds+j
MatchMatrix2 = sparse(zeros(NumEnds));

% There are basically two relevant angles, between each endpoint and the
% connecting vector (connVec). The cosine of one should be 1, the cosine of
% the other should be -1. Therefore, their product should be as close to -1
% as possible.

% There should also be a criteria for the minimum allowable angle between
% either of these. Starting with 45 degrees ('ODiffTol')

% MinScore = -0.5;    % cos(45)*cos(135) = -0.5

for j = 1:NumSegs %1155:1155
%     disp(j)
    for ep = 1:2
        if not(isempty(IMS.EndLib(j,ep).MatchEnds))
            
            ME_temp = IMS.EndLib(j,ep).MatchEnds;   % Too long to type
            MatchLinInds = (ME_temp(:,2)-1)*NumSegs+ME_temp(:,1);
            
            IMS.EndLib(j,ep).coScores = getCoScores(IMS,j,ep);      % Get cosine scores for each match
            IMS.EndLib(j,ep).MatchTable = ...
                sortrows([IMS.EndLib(j,ep).MatchEnds, MatchLinInds, IMS.EndLib(j,ep).coScores],4);
            
            % MatchTable:
            % | SegNum, EPNum (1 or 2), EPLinearIndex, coScore (desc.) |
            
            % Enter relevant quantities in match matrices 1 and 2 (1st and
            % 2nd choices)
            
            MatchMatrix1(IMS.EndLib(j,ep).EndIndex,IMS.EndLib(j,ep).MatchTable(1,3)) = 1;
            if size(IMS.EndLib(j,ep).MatchTable,1)>1
                MatchMatrix2(IMS.EndLib(j,ep).EndIndex,IMS.EndLib(j,ep).MatchTable(2,3)) = 1;
            end

        else
            IMS.EndLib(j,ep).coScores = [];
            IMS.EndLib(j,ep).MatchTable = [];
        end
    end
end


%% Match Maker
% First, grant mutual 1-1 matches. Then, grant mutual 1-2 matches IF
% neither end is already in a 1-1 match. Finally, grant 2-2 matches if
% neither is in a 1-1 or 1-2 match.

% However, 1-2 matches are tricky because they are not mutually exclusive,
% and can cause circular references. So we must rank the 1-2 matches by
% their overall scores, and iteratively add them to the match matrix, while
% excluding their members from being matched again.

waitbar(0.8,hwait,'Ranking Matches...')

MM11 = MatchMatrix1.*MatchMatrix1';
Unmatched11 = ~( repmat(sum(MM11,1),size(MM11,1),1) + repmat(sum(MM11,2),1,size(MM11,2)) );

MM12_init = MatchMatrix1.*MatchMatrix2';
MM12_init = MM12_init.*Unmatched11;
[ii jj val] = find(MM12_init);
MM12_subs=[ii, jj, val];
MM12_subs(:,4) = arrayfun(@(mm) IMS.EndLib(mm).MatchTable(1,4),MM12_subs(:,1));
MM12_rank = sortrows(MM12_subs,4);
MM12 = zeros(size(MM12_init));
for mi = 1:size(MM12_rank,1)
    mii = MM12_rank(mi,1); mij = MM12_rank(mi,2);
    if MM12_init(mii,mij) == 1
        MM12(mii,mij) = 1;
        MM12_init(mii,:) = 0; MM12_init(:,mii) = 0;
        MM12_init(mij,:) = 0; MM12_init(:,mij) = 0;
    end
end
MM12 = MM12+MM12';
Unmatched12 = ~( repmat(sum(MM12,1),size(MM12,1),1) + repmat(sum(MM12,2),1,size(MM12,2)) );

MM22 = MatchMatrix2.*MatchMatrix2';
MM22 = MM22.*(Unmatched11 & Unmatched12);
MM_Final = MM11 + MM12 + MM22;


%% Percolation
% This "stitches" the fiber segments together using a percolation algorithm
% with MM_Final as the graph matrix describing segment connections

% How this is going to work:
% Start with first dead end
% Add sister end to fiber
% If sister end has no match, remove ends from deadends and try the next
% dead end, increment count
% If sister has a match, switch SegJump to 0 and go to next iter
% When SegJump is 0:
% Find Match and add to fiber
% Make SegJump 1
% Continue

waitbar(0.9,hwait,'Percolating Segments...')
disp('Stitching Fibers...')

DeadEnds = find(sum(MM_Final,2)==0);
F = struct();
count = 1;
SegJump = 1;
CurrEnd = DeadEnds(1);
F(count).Fiber = [CurrEnd];
% save('FLDebug')

while not(isempty(DeadEnds))
    if SegJump
        SisterEnd = FindSister(CurrEnd,NumSegs);
        F(count).Fiber = [F(count).Fiber, SisterEnd];
        SisterMatch = FindMatch(SisterEnd,MM_Final);
        if isempty(SisterMatch)
            DeadEnds = DeadEnds(2:end);
            DeadEnds(find(DeadEnds==SisterEnd))=[];
            if isempty(DeadEnds)
                break
            end
            count = count+1;
            F(count).Fiber = DeadEnds(1);
            CurrEnd = DeadEnds(1);
        else
            SegJump = 0;
        end
    else
        F(count).Fiber = [F(count).Fiber SisterMatch];
        CurrEnd = SisterMatch;
        SegJump = 1;
    end
end


%% Get Fiber Labels

waitbar(0.9,hwait,'Labeling Fibers...')
disp('Labeling Fibers...')

% A Fiber is a list of segment end indices.
% We would like to obtain:
% FiberSegs: a list of segment labels in a given fiber, and
% FiberLabels: an image in which the pixels of each fiber are labeled

% Generate FiberSegs
for i = 1:length(F)
    Fiberi = F(i).Fiber;
    FirstEnds = Fiberi(1:2:end);
    FiberSegs = [];
    for e = FirstEnds
        FiberSegs = [FiberSegs IMS.EndLib(e).Label];
    end
    F(i).FiberSegs = FiberSegs;
end

% Generate FiberLabels
FiberLabels = zeros(m,n);
for i = 1:length(F)
    FiberSegsi = F(i).FiberSegs;
    NewFiber = MultiEquiv(SegLabels,FiberSegsi).*i;
    FiberLabels = FiberLabels + NewFiber;
end

IMS.Fibers = F;
IMS.SFib = SFib;
IMS.SegLabels = SegLabels;
IMS.FiberLabels = FiberLabels;

% save('FLDebug')

%% Active Contour Fit of Stitched Fibers

waitbar(0.95, hwait,'Final Fiber Fits')
IMS = fitAllFibers(IMS,settings);
close(hwait)

end

%% Sub Functions

function SisterEnd = FindSister(CurrEnd,NumSegs)

if CurrEnd > NumSegs
    SisterEnd = CurrEnd-NumSegs;
else
    SisterEnd = CurrEnd+NumSegs;
end

end

function SisterMatch = FindMatch(SisterEnd,MatchMatrix)

SisterMatch = find(MatchMatrix(SisterEnd,:));

end

function [Fnew,IncorpNew] = FiberPerc(F,IncorpList,MatchMatrix,count,i,prev)

% Append segment i to F(count).Fiber, follow any matches that don't go
% backwards onward into recursive call. Add segment i to incorplist.

Fnew(count).Fiber = [F(count).Fiber i];
IncorpNew = [IncorpList i];
[I,Matches] = find(MatchMatrix(i,:));

for j = Matches(Matches~=prev)
    [Fadd,IncorpAdd] = FiberPerc(Fnew,IncorpNew,MatchMatrix,count,j,i);
    Fnew(count).Fiber = union(Fnew(count).Fiber,Fadd(count).Fiber);
    IncorpNew = union(IncorpNew,IncorpAdd);
end

end

function Matches = KernSearch(Ends,RotKern,RotKStart,EP)

% This function basically has to take a rotated kernel, its starting point,
% and overlay it correctly with the right portion of the image, while
% accounting for the edges of the image. I have the details written on a
% piece of paper but it would be quite difficult to document this here. It
% works.

% It returns a list of indices of points that fall within the "field of
% view" of the end of a segment

[Km,Kn] = size(RotKern);
[Im,In] = size(Ends);
[RKi RKj] = find(RotKStart,1);
EPi = EP(1); EPj = EP(2);
SearchKern = RotKern-RotKStart;

Dt = RKi-1; Db = Km-RKi;
Dl = RKj-1; Dr = Kn-RKj;

IMt = max(1,EPi-Dt);
IMb = min(Im,EPi+Db);
IMl = max(1,EPj-Dl);
IMr = min(In,EPj+Dr);

Kt = RKi-(EPi-IMt);
Kb = RKi+(IMb-EPi);
Kl = RKj-(EPj-IMl);
Kr = RKj+(IMr-EPj);

MatchPts = SearchKern(Kt:Kb,Kl:Kr).*Ends(IMt:IMb,IMl:IMr);
[Mati,Matj] = find(sparse(MatchPts));
Matchi = Mati+IMt-1;
Matchj = Matj+IMl-1;

Matches = [Matchi, Matchj];

% save('MatchesDebug')

end

function MatchSegs = GetMatchSegs(ELS,SegLabels)

MatchSubs = ELS.MatchSubs;
[m n] = size(MatchSubs);

MatchSegs = zeros(m,1);

for i = 1:m
    MatchSegs(i) = SegLabels(MatchSubs(i,1),MatchSubs(i,2));
end

end

function MatchEnds = GetMatchEnds(EPMatch,Sub2End)

[m n] = size(EPMatch);
MatchEnds = zeros(m,2);

for i = 1:m
%     disp('----')
%     disp(EPMatch(i,:))
    MatchEnds(i,:) = Sub2End{EPMatch(i,1),EPMatch(i,2)};
end

end

function coScores = getCoScores(IMS,j,ep)

%% Cosine Scores
% cosine( searchVec1 , connVec ) * cosine( searchVec2 , connVec )

searchVec1 = IMS.EndLib(j,ep).SearchVec;
MatchEnds = IMS.EndLib(j,ep).MatchEnds;
[m n] = size(MatchEnds);

coScores = zeros(m,1);

for i = 1:m
    mj = MatchEnds(i,1); mep = MatchEnds(i,2);
    connVec = IMS.EndLib(mj,mep).EVCoord - IMS.EndLib(j,ep).EVCoord;        % connVec is in FiberApp xy space (x=j, y=i), same as search vecs
    searchVec_m = IMS.EndLib(mj,mep).SearchVec;
    coScores(i,1) = dot(searchVec1,connVec)/(norm(searchVec1)*norm(connVec)) *...
                    dot(connVec,searchVec_m)/(norm(connVec)*norm(searchVec_m));
end

end

function out = ImVec2RealAng(V)

[out,~] = cart2pol(V(1),-V(2));     % Find that vector's angle
out = out*180/pi;                   % Convert to degrees

end

function out = AngleDiff(A1,A2)

% AngleDiff finds the absolute value of the "angular distance" between two
% angles in degrees. In this function, AngleDiff(-89, 270) = 1, for example.

Reps = [abs(A1-A2),...
        abs(A1+360-A2),...
        abs(A1-(A2+360))];
    
out = min(Reps);

end

