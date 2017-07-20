function ims = SegMatch(ims)

%% Image Characteristics
% load(IMS)
w = ims.nmWid;                              % width of image in nm
NumSegs = length(ims.fibSegs);
NumEnds = 2*NumSegs;

%% Hard Coded:
SearchLat = ims.settings.searchLat;             % Now fraction of image width (was 90);
SearchLong = ims.settings.searchLong;           % was 200;

[m,n] = size(ims.segsInit);
nmPix = ims.nmPix;                              % size of a pixel in nm


%% Make Search Kernels

% First build the search kernel, which will look out from each endpoint of
% a segment for other endpoints

Kh = ceil(SearchLat/nmPix);                    % How far to search laterally (from the sides of the segment)
Kw = ceil(SearchLong/nmPix);                   % How far to search longitudinally (out in the direction of the segment)

Kernel = ones(Kh,Kw);                           % The search kernel

% KStart:
% | 0  0  0  0  0  0  0  0  0  0 |              The one is where the endpoint goes
% | 0  0  0  0  0  0  0  0  0  0 |              The 'c' is the center
% | 0  0  1  0  c  0  0  0  0  0 | ---------->  This is the search
% | 0  0  0  0  0  0  0  0  0  0 |              direction
% | 0  0  0  0  0  0  0  0  0  0 |

% In general I apply a two-pixel buffer to the back end of the kernel,
% because sometimes the nearest and best endpoint ends up slightly behind the
% searching endpoint. However, when applying imrotate, the interpolation
% sometimes drops this pixel. Therefore I store the "start pixel" as a
% vector from the center of the image, which is conserved during rotation.

KCent = [Kh/2,Kw/2];
KStart = [Kh/2,5];
KStart_rel = KStart-KCent;  % Vector from center of kernel to start pixel, in units of matrix indices

%% Run The Algo

disp('Running Search...')

for j = 1:NumSegs 
    for ep = 1:2
        
        %% Choosing a Search Angle
        % If it's endpoint 1, Search in the opposite direction of xy(2)-xy(1)
        % If it's endpoint 2, search in the direction of xy(end)-xy(end-1)
        % Keeping in mind that x = j and y = i
    
        if ep==1
            SearchVec = ims.fibSegs(j).xy(:,1) - ims.fibSegs(j).xy(:,2);
            SearchVec = SearchVec./norm(SearchVec);
        else
            SearchVec = ims.fibSegs(j).xy(:,end) - ims.fibSegs(j).xy(:,end-1);
            SearchVec = SearchVec./norm(SearchVec);
        end
        
        SearchAngle = ImVec2RealAng(SearchVec);             % i.e. a search vector of [1, -1] gives a real angle of +45 deg.
        ims.EndLib(j,ep).SearchVec = SearchVec;
        ims.EndLib(j,ep).SearchAngle = SearchAngle;
        
        %% Project Search Boxes to find Match Endpoints
        
        % Rotate kernel and find the subscript indices of the rotated
        % kernel on which it should overlap the searching endpoint
        % (RotKStart)
        RotKern = imrotate(Kernel,SearchAngle);                     % Rotate the kernel to that angle
        RotKStart_rel = [cosd(SearchAngle), -sind(SearchAngle);...  % Rotate the kernel start point to that angle
                         sind(SearchAngle), cosd(SearchAngle)]...
                         * KStart_rel';
        RotKCent = size(RotKern)./2;
        RotKStart = round(RotKStart_rel'+RotKCent);
        
        % Perform the kernel search on the endpoint binary image, starting
        % from the end of the *vectorized* segment
        EPMatch = KernSearch(ims.EndsImg,...
                             RotKern,...
                             RotKStart,...
                             round(ims.EndLib(j,ep).EVCoord(2:-1:1)'));
        
        % Clear out the existing segment from matches to avoid
        % self-matching
        EPMatch_clean = [];
        if all(size(EPMatch))
            for i=1:size(EPMatch,1)
                if ims.SegLabels(EPMatch(i,1),EPMatch(i,2))~=j
                    EPMatch_clean = [EPMatch_clean;...
                        EPMatch(i,:)];
                end
            end
        end
        EPMatch = EPMatch_clean;
        ims.EndLib(j,ep).MatchSubs = EPMatch;               % This is the matrix of subscript indices of endpoints that match with this endpoint

        if not(isempty(EPMatch))
            
            MatchSegs = zeros(size(EPMatch,1),1);           % Since it seems MatchEnds is stored as [seg,ep], this is pretty redundant, but whatever
            MatchEnds = zeros(size(EPMatch,1),2);
            for i = 1:size(EPMatch,1)
                MatchSegs(i) = ims.SegLabels(EPMatch(i,1),EPMatch(i,2));
                MatchEnds(i,:) = ims.Sub2End{EPMatch(i,1),EPMatch(i,2)};
            end
            ims.EndLib(j,ep).MatchSegs = MatchSegs;
            ims.EndLib(j,ep).MatchEnds = MatchEnds;
            
        else
            ims.EndLib(j,ep).MatchSegs = [];
            ims.EndLib(j,ep).MatchEnds = [];
        end
    end
end


%% Ranking Matches

% waitbar(0.7,hwait,'Scoring Matches...')
disp('Ranking Results...')

MatchMatrix1 = zeros(NumEnds);                       % EndIndex = (en-1)*NumEnds+j
MatchMatrix2 = zeros(NumEnds);

% There are basically two relevant angles, between each endpoint and the
% connecting vector (connVec). The cosine of one should be 1, the cosine of
% the other should be -1. Therefore, their product should be as close to -1
% as possible.

[ims.EndLib(:).MatchTableInit]=deal(zeros(0,4));    % Gotta instantiate this field so it can be checked by the first point while getting scores

for j = 1:NumSegs
    
    for ep = 1:2
        if not(isempty(ims.EndLib(j,ep).MatchEnds))
            
            ME_temp = ims.EndLib(j,ep).MatchEnds;   % Too long to type
            MatchLinInds = (ME_temp(:,2)-1)*NumSegs+ME_temp(:,1);
            
            ims.EndLib(j,ep).Scores = getScores(ims,j,ep);      % Get scores for each match - lower is better
            
            ims.EndLib(j,ep).MatchTable = ...
                                    sortrows([ims.EndLib(j,ep).MatchEnds,...
                                              MatchLinInds,...
                                              ims.EndLib(j,ep).Scores...
                                              ], 4);
            
            % Clear matches with scores over the max, but save for
            % debugging
            
            ims.EndLib(j,ep).MatchTableInit = ...
                ims.EndLib(j,ep).MatchTable;
            
            ims.EndLib(j,ep).MatchTable( ...
                ims.EndLib(j,ep).MatchTable(:,4) == -1, :) = [];

            if isempty(ims.EndLib(j,ep).MatchTable)
                continue
            end
            
            % MatchTable:
            % | SegNum, EPNum (1 or 2), EPLinearIndex, Score (ascending order) |
            
            % Enter relevant quantities in match matrices 1 and 2 (1st and
            % 2nd choices)
            
            MatchMatrix1( ims.EndLib(j,ep).EndIndex,...
                          ims.EndLib(j,ep).MatchTable(1,3))...
                          = 1;
                     
            if size(ims.EndLib(j,ep).MatchTable,1)>1
                MatchMatrix2( ims.EndLib(j,ep).EndIndex,...
                              ims.EndLib(j,ep).MatchTable(2,3))...
                              = 1;
            end

        else
            ims.EndLib(j,ep).Scores = [];
            ims.EndLib(j,ep).MatchTable = [];
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

% waitbar(0.8,hwait,'Ranking Matches...')

MM11 = MatchMatrix1.*MatchMatrix1';
Unmatched11 = ~( repmat(sum(MM11,1),size(MM11,1),1) + repmat(sum(MM11,2),1,size(MM11,2)) );

MM12_init = MatchMatrix1.*MatchMatrix2';
MM12_init = MM12_init.*Unmatched11;
[ii jj val] = find(MM12_init);
MM12_subs=[ii, jj, val];
MM12_subs(:,4) = arrayfun(@(mm) ims.EndLib(mm).MatchTable(1,4),MM12_subs(:,1));
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

% waitbar(0.9,hwait,'Percolating Segments...')
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

% waitbar(0.9,hwait,'Labeling Fibers...')
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
        FiberSegs = [FiberSegs ims.EndLib(e).Label];
    end
    F(i).FiberSegs = FiberSegs;
end

% Generate FiberLabels
FiberLabels = zeros(m,n);
for i = 1:length(F)
    FiberSegsi = F(i).FiberSegs;
    NewFiber = MultiEquiv(ims.SegLabels,FiberSegsi).*i;
    FiberLabels = FiberLabels + NewFiber;
end

ims.Fibers = F;
ims.FiberLabels = FiberLabels;

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
% accounting for the edges of the image.

% It returns a list of indices of points that fall within the "field of
% view" of the end of a segment

[Km,Kn] = size(RotKern);
[Im,In] = size(Ends);
RKi = RotKStart(1); RKj = RotKStart(2);
EPi = EP(1); EPj = EP(2);
% RotKern(RKi,RKj) = 0;

Dt = RKi-1; Db = Km-RKi;
Dl = RKj-1; Dr = Kn-RKj;

% Find where search kern hits image edges, if it does
IMt = max(1,EPi-Dt);
IMb = min(Im,EPi+Db);
IMl = max(1,EPj-Dl);
IMr = min(In,EPj+Dr);

% Crop search kern if it hits image edges
Kt = RKi-(EPi-IMt);
Kb = RKi+(IMb-EPi);
Kl = RKj-(EPj-IMl);
Kr = RKj+(IMr-EPj);

% Find endpoints covered by search kern
MatchPts = RotKern(Kt:Kb,Kl:Kr).*Ends(IMt:IMb,IMl:IMr);
[Mati,Matj] = find(MatchPts);
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

% function Scores = getScores(ims,j,ep)
% 
% %% Cosine Scores
% % cosine( searchVec1 , connVec ) * cosine( searchVec2 , connVec )
% 
% searchVec1 = ims.EndLib(j,ep).SearchVec;
% MatchEnds = ims.EndLib(j,ep).MatchEnds;
% [m n] = size(MatchEnds);
% 
% Scores = zeros(m,1);
% 
% for i = 1:m
%     mj = MatchEnds(i,1); mep = MatchEnds(i,2);
%     connVec = ims.EndLib(mj,mep).EVCoord - ims.EndLib(j,ep).EVCoord;        % connVec is in FiberApp xy space (x=j, y=i), same as search vecs
%     searchVec_m = ims.EndLib(mj,mep).SearchVec;
%     Scores(i,1) = dot(searchVec1,connVec)/(norm(searchVec1)*norm(connVec)) *...
%                     dot(connVec,searchVec_m)/(norm(connVec)*norm(searchVec_m)) +...
%                     ( (dot(searchVec1,searchVec_m)/(norm(searchVec1)*norm(searchVec_m))) > cosd(180-ims.settings.maxAngleDeg) )*1000;
% end
% 
% end
