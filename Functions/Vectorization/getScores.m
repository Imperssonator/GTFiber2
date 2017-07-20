function Scores = getScores(ims,j,ep)

MatchEnds = ims.EndLib(j,ep).MatchEnds;
[m n] = size(MatchEnds);

Scores = zeros(m,1);

for i = 1:m
    
    mj = MatchEnds(i,1); mep = MatchEnds(i,2);
    MatchsMatches = ims.EndLib(mj,mep).MatchTableInit(:,3); % The match's matches... only a temp variable, couldn't resist
    curEndInd = sub2ind(size(ims.EndLib),j,ep);
    
    if ismember(curEndInd, MatchsMatches)
        Scores(i,1) = ims.EndLib(mj,mep).MatchTableInit(MatchsMatches==curEndInd,4);
    else
        Scores(i,1) = ScoreMatch(ims,j,ep,mj,mep);
    end
end

end

function Score = ScoreMatch(ims,j,ep,mj,mep)

% Get "connecting vector" and its length (in nanometers)
connVec = ims.EndLib(mj,mep).EVCoord - ims.EndLib(j,ep).EVCoord;        % connVec is in FiberApp xy space (x=j, y=i), same as search vecs
gapLen = norm(connVec)*ims.nmPix;
LenScore = gapLen/ims.settings.stitchGap;

% Calculate Curvature Score
CurvScore = GetCurvScore(ims,j,ep,mj,mep,norm(connVec));

% Weighted Sum of length and curvature scores
Score = LenScore + CurvScore;

% Kill match if gap length or curvature greater than their max
% (negative scores are otherwise impossible)
if (LenScore > 1) || (CurvScore > 1)
    Score = -1;
end
  
end

function CurvScore = GetCurvScore(ims,j,ep,mj,mep,gapLenPix)

% Given the endpoint indices of two segment ends that are candidates for a
% fiber match, fit a fiber, F, to obtain a curvature score for the union
%
% The curvature score will be:
% (Maximum new curvature) / (ims.settings.maxCurv)
%
% Maximum new curvature is obtained by fitting a new fiber using the pixels
% of each pair of segments as an initialization, then evaluating the
% maximum curvature of all points that fall near the location of the stitch
% gap

% Get sorted pixel indices starting with segment j, then concatenate
% segment mj

% if ep==1
%     xy_init = ims.fibSegs(j).xy(:,end:-1:1);
% else
%     xy_init = ims.fibSegs(j).xy;
% end
% 
% if mep==1
%     xy_init = [xy_init,...
%                    ims.fibSegs(mj).xy];
% else
%     xy_init = [xy_init,...
%                    ims.fibSegs(mj).xy(:,end:-1:1)];
% end

if ep==1
    xy_init = flipud(ims.fibSegs(j).sortPixSubs(end:-1:1,:)');
else
    xy_init = flipud(ims.fibSegs(j).sortPixSubs');
end

if mep==1
    xy_init = [xy_init,...
                   flipud(ims.fibSegs(mj).sortPixSubs')];
else
    xy_init = [xy_init,...
                   flipud(ims.fibSegs(mj).sortPixSubs(end:-1:1,:)')];
end

fibImg = int16(ims.gray);

% Fit an active contour
xy = fitContour(fibImg,xy_init,ims.settings.fiberStep);

% Convert to real units
xy_nm = xy .* ims.nmPix;

if size(xy,2)>=3
    dxy=(xy_nm(:,3:end)-xy_nm(:,1:end-2))./(2*ims.settings.fiberStep_nm);
    ddxy=(xy_nm(:,3:end)-2.*xy_nm(:,2:end-1)+xy_nm(:,1:end-2))./(ims.settings.fiberStep_nm)^2;
    curv=(dxy(1,:).*ddxy(2,:)-dxy(2,:).*ddxy(1,:))./sum(dxy.^2,1).^(1.5);
    curv=[0,abs(curv),0];
else
    curv = zeros(1,size(xy,2));
end

% ._._._.  < segment j
% ._._.    < segment mj
% ._._._.-.-._._.   < joined, -.- = gap location
% length(xyj)=4, contour length=3
% length(xymj)=3, contour length = 2
% length(joined)=8, contour length = 7

% gapLoc = 5, which is mean(3,(7-2))+1
% This is gapLoc in terms of indices in joined.xy
% We would want to check the curvature of 4, 5, and 6 at least

% What if previous example joined like this:
% ._._._.--._._.   < joined, -- = gap location
% Now length(joined)=7, contour length = 6
% gapLoc = 4.5, which is mean(3,(6-2))+1
% And this makes sense.
% We would want to check the curvature of 4 and 5 at least

numSegsJoined = length(xy)-1;
numSegsj = length(ims.fibSegs(j).xy)-1;
numSegsmj = length(ims.fibSegs(mj).xy)-1;
gapLoc = mean([...                                  % location of gap as number of vectors into fiber
               numSegsj,...
               (numSegsJoined - numSegsmj)])...
         + 1;

% look at points no less than fiberStep away from gap center (ensures at least 2 pts are checked)
% gapCheck has units of "indices in xy"
gapCheck = max(gapLenPix,ims.settings.fiberStep) / ims.settings.fiberStep;
gapCheckInds = (ceil(gapLoc-gapCheck) : floor(gapLoc+gapCheck));
gapCheckInds(gapCheckInds<1)=1;
gapCheckInds(gapCheckInds>length(xy))=length(xy);

CurvScore = max(curv(gapCheckInds)) / ims.settings.maxCurv;

end

function xy = fitContour(fibImg,xy_init,fiberStep)

% ij_init = column vector of i,j coords of pixels in order
% fibImg = int16 grayscale image
% fiberStep = discretization length in PIXELS

[gradX, gradY] = gradient2Dx2(fibImg);

xy = distributePoints(xy_init,fiberStep);
a = 0;
b = 5;
g = 10;
k1 = 20;
k2 = 10;
fiberIntensity = 255;

% Apply fitting algorithm FA.iterations times
for k = 1:3
    % Construct a matrix M for the internal energy contribution
    n = length(xy);
    
    m1 = [3   2 1:n-2         ;
          2     1:n-1         ;
          1:n           ;
          2:n   n-1     ;
          3:n   n-1 n-2];
    
    m2 = cumsum(ones(5,n), 2);
    
    col = [b; -a-4*b; 2*a+6*b+g; -a-4*b; b];
    m3 = col(:, ones(n,1));
    
    m3(:,1:2) = [b       0        ;
                 -2*b    -a-2*b   ;
                 a+2*b+g 2*a+5*b+g;
                 -a-2*b  -a-4*b   ;
                 b       b       ];
    
    m3(:,end-1:end) = [b         b      ;
                       -a-4*b    -a-2*b ;
                       2*a+5*b+g a+2*b+g;
                       -a-2*b    -2*b   ;
                       0         b     ];
    
    M = sparse(m1, m2, m3, n, n);
    
    % Vector field of external energy
    % (interp2 MATLAB function requests double input, which consume too
    % much memory for big AFM images)
    % gradX and gradY are doubled! (see gradient2Dx2)
    vfx = interp2D(gradX, xy(1,:), xy(2,:), 'cubic')/2;
    vfy = interp2D(gradY, xy(1,:), xy(2,:), 'cubic')/2;   
    vf = k1*[vfx; vfy]/fiberIntensity;
    
    % Normalized vectors of fiber ends
    v_start = xy(:,1) - xy(:,2);
    v_start = v_start/sqrt(sum(v_start.^2));
    v_end = xy(:,end) - xy(:,end-1);
    v_end = v_end/sqrt(sum(v_end.^2));
    
    % Image intensities at the ends
    int_start = interp2D(fibImg, xy(1,1), xy(2,1), 'cubic');
    int_end = interp2D(fibImg, xy(1,end), xy(2,end), 'cubic');
    
    % Full external energy
    vf(:,1) = vf(:,1) + k2*v_start*int_start/fiberIntensity;
    vf(:,end) = vf(:,end) + k2*v_end*int_end/fiberIntensity;
    
    % Next position
    xy = distributePoints((g*xy+vf)/M, fiberStep);
end

end