function ims = orcorr(ims)

% Generate lists of segment vectors and one of their endpoints
fiberXY = {ims.Fibers(:).xy_nm};

vecs = cellfun(@(pts) diff(pts,1,2),fiberXY,'UniformOutput',false);
vecPts = cellfun(@(pts) pts(:,1:end-1),fiberXY,'UniformOutput',false);

vecList = [vecs{:}];
ptList = [vecPts{:}];

N=size(vecList,2);

% Number of bins equal to image width in pixels because why not
numBins = ims.pixWid;

% Maximum distance considered should be the image diagonal
maxDist = norm(size(ims.gray)*ims.nmPix);

% hh is a histogram with rows:
% 1: number of vector pairs at distance r
% 2: <cos(2*theta)> for all vector pairs at distance r
% 3: distance r
hh = zeros(3,numBins);

% pre-compute product of vector magnitudes for dot product
mag = ims.settings.fiberStep_nm^2;

% Iterate over all vector pairs, excluding self-pairs
hwait=waitbar(0,'Calculating Orientation Correlation Function...');
for i = 1:N-1
%     disp(N-i)
    waitbar(i/N,hwait)
    for j = i+1:N
        
        % Compute distance
        dd = sqrt(sum((ptList(:,i)-ptList(:,j)).^2));
        
        % Convert to bin number
        bb = ceil(dd/maxDist*numBins);
        
        % Increment histogram count
        hh(1,bb)=hh(1,bb)+1;
        
        % Compute rolling average of angle cosine
        hh(2,bb)= hh(1,bb)/(hh(1,bb)+1)*hh(2,bb)...
                  + 1/(hh(1,bb)+1)...
                    * cos(2*acos((vecList(1,i)*vecList(1,j)+vecList(2,i)*vecList(2,j))/mag));
    end
end

for i = 1:numBins
    hh(3,i)=i*maxDist/numBins;
end

maxPlotR = ims.nmWid/2;
maxPlotBin = ceil(maxPlotR*numBins/maxDist);

[~, minPlotBin] = max(hh(2,:),[],2);
R = hh(3,minPlotBin:maxPlotBin);
Y = hh(2,minPlotBin:maxPlotBin);

% Calculate correlation length
corr_len = lsqcurvefit(@(b,x) exp(-x./b),ims.nmWid/5,R,Y);

ims.OrCorr = struct();
ims.OrCorr.hist = hh;
ims.OrCorr.corrLen = corr_len;

end