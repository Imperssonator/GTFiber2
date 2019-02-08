function ims = orcorr2d_noVec(ims)

% Get image size in nanometers
pixSize = size(ims.gray);
nmSize = pixSize.*ims.nmPix;

% Generate list of angle-containing pixel locations
[XX,YY] = meshgrid(linspace(0,ims.nmWid,ims.pixWid), linspace(0,ims.nmWid,ims.pixWid));
BW = ims.skelTrim;
AM = ims.AngMap;
ptList = [-YY(BW), XX(BW)]; % Putting this backwards because bin coordinates are [i,j]...
AngList = AM(BW);

% Make bins equal to image size + 1 in x and y directions
% But numBins is half the number of bins (only in positive directions)
numBins = round(size(ims.gray)/2);
db = nmSize./numBins;

% hh is a 2D histogram with bins for each pixel of the image and
% three entries at each bin:
% 1: number of vector pairs in that spatial bin
% 2: <cos(2*theta)> for all vector pairs in that bin
hh = zeros([2.*numBins+1, 2]);

% Iterate over all vector pairs, excluding self-pairs
hwait=waitbar(0,'Calculating Orientation Correlation Array...');

N=size(ptList,1);
disp(N)

% We need to do this stochastically sometimes
N_rand = min(N,10000); % Anything more than 10,000 points will take forever 
rand_inds = randsample(N, N_rand, false);
ptList_rand = ptList(rand_inds,:);
AngList_rand = AngList(rand_inds);

for i = 1:N_rand
    disp(N_rand-i)
    waitbar(i/N_rand, hwait)
    
    for j = i+1:N_rand
        % Compute distance vector
        dd = ptList_rand(i,:)-ptList_rand(j,:);
        
        % Convert to bin number
        bb = round(dd./db) + numBins + 1;
        bi = bb(1); bj = bb(2);
%         if bi<1; bi=1; end
%         if bj<1; bj=1; end
%         if bi>2*numBins(1)+1; bi=2*numBins(1)+1; end
%         if bj>2*numBins(2)+1; bj=2*numBins(2)+1; end
        
        % Increment histogram count
        hh(bi,bj,1)=hh(bi,bj,1)+1;
        
        % Compute rolling average of angle cosine
        hh(bi,bj,2)= hh(bi,bj,1)/(hh(bi,bj,1)+1)*hh(bi,bj,2)...
                  + 1/(hh(bi,bj,1)+1)...
                    * cosd((AngList_rand(i)-AngList_rand(j)))^2;        
    end
end

% Rotate the matrix by 180 and average with itself to mirror the
% orientations
hh(:,:,1) = hh(:,:,1) + rot90(hh(:,:,1),2);
hh(:,:,2) = (hh(:,:,2) + rot90(hh(:,:,2),2)) ./ 2;

% Enforce Correlation of 1 at center point and give it a count
hh(numBins,numBins,1) = 1;
hh(numBins,numBins,2) = 1;

[X, Y] = meshgrid((-nmSize(1):db(1):nmSize(1)),(nmSize(2):-db(2):-nmSize(2)));
hh(:,:,3) = X;
hh(:,:,4) = Y;

corrLen_x = fminsearch(@(b) sum((exp(-hh(numBins(1)+1,numBins(2)+1:end,3)./b) - ... Exponential decay
                                 hh(numBins(1)+1,numBins(2)+1:end,2)).^2), ... Sum squared error
                       ims.nmWid/5); % Guess 1/5 image width as decay length
                   
corrLen_y = fminsearch(@(b) sum((exp(-hh(1:numBins(1)+1,numBins(2)+1,4)./b) - ... Exponential decay
                                 hh(1:numBins(1)+1,numBins(2)+1,2)).^2), ... Sum squared error
                       ims.nmWid/5); % Guess 1/5 image width as decay length

ims.OrCorr2D = struct('array', hh, ...
                      'corrLen_x', corrLen_x, ...
                      'corrLen_y', corrLen_y);

close(hwait)

end