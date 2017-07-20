function ims = BreakHighCurvSegs(ims)

% Start with Initial segs
segsBroken = ims.segsInit;

for i = 1:length(ims.fibSegs)
    
    % Get subscript indices of sorted pixel path
    sps = ims.fibSegs(i).sortPixSubs;
    
    % Compute the length between each pixel (basically a list of 1's and
    % 1.4's), because each pixel is not necessarily an equal unit of fiber
    % length
    sortPixDist = sqrt(sum((sps(2:end,:)-sps(1:end-1,:)).^2,2));
    sortPixPathLen = [0; cumsum(sortPixDist)];
    
    % Convert to fractional path length
    sortPixPathLen = sortPixPathLen/max(sortPixPathLen);
    
    % Convert to path bins that correspond to fiber vertices
    pathBins = round(sortPixPathLen*(length(ims.fibSegs(i).xy))-1)+1;
    
    % Remove pixels from segsBroken if they belong to high curvature path
    % segments
    high_curv_bins = find(ims.fibSegs(i).curv > ims.settings.maxCurv);
    high_curv_pix = ~~MultiEquiv(pathBins,high_curv_bins);
    
    segsBroken(ims.fibSegs(i).sortPixInds(high_curv_pix))=0;

end

ims.segsInit = segsBroken;

% And then we restart the whole process!

end