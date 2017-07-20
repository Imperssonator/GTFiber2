function ims = CleanSkeleton(ims)

% Starting from the initial skeleton, remove branches, groom, and clean it
% until it consists only of unbranched, isolated segments.

% Closing to remove small holes
% skelClose = bwmorph(imclose(ims.skel,strel('disk',1)),'skeleton',Inf);

% Remove branches up to a certain size
ims.skelTrim = RemoveBranches(ims.skel,ims.settings.maxBranchSize);

% Find branch points, dilate them, remove the dilated points to separate
% all segments
branchPts = bwmorph(ims.skelTrim,'branchpoints');
bigBranchPts = imdilate(branchPts,ones(3));
skelNoBranch = ims.skelTrim&~bigBranchPts;

% Remove segments that are shorter than the max branch size
ims.segsInit = bwareaopen(skelNoBranch,ims.settings.maxBranchSize,8);

% Some segments may still have holes, this removes them.

ims.segsInit = RemoveNonLines(ims.segsInit);

end