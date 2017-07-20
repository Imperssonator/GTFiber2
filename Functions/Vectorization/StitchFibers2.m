function ims = StitchFibers2(ims)

ims = CleanSkeleton(ims);

% 2. Label segments and Make EndLib
hwait = waitbar(0,'Building Endpoint Library...');
ims = BuildEndLib(ims);

% 3. Vectorize Segments - once we are here, we do not go back
waitbar(0.2,hwait,'Vectorizing Segments...');
ims = fitAllSegments(ims);

% 4. Break high curvatures
waitbar(0.3,hwait,'Breaking High Curvatures...');
ims = BreakHighCurvSegs(ims);

% 4.1 Rebuild EndLib
waitbar(0.4,hwait,'Rebuilding Endpoint Library...');
ims = BuildEndLib(ims);

% 4.2 Re-fit new segments
waitbar(0.5,hwait,'Re-vectorizing Segments...');
ims = fitAllSegments(ims);

% % 5. Remove Non-conforming widths (measure widths)
% waitbar(0.5,hwait,'Removing Out-of-range Widths...');
% ims = SegWidths(ims);
% ims = RemoveBadWidths(ims);

% % 5. Rebuild EndLib
% waitbar(0.6,hwait,'Rebuilding Endpoint Library...');
% ims = BuildEndLib(ims);
% ims = fitAllSegments(ims);

% 6. Match Search, Ranking/Pairing, Percolation
waitbar(0.7,hwait,'Segment Matching and Percolation...');
ims = SegMatch(ims);

% 7. Vectorize Fibers
waitbar(0.8,hwait,'Vectorizing Final Fibers...');
ims = fitAllFibers(ims);
ims = RemoveShortFibers(ims);

% 8. Calculate op2d and fiber length and width distributions
waitbar(0.9,hwait,'Analyzing Structure...');
ims = op2d_vec(ims);
ims = calc_orient_dist(ims);
ims = FiberLengths(ims,0);
ims = FiberWidths(ims);

% save('last_result','ims')
close(hwait)

end