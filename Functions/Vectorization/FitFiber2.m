function [xy,curv] = FitFiber2(ims,fibSegs)

% Basically taking the active contours algorithm from FiberApp and giving it
% an extremely good initial guess - just a list of the pixels that are the
% backbone of a fiber segment and an image that is the binary ground truth
% of the fiber segment

if length(fibSegs)>1
    sortPixInds = sort_segs(fibSegs);
else
    sortPixInds = fibSegs(1).sortPixInds;
end

sortPixSubs = ind2subv(size(ims.SegLabels),sortPixInds);

fibImg = int16(ims.gray);
[gradX, gradY] = gradient2Dx2(fibImg);

fiberStep = ims.settings.fiberStep;  % Number of pixels a discrete step should take

% Short names for fiber tracking parameters and data
xy_col = sortPixSubs;    % column vector of i,j coords of pixels in order
xy = flipud(xy_col');
xy = distributePoints(xy,fiberStep);
a = 0;
b = 10;
g = 10;
k1 = 20;
k2 = 0;

% meanInt = mean(ims.gray(:));
bgInt = mean(ims.gray(~ims.CEDclean));
fiberIntensity = mean(ims.gray(~~ims.CEDclean));

% Apply fitting algorithm FA.iterations times
for k = 1:50
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
    vf = k1*[vfx; vfy]/(fiberIntensity-bgInt);
    
    % Normalized vectors of fiber ends
    v_start = xy(:,1) - xy(:,2);
    v_start = v_start/sqrt(sum(v_start.^2));
    v_end = xy(:,end) - xy(:,end-1);
    v_end = v_end/sqrt(sum(v_end.^2));
    
    % Image intensities at the ends
    int_start = interp2D(fibImg, xy(1,1), xy(2,1), 'cubic');
    int_end = interp2D(fibImg, xy(1,end), xy(2,end), 'cubic');
    
    % Full external energy
    vf(:,1) = vf(:,1) + k2*v_start*(int_start)/(fiberIntensity-bgInt);
    vf(:,end) = vf(:,end) + k2*v_end*(int_end)/(fiberIntensity-bgInt);
    
    % Next position
    xy = distributePoints((g*xy+vf)/M, fiberStep);
end

if size(xy,2)>=3
    dxy=(xy(:,3:end)-xy(:,1:end-2))./(2*ims.settings.fiberStep);
    ddxy=(xy(:,3:end)-2.*xy(:,2:end-1)+xy(:,1:end-2))./(ims.settings.fiberStep)^2;
    curv=(dxy(1,:).*ddxy(2,:)-dxy(2,:).*ddxy(1,:))./sum(dxy.^2,1).^(1.5);
    curv=[0,abs(curv),0];
else
    curv = zeros(1,size(xy,2));
end


end

function sortPixInds = sort_segs(fibSegs)

xy_ends = [];
seg_sort = [];
for i = 1:length(fibSegs);
    xy_ends = [xy_ends;
        fibSegs(i).sortPixSubs(1,:);
        fibSegs(i).sortPixSubs(end,:)];
    seg_sort = [seg_sort;
        [i 1];
        [i 2]];
end
xy_ends=flipud(xy_ends');

xy_arc = zeros(size(xy_ends,2),1);
xy_seek = fibSegs(1).xy_seek;
for i = 1:length(xy_arc);
    xy_dist = sum((xy_seek-repmat(xy_ends(:,i),1,size(xy_seek,2))).^2,1);
    [~, xy_arc(i)] = min(xy_dist);
end

[ends_sort,orig_ends] = sort(xy_arc);

seg_sort = seg_sort(orig_ends,:);

sortPixInds = [];
for i = 1:2:size(seg_sort,1);
    if seg_sort(i,2)==1;
        sortPixInds = [sortPixInds;
            fibSegs(seg_sort(i,1)).sortPixInds];
    else
        sortPixInds = [sortPixInds;
            flipud(fibSegs(seg_sort(i,1)).sortPixInds)];
    end
end

end

