function ims = fitAllFibers(ims)

numFibs = max(max(ims.FiberLabels));

% The first goal here is to take individual segments and turn them into a
% list of pixels where the list increases by distance from one end of the
% segment, which is quantified by the bwdistgeodesic function



for f = 1:numFibs
%     disp(i)
    ims = FitFiber(ims,f);
    ims.Fibers(f).xy_nm = ims.Fibers(f).xy .* ims.nmPix;
%     high_curve = [high_curve; sorted_seg(ims.fibSegs(i).curv>1e-03,1)];
end

% high_curve_im = zeros(size(ims.SegLabels));
% high_curve_im(high_curve) = 1;
% imtool(high_curve_im)

end

function ims = FitFiber(ims,fibNum)

% Basically taking the active contours algorithm from FiberApp and giving it
% an extremely good initial guess - just a list of the pixels that are the
% backbone of a fiber segment and an image that is the binary ground truth
% of the fiber segment

FiberSegs = ims.Fibers(fibNum).FiberSegs;                   % Labels of segments in fiber
numFibSegs = length(FiberSegs);
StartInds = ims.Fibers(fibNum).Fiber(1:2:end);              % Indices of the first endpoint of each segment in the chain
sortPixInds = [];

for i = 1:numFibSegs
    
    if StartInds(i) <= length(ims.EndLib)                           % If the start ind is the 'first' index of that segment
        sortPixInds = [sortPixInds;                                 % Take the sorted pixel indices as they are
                       ims.fibSegs(FiberSegs(i)).sortPixInds];
    else                                                            % Otherwise, flip them upside down
        sortPixInds = [sortPixInds;
                       ims.fibSegs(FiberSegs(i)).sortPixInds(end:-1:1)];
    end
    
end

ims.Fibers(fibNum).sortPixInds = sortPixInds;
ims.Fibers(fibNum).sortPixSubs = ind2subv(size(ims.SegLabels),sortPixInds);

fibImg = int16(ims.gray);
[gradX, gradY] = gradient2Dx2(fibImg);

fiberStep = ims.settings.fiberStep;  % Number of pixels a discrete step should take
    
% Short names for fiber tracking parameters and data
xy_col = ims.Fibers(fibNum).sortPixSubs;    % column vector of i,j coords of pixels in order
xy = flipud(xy_col');
xy = distributePoints(xy,fiberStep);
a = 0;
b = 5;
g = 10;
k1 = 20;
k2 = 0;
fiberIntensity = 255;

% Apply fitting algorithm FA.iterations times
for k = 1:20
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

% Save fiber data
ims.Fibers(fibNum).length =...
    (size(xy,2)-1)* ims.settings.fiberStep_nm;
ims.Fibers(fibNum).curv = curv;
ims.Fibers(fibNum).xy = xy;
ims.Fibers(fibNum).xy_nm = xy_nm;

end

