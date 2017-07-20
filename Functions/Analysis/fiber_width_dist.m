function FWD = fiber_width_dist(ims,settings)

% ORIGIN IS (0,0) AT TOP LEFT CORNER for both pixel and nm space
% X INCREASES DOWN
% Y INCREASES RIGHT

numSamps = min(settings.fibWidSamps,sum(sum(ims.skelTrim)));
angMap = ims.AngMap;
skel = ims.skelTrim;
nmPix = ims.nmPix;
thresh = ims.CEDclean;  % Thresholded black and white image
% imtool(thresh)

[m,n] = size(angMap);

skelList = find(skel);  % Linear indices of pixels in skeletons
skelSamp = datasample(skelList,numSamps,'Replace',false);
skelSamp = skelSamp(~isnan(angMap(skelSamp)));  % Remove samples where the skeleton was on a weird artifact
numSamps = length(skelSamp);

FWD = zeros(size(skelSamp));

for i = 1:numSamps
    
    
    % Set up search vectors and origin
    
    s = skelSamp(i);
    disp(s)
    ang_i = angMap(s);  % Fiber orientation at this skeletal pixel, 0 = horizontal
    wvec = [cosd(ang_i), sind(ang_i)];    % unit vector pointing along width path... the x->y and y->-x axis flip makes this line work.......
    orig_pix = ind2subv([m,n],s) - [0.5, 0.5];  % center of skeletal pixel in pixel space
    
    % parameterize the width path and list the times where it hits
    % gridlines
    
    numPts = 100;
    n2 = 2*(1:numPts)'-1;
    T = [n2*sign(wvec(1))/(2*wvec(1)), n2*sign(wvec(2))/(2*wvec(2))];   % Somehow I derived this...
    Tx = [T(:,1), (1:numPts)', ones(numPts,1)];
    Ty = [T(:,2), (1:numPts)', ones(numPts,1)+1];
    Tlist = [Tx;Ty];
    Tsort = sortrows(Tlist,1);
    
    % go through list of gridline hits and if the pixel beyond that
    % gridline is black, stop and record width
    
    for t = 1:length(Tsort)
        gridpt = orig_pix + Tsort(t,1) * wvec;
        xory = Tsort(t,3);      % =1 if x grid, =2 if y grid
        if xory == 1
            next_pt = [round(gridpt(1))+(wvec(1)>0), ceil(gridpt(2))]; % If moving down, check next pixel down, if up, check up
        else
            next_pt = [ceil(gridpt(1)), round(gridpt(2))+(wvec(2)>0)]; % If moving right, check next pixel to right, if left, check left
        end
        
%         disp(next_pt)
        if next_pt(1)>m || next_pt(1)<1 || next_pt(2)>n || next_pt(2)<1
            check = 1;
        elseif thresh(next_pt(1),next_pt(2)) == 0
            check = 1;
        else
            check = 0;
        end
            
        if check
            pos_bound = gridpt;
            w_pos = norm(pos_bound-orig_pix,2);
            break
        end
    end
    
    % Now change direction of search and go again
    wvec = -wvec;
    
    % parameterize the width path and list the times where it hits
    % gridlines
    
    numPts = 100;
    n2 = 2*(1:numPts)'-1;
    T = [n2*sign(wvec(1))/(2*wvec(1)), n2*sign(wvec(2))/(2*wvec(2))];
    Tx = [T(:,1), (1:numPts)', ones(numPts,1)];
    Ty = [T(:,2), (1:numPts)', ones(numPts,1)+1];
    Tlist = [Tx;Ty];
    Tsort = sortrows(Tlist,1);
    
    % go through list of gridline hits and if the pixel beyond that
    % gridline is black, stop and record width
    
    for t = 1:length(Tsort)
        gridpt = orig_pix + Tsort(t,1) * wvec;
        xory = Tsort(t,3);      % =1 if x grid, =2 if y grid
        if xory == 1
            next_pt = [round(gridpt(1))+(wvec(1)>0), ceil(gridpt(2))]; % If moving down, check next pixel down, if up, check up
        else
            next_pt = [ceil(gridpt(1)), round(gridpt(2))+(wvec(2)>0)]; % If moving right, check next pixel to right, if left, check left
        end
        
        if next_pt(1)>m || next_pt(1)<1 || next_pt(2)>n || next_pt(2)<1
            check = 1;
        elseif thresh(next_pt(1),next_pt(2)) == 0
            check = 1;
        else
            check = 0;
        end
            
        if check
            neg_bound = gridpt;
            w_neg = norm(neg_bound-orig_pix,2);
            break
        end
    end
    
    FWD(i) = (w_pos+w_neg)*nmPix;
end

% save('widtest')

end
    
    