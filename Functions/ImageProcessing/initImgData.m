function ims = initImgData(imfile)

[img, c_map] = imread(imfile);
if size(img,3)>3
    img = img(:,:,1:3);
end

if ispc
    sep = '\';
else
    sep = '/';
end

LS = findLastSlash(imfile);
LD = findLastDot(imfile);
imName = imfile(LS+1:LD-1); % just the name without extension
imPath = imfile;    % the whole path, with extension
imNamePath = imfile(1:LD-1);    % whole path without extension

ims = struct();
ims.img = img;
ims.c_map = c_map;
ims.imName = imName;
ims.imPath = imPath;
ims.imNamePath = imNamePath;
ims.figSavePath = [imNamePath, sep, imName];  % .../<directory with image>/<image name>/<image name>

% generate grayscale
if ndims(img)>2
    ims.gray = rgb2gray(img);
else
    ims.gray = img;
end

% ims = YBSimpleSeg(ims);

end

function ims = YBSimpleSeg(ims)

%YBSeg Yanowitz-Bruckstein image segmentation with fiber unentanglement
%   SP is the structure path
%   File is the file path from current active dir
%   Dim is the image dimension in nm

ims.edge = edge(ims.gray,'canny');         % Apply a Canny edge finder - 1's at the edges, 0's elsewhere
ims.grayDouble = double(ims.gray);                     % Turn the grey image into double prec.
ims.edgeDouble = double(ims.edge);                     % Turn the edge image into double prec.
ims.initThresh = ims.grayDouble.*ims.edgeDouble;                        % Fill in the grey values of the edge pixels in a new image file                      

ims.threshSurf = YBiter(ims.initThresh);                    % Perform Yanowitz-Bruckstein surface interpolation to create threshold surface from edge gray values

ims.yb_bw = ims.gray>ims.threshSurf;                          % Segment the image; pixels above threshold surface are white, if not, black
ims.std_bw = im2bw(ims.img);
ims.yb_gray = mat2gray(ims.grayDouble-ims.threshSurf);

end

function Vf = YBiter(V0)

%YBiter Yanowitz/Bruckstein surface interpolation

w = 1;
[m,n] = size(V0);
InitVal = mean(V0(V0~=0)); % Start non-edges as the mean of the edges because 0's aint' working
Vupdate = V0==0;           % A logical array of pixels to update on each iteration
Vupdate(1,:) = 0;
Vupdate(:,1) = 0;
Vupdate(m,:) = 0;
Vupdate(:,n) = 0;

V0(V0==0)=InitVal;         % Put the average edge values in the non-edge cells

Vnew = zeros(m,n);         % Initialize the updated threshold surface
Vold = V0;

iter = 0;
maxiter = 40;

% tic
while iter<maxiter
    iter = iter+1;
%     disp(iter)
    
    Lap = del2(Vold);
    Vnew = Vold + Vupdate .* (w .* Lap);
    Vold = Vnew;
    
end
% toc
Vf = Vnew;

end