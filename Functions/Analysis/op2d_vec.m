function ims = op2d_vec(ims)

% A super trimmed-down version of FiberApp's OrderParameter2D that uses
% default settings and the full image with no random images.

% Adds fibLengthDensity as a field to ims (in 1/um) and
% op2d, a structure with fields:
%   Sfull: the full image fitted 2d order parameter
%   decayLen: the fitted decay length parameter, aka lambda_C
%   xdata: the vector of frame sizes used to calculate S2D
%   S_im: the vector of S2D values at all frame sizes in xdata

step = ims.settings.fiberStep_nm;
xy = {ims.Fibers(:).xy_nm};  % This is how FiberApp handles the xy data
gridStep = 4*step;

% Check processing area values
p_all = cellfun(@get_points, xy, 'UniformOutput', false);
p_all = [p_all{:}]; % combine points of all fibers
xMin = floor(min(p_all(1,:)));
xMax = ceil(max(p_all(1,:)));
yMin = floor(min(p_all(2,:)));
yMax = ceil(max(p_all(2,:)));

% Image size in grid steps (cells)
numCellX = floor((xMax-xMin)/gridStep);
numCellY = floor((yMax-yMin)/gridStep);
procStepNum = min([numCellX, numCellY]);
procLength = procStepNum*gridStep;

% Pick data within the range
% p_im - coordinates of vectors' start points
% v_im - vector coordinates
xy_im = cellfun(@(xy)get_points_in_range(xy, xMin, xMax, yMin, yMax), xy, ...
    'UniformOutput', false);
xy_im = [xy_im{:}];

p_im = cellfun(@get_points, xy_im, 'UniformOutput', false);
v_im = cellfun(@get_vectors, xy_im, 'UniformOutput', false);
p_im = [p_im{:}];
v_im = [v_im{:}];

% Calculate fiber length density (in 1/um)
fibLengthDensity = 1000*size(p_im, 2)*step/((xMax-xMin)*(yMax-yMin));
ims.fibLengthDensity = fibLengthDensity;

% Calculate Order Parameter 2D vs Box size
[xdata, S_im] = OP2D(p_im, v_im, xMin, yMin, ...
    gridStep, procStepNum, numCellX, numCellY, 0);

% Fit data
BETA = lsqnonlin(@(B) B(1)+(1-B(1)).*exp(-xdata./B(2)) - S_im,...   % fitting function (should be all 0's)
                 [0.5,1000],... % initial guess
                 [0 1E-3],...   % lower bound
                 [1 Inf]);      % upper bound

ims.op2d = struct('Sfull',BETA(1),...
                  'decayLen',BETA(2),...
                  'xdata',xdata,...
                  'S_im',S_im);

end


function p = get_points_in_range(xy, xMin, xMax, yMin, yMax)
ind = diff([0, xMin <= xy(1,:) & xy(1,:) <= xMax & yMin <= xy(2,:) & xy(2,:) <= yMax, 0]);
in = find(ind == 1);
out = find(ind == -1);
p = cell(1, length(in));
for l = 1:length(in)
    p{l} = xy(:, in(l):out(l)-1);
end

end

function p = get_points(xy)
p = xy(:,1:end-1);

end

function v = get_vectors(xy)
v = diff(xy, 1, 2);
l = sqrt(sum(v.^2));
v = v./[l; l];

end

function [p_rand, v_rand] = move_and_rotate(xy, xc, yc, angle)
xy = xy{1};

% Center at (0; 0)
xy = bsxfun(@minus, xy, xy(:,1));

% Move to (xc; yc) and turn by angle
p_rand = [xc + xy(1,:).*cos(angle) - xy(2,:).*sin(angle);
    yc + xy(1,:).*sin(angle) + xy(2,:).*cos(angle)];

% Get vectors and points
v_rand = diff(p_rand, 1, 2);
l = sqrt(sum(v_rand.^2));
v_rand = v_rand./[l; l];

p_rand = p_rand(:,1:end-1);

end

function [xdata, ydata] = OP2D(p_all, v_all, xMin, yMin, ...
    gridStep, procStepNum, numCellX, numCellY, isCircleArea)
%ORDERPARAMETER2D Calculate Order Parameter 2D vs Box size

% Arrange and calculate data in cells
SIZE = [numCellY numCellX];
A = zeros(SIZE);
B = zeros(SIZE);
N = zeros(SIZE, 'uint32');

p_all_x = p_all(1,:);
p_all_y = p_all(2,:);
% First cycle for horizontal slicing
for k = 1:numCellY
    sY = yMin + (k-1)*gridStep;
    ind = find(sY <= p_all_y & p_all_y < sY + gridStep);
    if isempty(ind); continue; end
    
    p_slice_x = p_all_x(ind);
    v_slice = v_all(:,ind);
    
    % Second cycle for vertical slicing
    for l = 1:numCellX
        sX = xMin + (l-1)*gridStep;
        ind = find(sX <= p_slice_x & p_slice_x < sX + gridStep);
        if isempty(ind); continue; end
        
        v_block = v_slice(:,ind);
        A(k,l) = sum(v_block(1,:).^2);
        B(k,l) = sum(prod(v_block));
        N(k,l) = size(v_block,2);
    end
end

xdata = gridStep*(1:procStepNum);
ydata = zeros(1, procStepNum);

% First step (k = 1)
Ah = A;
Av = zeros(SIZE);
a = A;

Bh = B;
Bv = zeros(SIZE);
b = B;

Nh = N;
Nv = zeros(SIZE, 'uint32');
n = N;

% In case of circle area
if isCircleArea
    % These values contain information about corners of the square box
    % 1 - top left corner
    % 2 - top right corner
    % 3 - bottom left corner
    % 4 - bottom right corner
    
    cir1a = zeros(SIZE);
    cir2a = zeros(SIZE);
    cir3a = zeros(SIZE);
    cir4a = zeros(SIZE);
    
    cir1b = zeros(SIZE);
    cir2b = zeros(SIZE);
    cir3b = zeros(SIZE);
    cir4b = zeros(SIZE);
    
    cir1n = zeros(SIZE, 'uint32');
    cir2n = zeros(SIZE, 'uint32');
    cir3n = zeros(SIZE, 'uint32');
    cir4n = zeros(SIZE, 'uint32');
    
    % Construct a matrix for a quarter of the box. Each element of this
    % matrix contain a circle diameter, at which this element doesn't
    % belong to the circle anymore
    [xx, yy] = meshgrid(1:floor(procStepNum/2));
    notInCircle = uint16(ceil(2*(xx+yy-1)+sqrt(4*(2*xx.*yy-xx-yy)+2)));
end

val = sqrt((2*a-double(n)).^2+(2*b).^2)./double(n);
% ----------
% val(n==0) = 0;
% ydata(1) = mean(val(:));
% ----------
ydata(1) = mean(val(n~=0));

% Remaining steps (k = 2:procStepNum)
for k = 2:procStepNum
    Ah = Ah(2:end,1:end-1) + A(k:end,k:end);
    Av = Av(1:end-1,2:end) + A(k-1:end-1,k:end);
    a = a(1:end-1,1:end-1) + Av + Ah;
    
    Bh = Bh(2:end,1:end-1) + B(k:end,k:end);
    Bv = Bv(1:end-1,2:end) + B(k-1:end-1,k:end);
    b = b(1:end-1,1:end-1) + Bv + Bh;
    
    Nh = Nh(2:end,1:end-1) + N(k:end,k:end);
    Nv = Nv(1:end-1,2:end) + N(k-1:end-1,k:end);
    n = n(1:end-1,1:end-1) + Nv + Nh;
    
    % In case of circle area
    if isCircleArea
        cir1a = cir1a(1:end-1, 1:end-1);
        cir2a = cir2a(1:end-1, 2:end);
        cir3a = cir3a(2:end, 1:end-1);
        cir4a = cir4a(2:end, 2:end);
        
        cir1b = cir1b(1:end-1, 1:end-1);
        cir2b = cir2b(1:end-1, 2:end);
        cir3b = cir3b(2:end, 1:end-1);
        cir4b = cir4b(2:end, 2:end);
        
        cir1n = cir1n(1:end-1, 1:end-1);
        cir2n = cir2n(1:end-1, 2:end);
        cir3n = cir3n(2:end, 1:end-1);
        cir4n = cir4n(2:end, 2:end);
        
        [i, j] = find(notInCircle == k);
        for l = 1:length(i)
            ind1i = i(l):i(l)+SIZE(1)-k;
            ind1j = j(l):j(l)+SIZE(2)-k;
            ind2i = 1-i(l)+k:1-i(l)+SIZE(1);
            ind2j = 1-j(l)+k:1-j(l)+SIZE(2);
            
            cir1a = cir1a + A(ind1i, ind1j);
            cir2a = cir2a + A(ind1i, ind2j);
            cir3a = cir3a + A(ind2i, ind1j);
            cir4a = cir4a + A(ind2i, ind2j);
            
            cir1b = cir1b + B(ind1i, ind1j);
            cir2b = cir2b + B(ind1i, ind2j);
            cir3b = cir3b + B(ind2i, ind1j);
            cir4b = cir4b + B(ind2i, ind2j);
            
            cir1n = cir1n + N(ind1i, ind1j);
            cir2n = cir2n + N(ind1i, ind2j);
            cir3n = cir3n + N(ind2i, ind1j);
            cir4n = cir4n + N(ind2i, ind2j);
        end
        
        aa = a - cir1a - cir2a - cir3a - cir4a;
        bb = b - cir1b - cir2b - cir3b - cir4b;
        nn = double(n - cir1n - cir2n - cir3n - cir4n);
    else
        aa = a;
        bb = b;
        nn = double(n);
    end
    
    val = sqrt((2*aa-nn).^2+(2*bb).^2)./nn;
    % ----------
%     val(nn==0) = 0;
%     ydata(k) = mean(val(:));
    % ----------
    ydata(k) = mean(val(nn~=0));
end

xdata = [0, xdata];
ydata = [1, ydata];

end

