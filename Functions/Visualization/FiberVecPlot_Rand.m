function out = FiberVecPlot_Rand(ims)

step = ims.settings.fiberStep_nm / 2;
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

xc = random('Uniform', xMin, xMax, size(xy_im));
yc = random('Uniform', yMin, yMax, size(xy_im));
angle = random('Uniform', 0, 2*pi, size(xy_im));

[p_rand, v_rand] = arrayfun(@move_and_rotate, xy_im, xc, yc, angle, ...
    'UniformOutput', false);

[h,w] = size(ims.gray);
f1 = figure();
f1.Position = [1 1 1+w 1+h];
f1.Color = [1 1 1];
ha = axes('parent',f1);
hold on

for i = 1:length(p_rand)
XYi = p_rand{i};
plot(ha,XYi(1,:),XYi(2,:),'-b','LineWidth',1)
% plot(ha,XYi(1,:),XYi(2,:),'ob','MarkerSize',8)
end
% axis equal
set(ha,'Ydir','reverse')
ax = ha;
ax.XLim = [0 w*ims.nmPix];
ax.YLim = [0 h*ims.nmPix];
ax.Visible = 'off';
axis equal
ax.Position = [0 0 1 1];


function p = get_points_in_range(xy, xMin, xMax, yMin, yMax)
ind = diff([0, xMin <= xy(1,:) & xy(1,:) <= xMax & yMin <= xy(2,:) & xy(2,:) <= yMax, 0]);
in = find(ind == 1);
out = find(ind == -1);
p = cell(1, length(in));
for l = 1:length(in)
    p{l} = xy(:, in(l):out(l)-1);
end

function p = get_points(xy)
p = xy(:,1:end-1);

function v = get_vectors(xy)
v = diff(xy, 1, 2);
l = sqrt(sum(v.^2));
v = v./[l; l];

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