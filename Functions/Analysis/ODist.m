function [n, centers, S, MeanOrient] = ODist(AngMap,angleStep,Label,varargin)

if length(varargin)>0
    noFig = varargin{1};
else
    noFig = 0;
end

% Hard Coded Stuff
% ___________________

dirScale = 0.67; % what fraction of max bin should director segment length be
axGray = 0.84; % how light should the radial/angular axes be
outerGray = 0.4;

% ___________________

% Orientations to Unit Vectors
AMx = cosd(AngMap); amxlin = AMx(:)';
AMy = -sind(AngMap); amylin = AMy(:)';
amxclean = amxlin(~isnan(amxlin));
amyclean = amylin(~isnan(amylin));
vect = [amxclean; amyclean];

% Calculate 2D order parameter (S)
A = sum(vect(1,:).^2);      
B = sum(prod(vect));
N = size(vect, 2);
S = sqrt((2*A-N)^2+4*B^2)/N;

% Calculate orientation distribution
vect(:, vect(2,:)>0) = - vect(:, vect(2,:)>0); % Turn the coordinate system from informatics into geometrical 
orientation = acos(vect(1,:));
% disp(vect)

% Calculate Director
JJJ = zeros(length(orientation),1,4);
JJJ(:,1,1) = cos(orientation).^2;
JJJ(:,1,2) = cos(orientation).*sin(orientation);
JJJ(:,1,3) = cos(orientation).*sin(orientation);
JJJ(:,1,4) = sin(orientation).^2;
MeanJ = mean(JJJ,1);
MeanOrient = RecoverAnglesOD(MeanJ);
% disp(MeanOrient)

% Check angleStep value
stepNum = round(180/angleStep);
angleStep = 180/stepNum;

% Calculate for 0-pi range
n = histc(orientation, linspace(0, pi, stepNum+1));
n(1) = n(1) + n(end);
n(end) = []; % remove orientation = pi values

% if polarCoord % according to the selected coordinate system
% Reflect through the origin to the full 360 deg range

step = pi*angleStep/180; % angle step in rad
centers = - step/2 + linspace(0, 2*pi, 2*stepNum+1) ;
n = [n(end), n, n];
maxn = max(n);
dirLen = dirScale*maxn;

if noFig
    return
end

% Plot distribution in a new figure
ofig = figure('NumberTitle', 'off', 'Name', ['ODist ' datestr(now, 'HH:MM:SS dd/mm/yy')]);
ax1 = gca;

% Create Polar Plot
% First create outer limit
% P = polar(ax1,centers, 1500 * ones(size(centers)));
% set(P, 'Visible', 'off')
% hold on
% Then plot
polar(ax1,centers, n, '-b');
hold on
polar(ax1,[MeanOrient*pi/180, 0, MeanOrient*pi/180+pi],[dirLen, 0, dirLen],'-k');

% Fix the axis labels
htext = findall(ofig,'Type','text');
for t = 13:length(htext)
    htext(t).String = '';
end
for t = 1:12
    htext(t).FontSize = 20;
end

if ispc
    s_font=28;
else
    s_font=40;
end
% Add figure label and S2D label
text('Units', 'normalized', 'Position', [-0.19 0.10], ...
    'BackgroundColor', [1 1 1], ...
    'String', ['S_{2D} = ', char(10), ' ', num2str(S, 2)],...
    'FontSize', s_font);
centers = 180*centers/pi; % recalculate into deg in case of saving to a text file

% Fix line widths
hlines = findall(ofig,'Type','line');
for i = 1:length(hlines)
set(hlines(i),'LineWidth',3);
end

ymax = [];
for i = 1:length(hlines)
    if size(hlines(i).YData,2)>100
        ymax = [ymax; [i, max(hlines(i).YData)]];
    end
end
[highy, outmax] = max(ymax(:,2));
outline = ymax(outmax,1);

% Fix line colors
set(hlines(3:end),'Color',axGray.* [1 1 1]);
set(hlines(outline),'Color',outerGray.*[1 1 1]);

% Reposition and scale figure
ofig.Position = [123 147 581 559];

end

function v = get_vectors(xy)
v = diff(xy, 1, 2);
l = sqrt(sum(v.^2));
v = v./[l; l];

end

function [AngMap, Dirs] = RecoverAnglesOD(J)

% disp('Finding Angles...')

[m, n] = size(J(:,:,1));
JE = zeros(m,n,2);
JV = zeros(m,n,4);
for i = 1:m
    for j = 1:n
        [V,D] = eig(reshape(J(i,j,:),[2 2]));
        JE(i,j,:) = diag(D);
        JV(i,j,:) = V(:);
    end
end
JESort = sort(JE,3,'descend');
JVSort = zeros(size(JV));
for i = 1:m
    for j = 1:n
        if JESort(i,j,1)~=JE(i,j,1)
            JVSort(i,j,:) = [JV(i,j,3) JV(i,j,4) JV(i,j,1) JV(i,j,2)];
        else
            JVSort(i,j,:) = JV(i,j,:);
        end
    end
end

% Coher = ((JESort(:,:,1)-JESort(:,:,2))./(JESort(:,:,1)+JESort(:,:,2))).^2;
Dirs = JVSort(:,:,1:2);
AngMap = atan2d(JVSort(:,:,2),JVSort(:,:,1));      % last evec is 'coherence orientation'


end
