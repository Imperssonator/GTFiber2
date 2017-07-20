function [S, MO, cent, n, Frames, Sfull, Smod, BETA] = getS2Dauto(img,settings)

addpath('Functions')
addpath('Functions/coherencefilter_version5b')

% Make a dummy figure so the filter will be happy
auto_handles = struct();

% Initialize image data structure and pixelize (convert from nm to pix) setting values
ims = initImgData(img); % initialize the image structure, which generates a grayscale and some simple thresholded stuff to start with
[settings, ims] = pix_settings(settings,ims);
if settings.figSave
    mkdir(ims.imNamePath);  % make the directory to save the results figures
end

% Run the filter bank at the current settings
auto_handles.ims = ims;
auto_handles.settings = settings;
auto_handles = main_filter(auto_handles);
ims = auto_handles.ims;

% Display the Angle Color Map for pixels that are in the skeleton
% disp('Generating Figures...')
[S, MO, cent, n] = AngleColorMapInt(ims,settings);
if settings.fullOP
%     [Frames, Sfull, Smod] = op2d_am(ims, 0, maxFrame, settings.frameStep, settings.gridStep, 0);
    [Frames, Sfull, Smod, BETA] = op2d_am(ims, settings);
end

end

function [S, MO, cent, n] = AngleColorMapInt(ims,settings) %AngMap,BW,varargin)

% Take AngMap (orientations at every pixel) and BW (1's for fibers, 0's
% otherwise), and produce an image with black pixels for non-fibers, and
% colors for fibers corresponding to orientation on a double-wrapped color
% wheel

sat = 0.95;

AngMap = ims.AngMap;
BW = ims.skelTrim;
AngMapNaN = AngMap;
AngMapNaN(~BW) = NaN;
AngMapNaN(AngMapNaN<0) = AngMapNaN(AngMapNaN<0)+180;
[m, n] = size(AngMap);
AngIm = zeros(m,n,3);
for i = 1:m
    for j = 1:n
        if isnan(AngMapNaN(i,j))
            AngIm(i,j,:) = [0 0 0];
        else
            AngIm(i,j,:) = hsv2rgb([2*AngMapNaN(i,j)/360, 1, sat]);
        end
    end
end

AngImDiskFilt = imfilter(AngIm,fspecial('disk',2),'replicate');
AngImHSV = rgb2hsv(AngImDiskFilt);
AngImHSV(:,:,3) = ones(size(AngImHSV,1),size(AngImHSV,2)).*AngImHSV(:,:,3)>0.1;

AngIm = imgaussian(hsv2rgb(AngImHSV),0.7);

if settings.figSwitch
    figure; imshow(AngIm)
end
if settings.figSave
    imwrite(AngIm,[ims.figSavePath, '_AM', '.tif']);
    imwrite(ims.CEDclean,[ims.figSavePath, '_BW', '.tif']);
end

[n, cent, S, MO] = ODistInt(AngMapNaN,5,1,ims,settings);

end

function [n, centersd, S, MeanOrient] = ODistInt(AngMap,angleStep,Label,ims,settings)

% Portions of this function are subject to the following license agreement:
% Copyright (c)
% 2011-2014, ETH Zurich (Switzerland)
% 2015-2016, Ivan Usov
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 
% 2. Publication of any type of results based on the use of this open source code legally requires citation to the original publication: Usov, I and Mezzenga, R. FiberApp: an Open-source Software for Tracking and Analyzing Polymers, Filaments, Biomacromolecules, and Fibrous Objects. Macromolecules, 2015, 49, 1269-1280.
% 
% 3. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


% Hard Coded Stuff
% ___________________
dirScale = 0.67; % what fraction of max bin should director segment length be
axGray = 0.84; % how light should the radial/angular axes be
outerGray = 0.4;  % brightness of outer radial limit of polar plot
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
centersd = 180*centers/pi; % recalculate into deg in case of saving to a text file

n = [n(end), n, n];
maxn = max(n);
dirLen = dirScale*maxn;

if ~settings.figSwitch && ~settings.figSave
    return
end

% Plot distribution in a new figure
ofig = figure('NumberTitle', 'off', 'Name', ['ODist ' datestr(now, 'HH:MM:SS dd/mm/yy')]);
ax1 = gca;
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

if settings.figSave
    hgexport(ofig, [ims.figSavePath, '_OD', '.tif'],  ...
        hgexport('factorystyle'), 'Format', 'tiff');
    close(ofig)
end

end

function v = get_vectors(xy)
v = diff(xy, 1, 2);
l = sqrt(sum(v.^2));
v = v./[l; l];

end

function [AngMap, Dirs] = RecoverAnglesOD(J)

% Take an [m x n x 4] structure tensor and return an m x n matrix of angles
% "AngMap" and an [m x n x 2] array of vectors "Dirs" corresponding to those
% angles

% disp('Finding Angles...')
% tic
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
% toc

end



