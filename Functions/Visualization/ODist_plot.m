function hf = ODist_plot(ims,figSave)

% Produce the radar-style polar plot of the orientation distribution from
% and ims structure that has had its ODist calculated with calc_orient_dist

if exist('figSave','var')~=1
    figSave=0;
end

% Hard Coded Stuff
% ___________________

dirScale = 0.5; % what fraction of max bin should director segment length be
axGray = 0.84; % how light should the radial/angular axes be
outerGray = 0.4;
lab180scale = 1.2;
% ___________________

% Pull distribution and director orientation out of ims
n = ims.ODist.n;
centers_polar = ims.ODist.centers_polar;
MeanOrient = ims.ODist.director;
maxn = max(n);
dirLen = dirScale*maxn;

% Make polar plot with director
hf=figure;
ax1 = gca;
polar(ax1, centers_polar, n, '-b');
hold(ax1,'on')
polar(ax1,[MeanOrient*pi/180, 0, MeanOrient*pi/180+pi],[dirLen, 0, dirLen],'-k');

% Fix the axis labels
pos180 = ax1.XTick(1) * lab180scale;

htext = findall(hf,'Type','text');
keep_inds = [find(cellfun(@(x) strcmp('0',x),{htext(:).String}),1);...
             find(cellfun(@(x) strcmp('90',x),{htext(:).String}),1);...
             find(cellfun(@(x) strcmp('180',x),{htext(:).String}),1);...
             find(cellfun(@(x) strcmp('270',x),{htext(:).String}),1)...
             ];

for t = 1:length(htext)
    if not(ismember(t,keep_inds))
        htext(t).String = '';
    end
    htext(t).FontSize = 24;
    if t==keep_inds(3) % if we're dealing with the 180 label
        htext(t).Position(1) = pos180;
    end
end

% Fix line widths
hlines = findall(hf,'Type','line');
for i = 1:length(hlines)
set(hlines(i),'LineWidth',2);
end

% Draw outline darker
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
% set(hf,'color','none');

% Reposition and scale figure
hf.Position = [600 300 400 400];

if figSave
    hgexport(hf, [ims.figSavePath, '_OD', '.tif'],  ...
        hgexport('factorystyle'), 'Format', 'tiff');
    close(hf)
end

% text('Units', 'normalized', 'Position', [-0.09 0.16], ...
%     'BackgroundColor', [1 1 1], ...
%     'String', ['\itS\rm_{full} = ' num2str(ims.op2d.Sfull, 2)]);

end

