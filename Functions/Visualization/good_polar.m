function [hf, Fim] = good_polar(angles,n)

% Produce a good-looking radar plot given bins and counts

% Hard Coded Stuff
% ___________________
axGray = 0.84; % how light should the radial/angular axes be
outerGray = 0.4;
lab180scale = 1.5;
% ___________________

maxn = max(n);

% Make polar plot with director
hf=figure;
ax1 = gca;
polar(ax1, angles, n, '-b');
hold(ax1,'on')


% Fix the axis labels
pos180 = ax1.XTick(1) * lab180scale;

htext = findall(hf,'Type','text');
% keep_inds = [find(cellfun(@(x) strcmp('0',x),{htext(:).String}),1);...
%              find(cellfun(@(x) strcmp('90',x),{htext(:).String}),1);...
%              find(cellfun(@(x) strcmp('180',x),{htext(:).String}),1);...
%              find(cellfun(@(x) strcmp('270',x),{htext(:).String}),1)...
%              ];
keep_inds=[];

for t = 1:length(htext)
    if not(ismember(t,keep_inds))
        htext(t).String = '';
    end
    htext(t).FontSize = 24;
    if length(keep_inds)>2
        if t==keep_inds(3) % if we're dealing with the 180 label
            htext(t).Position(1) = pos180;
        end
    end
end

% Fix line widths
hlines = findall(hf,'Type','line');
for i = 1:length(hlines)
set(hlines(i),'LineWidth',1.5);
end

% Draw outline darker
ymax = [];
for i = 1:length(hlines)
    if size(hlines(i).YData,2)>100
        ymax = [ymax; [i, max(hlines(i).YData)]];
    end
end
[~, outmax] = max(ymax(:,2));
outline = ymax(outmax,1);

% Fix line colors
set(hlines(3:end),'Color',axGray.* [1 1 1]);
set(hlines(outline),'Color',outerGray.*[1 1 1]);
% set(hf,'color','none');

% Reposition and scale figure
hf.Position = [600 300 400 400];

end

