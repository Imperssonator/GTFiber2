function [hScat,ax] = GoodBox(GroupSort,MobSort)

[xLabs, groupSplits, xVals] = unique(GroupSort,'stable');

% Find the median of each group and re-sort the group labels based on this

Medians=zeros(length(xLabs),1);
for g = 1:length(xLabs)
    SetV = MobSort(xVals==g,1);
    Medians(g) = median(SetV);
end
ReSort = [(1:length(xLabs))',Medians];
ReSort=sortrows(ReSort,-2);

xLabs_old=xLabs;
xVals_old=xVals;
for i = 1:length(xLabs)
    xLabs{i} = xLabs_old{ReSort(i,1)};
    xVals(xVals_old==i) = find(ReSort(:,1)==i);
end

xValsJog = xVals+(rand(length(xVals),1)-0.5)/5;
if ~iscell(xLabs)
    xLabs = mat2cell(xLabs,size(xLabs,1),size(xLabs,2));
end

% Create figure
f1=figure;
hold on

% Create scatterplot
hScat = scatter(xValsJog,MobSort); drawnow;
ax = gca;

% Separate Groups and plot medians / 25th and 75th percentiles as lines
for g = 1:max(xVals)
    SetV = MobSort(xVals==g,1);
    SetX = g;
    med = median(SetV);
    p25 = prctile(SetV,25);
    p75 = prctile(SetV,75);
    boxWid = 0.25;
    plot(ax,[SetX-boxWid, SetX+boxWid], [med, med], '-k')
    plot(ax,[SetX-boxWid, SetX+boxWid], [p25, p25], '-k')
    plot(ax,[SetX-boxWid, SetX+boxWid], [p75, p75], '-k')
    plot(ax,[SetX-boxWid, SetX-boxWid], [p25, p75], '-k')
    plot(ax,[SetX+boxWid, SetX+boxWid], [p25, p75], '-k')
end

ax.XLim = [0 max(xVals)+1];
ax.XTickMode = 'manual';
ax.XTick = (1:max(xVals));
ax.XTickLabels = xLabs;
ax.FontSize=20;
% ax.YScale = 'log';
ax.Box = 'on';
ax.XTickLabelRotation = 40;

% Make Transparent Labels
hScat.SizeData = 100;
hScat.MarkerEdgeColor = uint8([70, 70, 70]);
hScat.MarkerEdgeAlpha = 200/255;
hScat.MarkerFaceColor = uint8([120, 120, 120]);
hScat.MarkerFaceAlpha = 150/255;

hf=gcf;
hf.Position = [440 479 512 319];

end