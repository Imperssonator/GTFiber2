function [hf, Fim] = plotOrCorr(ims,figSave)

% Hard coded figure settings that look nice
if ispc
    font=16;
    flfont=14;
    flpos=[0.7, 0.9];
    position = [227  413  599  569];
else
    font=20;
    flfont=16;
    flpos=[0.6, 0.85];
    position = [390   143   500   450];
end
marker = 7;
markerline = 1;
line = 1.25;
lenscale = 1000;
edgedark = 0;
edgewidth = 0.75;

hh = ims.OrCorr.hist;
numBins = size(hh,2);
maxDist = hh(3,end);

maxPlotR = ims.nmWid/2;
maxPlotBin = ceil(maxPlotR*numBins/maxDist);

[~, minPlotBin] = max(hh(2,:),[],2);
R = hh(3,minPlotBin:maxPlotBin);
Y = hh(2,minPlotBin:maxPlotBin);

% Calculate correlation length
corr_len = ims.OrCorr.corrLen;
disp(corr_len)

% Build plot
hf=figure;
ax=gca; hold(ax,'on')

hp = plot(ax,R,Y,'ok');
hfit = plot(ax,R,exp(-R./corr_len),'-b');
hline = plot(ax,[corr_len, corr_len],[0,1],'-r');

hp.MarkerSize = marker;
hp.LineWidth = markerline;
hfit.LineWidth = line;
hLine.LineWidth = line;

xlabel('Separation, r (nm)')
ylabel('<cos(2\theta)>')
ax.FontSize = font;
ax.XLim = [0 ceil(maxPlotR)];
ax.Box = 'on';
ax.LineWidth = 0.75;
ax.PlotBoxAspectRatio = [1 1 1];
ax.TickLength = [0.015, 0.015];

htex = text('Units', 'normalized', 'Position', flpos, ...
    'BackgroundColor', [1 1 1], ...
    'String', {['{\it\kappa} = ', num2str(corr_len,4), ' nm']}, ...
    'FontSize', flfont,...
    'EdgeColor', edgedark*[1 1 1],...
    'LineWidth', edgewidth);

hf.Position=position;

F = getframe(hf);
Fim = F.cdata;

if figSave
    hgexport(hf, [ims.figSavePath, '_OrCorr', '.tif'],  ...
        hgexport('factorystyle'), 'Format', 'tiff');
    close(hf)
end

end