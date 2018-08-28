function [hf, Fim] = plotS2D(ims,figSave)

% Hard coded figure settings that look nice
if ispc
    font=16;
    flfont=14;
    flpos=[0.7, 0.85];
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

% Initialize figure
hf=figure;

ax = gca;
hold(ax,'on');

frames = ims.op2d_fixed.xdata ./ 1000;

area_rand = area(ax,frames, ims.op2d_fixed.S_rand);
area_rand.FaceColor = [1 0 0];
area_rand.FaceAlpha = 0.3;
area_im = area(ax,frames, ims.op2d_fixed.S_im);
area_im.FaceColor = [0 0 1];
area_im.FaceAlpha = 0.3;

% Rescale data and fit from um to nm and plot
frames = ims.op2d_fixed.xdata ./ 1000;
fineFrames = linspace(frames(1),frames(end),500);
BETA=[ims.op2d_fixed.Sfull; ims.op2d_fixed.decayLen / 1000];
fitY = BETA(1)+(1-BETA(1)).*exp(-fineFrames./(2*BETA(2)));  % Use fineFrames to make figure look good

p1 = plot(ax,frames,ims.op2d_fixed.S_im,'ok'); %, ...
p2 = plot(ax,fineFrames,fitY,'-b');

% Make things look nice
xlabel('Frame Size (µm)');
ylabel('{\itS}_{2D}')
ax.FontSize = font;
ax.XLim = [0 ceil(frames(end))];
ax.Box = 'on';
ax.LineWidth = 0.75;
ax.PlotBoxAspectRatio = [1 1 1];
ax.YLim = [-0.5 0.5];
ax.TickLength = [0.015, 0.015];
p1.LineWidth = markerline;
p1.MarkerSize = marker;
p2.LineWidth = line;
p2.Color = [0 0 1];
% p(2).LineWidth = line; p(2).Color = [1, 0, 0];

htex = text('Units', 'normalized', 'Position', flpos, ...
    'BackgroundColor', [1 1 1], ...
    'String', {['{\it\lambda}_{C} = ', num2str(BETA(2)*1000,4), ' nm'],...
               ['{\itS}_{full} = ', num2str(BETA(1),2)],...
               ['{\ita} = ', num2str(ims.op2d_fixed.a,2)]},...
    'FontSize', flfont,...
    'EdgeColor', edgedark*[1 1 1],...
    'LineWidth', edgewidth);

hf.Position=position;

F = getframe(hf);
Fim = F.cdata;

if figSave
    fig_file = [ims.figSavePath, '_OP2D', '.tif'];
    ensure_dir(fig_file);
    hgexport(hf, fig_file,  ...
        hgexport('factorystyle'), 'Format', 'tiff');
    close(hf)
end

end
