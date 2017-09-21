function [handles, Fim] = FiberVecPlot_stitch(handles,figSave)

% Given imageData filepath, plot the fibers.

if exist('figSave','var')~=1
    figSave=0;
end

XY = {handles.ims.Fibers(:).xy};
w = size(handles.ims.img,2);
h = size(handles.ims.img,1);

f1 = figure('Visible','on');
if h>512
    f1.Position = [1 1 w/h*512 512];
else
    f1.Position = [1 1 w h];
end
f1.Color = [1 1 1];
ha = axes('parent',f1);
hold on

for i = 1:length(XY)
    color = [0 0 1]; %rand(1,3);
    XYi = XY{i};
    plot(ha,XYi(1,:),XYi(2,:),'Color',color,'LineWidth',1)
end
% axis equal
set(ha,'Ydir','reverse')
ax = ha;
ax.XLim = [0 w];
ax.YLim = [0 h];
ax.Position = [0 0 1 1];

ax.Visible = 'off';
ax.PlotBoxAspectRatio = [1 1 1];
%
% F = getframe(f1);
% Fim = F.cdata;
% Fres = imresize(Fim,[h, w]);

handles = plotGT(ax,handles,'fiber_axes');

F = getframe(f1);
Fim = F.cdata;

if figSave
    Fres = imresize(Fim,[h, w]);
    imwrite(Fres, [handles.ims.figSavePath, '_FiberPlot.tif']);
end

close(f1)

end
