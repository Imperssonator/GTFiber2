function FiberIm = ContourImage(ims,segNum,LineWidth)

% Given ims and a segment number, make a binary image of the pixels that
% its "seeking" contour hits.

if exist('LineWidth','var')~=1
    LineWidth = 3;
end

w = size(ims.img,2);
h = size(ims.img,1);

f1 = figure('Visible','off');
f1.Position = [1 1 1+w 1+h];
f1.Color = [1 1 1];
ha = axes('parent',f1);
hold on

for s = 1:length(segNum)
    XYi = ims.fibSegs(segNum(s)).xy_seek;
    plot(ha,XYi(1,:),XYi(2,:),'-k','LineWidth',LineWidth)
end

% axis equal
set(ha,'Ydir','reverse')
ax = ha;
ax.XLim = [0 w];
ax.YLim = [0 h];

ax.Visible = 'off';
% ax.PlotBoxAspectRatio = [1 1 1];
ax.Position = [0 0 1 1];

F = getframe(f1);
Fim = F.cdata;
FiberIm = ~im2bw(imresize(Fim,[h, w]),0.3);
close(f1)
% imtool(FiberIm)

end