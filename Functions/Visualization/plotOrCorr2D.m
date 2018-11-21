function [hf, Fim] = plotOrCorr2D(ims,figSave)

% create Light Bartlein orange-white-purple diverging scheme
LB=flipud(lbmap(256,'BrownBlue'));

hf = figure;
ax = axes('Parent',hf);
hold(ax,'on')

H = ims.OrCorr2D.array;
X = H(:,:,3)./1000;
Y = H(:,:,4)./1000;
C = H(:,:,2);

colormap(LB);
pc = pcolor(ax,X,Y,C);
cb = colorbar();
caxis([-1, 1]);
set(cb,'Limits',[-1,1]);

xlabel('dx (µm)');
ylabel('dy (µm)');
ylabel(cb, '<cos(2\theta)>')
set(ax,'FontSize',16);
set(pc,'LineStyle','none');


xbounds = [X(1,round(0.25*size(X,2))), X(1,round(0.75*size(X,2)))];
ybounds = [Y(round(0.75*size(Y,1)),1), Y(round(0.25*size(Y,1)),1)];
xlim(ax,xbounds);
ylim(ax,ybounds);

hf.Position = [560 569 516 379];

F = getframe(hf);
Fim = F.cdata;

if figSave
    fig_file = [ims.figSavePath, '_OC2D', '.tif'];
    ensure_dir(fig_file);
    hgexport(hf, fig_file,  ...
        hgexport('factorystyle'), 'Format', 'tiff');
    close(hf)
end

end
