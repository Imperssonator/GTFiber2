function [] = FiberColors(ims,segorfib,alpha,figSave)

% Given imageData filepath, plot the fibers.

if exist('figSave','var')~=1
    figSave=0;
end

figure;
ax = gca;
hold on

switch segorfib
    case 'seg'
        
        for i = 1:length(ims.fibSegs)
            color = rand(1,3);
            plot(ax,ims.fibSegs(i).xy(1,:),ims.fibSegs(i).xy(2,:),'Color',color,'LineWidth',1.5);
        end
        
    case 'fib'
        for i = 1:length(ims.Fibers)
            color = rand(1,3);
            plot(ax,ims.Fibers(i).xy(1,:),ims.Fibers(i).xy(2,:),'Color',color,'LineWidth',1.5);
        end
        
end

hf=gcf;
[h, w] = size(ims.gray);
hf.Position = [1 1 1+w 1+h];
ax.XLim = [0 w];
ax.YLim = [0 h];
ax.Position = [0 0 1 1];

if figSave
    F = getframe(hf);
    Fim = F.cdata;
    Fres = imresize(Fim,[h, w]);
    imwrite(Fres, [ims.figSavePath, '_FibVec_alpha', num2str(alpha), '.tif']);
    close(hf)
end


end
