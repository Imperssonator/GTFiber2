function [] = FiberAlpha(ims,segorfib,alpha,figSave)

% Given imageData filepath, plot the fibers.

if exist('figSave','var')~=1
    figSave=0;
end

figure;
him = imshow(ims.gray);
him.AlphaData = alpha;
ax = gca;
hold on

switch segorfib
    case 'seg'
        
        for i = 1:length(ims.fibSegs)
            color = [0 0 1]; %rand(1,3);
            plot(ax,ims.fibSegs(i).xy(1,:),ims.fibSegs(i).xy(2,:),'Color',color,'LineWidth',1.5);
        end
        
    case 'fib'
        for i = 1:length(ims.Fibers)
            color = [0 0 1]; %rand(1,3);
            plot(ax,ims.Fibers(i).xy(1,:),ims.Fibers(i).xy(2,:),'Color',color,'LineWidth',1.5);
        end
        
end

ax.Position = [0 0 1 1];
hf=gcf;
[h, w] = size(ims.gray);
hf.Position = [1 1 1+w 1+h];

if figSave
    F = getframe(hf);
    Fim = F.cdata;
    Fres = imresize(Fim,[h, w]);
    imwrite(Fres, [ims.figSavePath, '_FibVec_alpha', num2str(alpha), '.tif']);
    close(hf)
end


end
