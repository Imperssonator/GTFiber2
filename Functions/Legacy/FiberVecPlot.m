function [] = FiberVecPlot(ims,gt_fa)

% Given imageData filepath, plot the fibers.

switch gt_fa
    case 'gt'
        
        XY = {ims.Fibers(:).xy};
        w = size(ims.img,2);
        h = size(ims.img,1);
        
        f1 = figure;
        f1.Position = [1 1 1+w 1+h];
        f1.Color = [1 1 1];
        ha = axes('parent',f1);
        hold on
        
        for i = 1:length(XY)
            XYi = XY{i};
            plot(ha,XYi(1,:),XYi(2,:),'-b','LineWidth',1.5)
        end
        axis equal
        set(ha,'Ydir','reverse')
        ax = ha;
        ax.XLim = [0 w];
        ax.YLim = [0 h];
        ax.Visible = 'off';
        % ax.PlotBoxAspectRatio = [1 1 1];
        ax.Position = [0 0 1 1];
        
    case 'fa'
        load(ims)
        f1 = figure;
        hold on
        
        XY = imageData.xy;
        w = imageData.sizeX;
        h = imageData.sizeY;
        
        for i = 1:length(XY)
            XYi = XY{i};
            plot(XYi(1,:),XYi(2,:),'-b','LineWidth',1.5)
        end
        axis equal
        set(gca,'Ydir','reverse')
        ax = gca;
        ax.Visible = 'off';
        ax.XLim = [0 w];
        ax.YLim = [0 h];
        ax.PlotBoxAspectRatio = [1 1 1];
        f1.Position = [1 1 1+w 1+h];
        f1.Color = [1 1 1];
        ax.Position = [0 0 1 1];
end

end
