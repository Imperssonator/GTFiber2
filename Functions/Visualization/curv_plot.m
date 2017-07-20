function [him, hs, hax, hf, hc] = curv_plot(img,xy,curv,lims)

if exist('lims')~=1
    lims = [min(curv),max(curv)];
end

buff = 10;
zoom = 5;
MarkerSize=30;
bbox = [min(xy(1,:))-buff,...
        max(xy(1,:))+buff,...
        min(xy(2,:))-buff,...
        max(xy(2,:))+buff];


m=64;
cmap=colormap([gray(256);parula(m)]);

figure;
him=imshow(imadjust(img),cmap);

% Plot the points of the fibers, but also plot two points that are the min
% and max curv values so the colormap doesn't try to rescale and screw things up
hold on
hs=scatter([1 1],[1 1],1,[257,257+m-1]);
hs(2)=scatter(xy(1,:),xy(2,:),MarkerSize,'LineWidth',1);
hold off

hs(2).MarkerFaceColor = 'flat';
him.AlphaData = 0.4;

% This is kind of weird, but you basically want to map from [cmin cmax] to
% [1 64], and cmin is usually 0.
cmin = lims(1); cmax = lims(2);
curv_color = min(m,round((m-1)*(curv-cmin)/(cmax-cmin))+1);
hs(2).CData=curv_color+256;

hc=colorbar();
hc.Limits=[257, 257+m-1];
hc.Ticks = linspace(257,256+m,5);
hc.TickLabels = ...
    cellfun(@(x) num2str(x*1000,2),...  % Multiply by 1000 to scale to 1/um
            mat2cell(linspace(cmin,cmax,5)',...
                     ones(size(linspace(cmin,cmax,5)',1),1),...
                     1),...
            'UniformOutput',false);
hc.FontSize=16;

hf=gcf;
hax=gca;


% hax.XLim = bbox(1:2);
% hax.YLim = bbox(3:4);
% hax.Position = [0 0 1 1];
% hf.Position(3:4) = [(bbox(2)-bbox(1))*zoom,...
%                     (bbox(4)-bbox(3))*zoom];

end