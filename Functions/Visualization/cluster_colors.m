function [] = cluster_colors(ims,VV)

% Given imageData filepath, plot the fibers.

figure;
him = imshow(ims.gray);
him.AlphaData = 0.6;
ax = gca;
hold on

for i = 1:max(VV)
    color = rand(1,3);
    segs = find(VV==i);
    for s = 1:length(segs)
        plot(ax,ims.fibSegs(segs(s)).xy(1,:),ims.fibSegs(segs(s)).xy(2,:),'Color',color,'LineWidth',2);
    end
end

ax.Position = [0 0 1 1];


end
