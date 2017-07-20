function [] = LogNormHist(dist,bins)

navy = [31; 119; 180]./255;
figure;
histfit(dist,bins,'lognormal')
ax=gca;
ax.FontSize=22;
ax.XLim=[0 3000];
xlabel('Fiber Length (nm)')
ylabel('Number of Fibers')
ax.Children(1).Color=[0 0 0];
ax.Children(2).FaceColor=navy;
ax.Children(2).FaceAlpha=0.5;

end