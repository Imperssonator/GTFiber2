load('s2t48b2');
ims_son=ims;
load('u8t48b2');
ims_uv=ims;

uv_fit = fitdist(ims_uv.FWD','normal');
son_fit = fitdist(ims_son.FWD','normal');
yu = pdf(uv_fit,X).*length(ims_uv.FWD);
ys = pdf(son_fit,X).*length(ims_son.FWD);

navy = [31; 119; 180]./255;
green = [44; 160; 44]./255;
red = [214; 39; 40]./255;

figure;ax=gca;hold on
hu=histogram(ax,ims_uv.FWD,10:1:40);
hs=histogram(ax,ims_son.FWD,10:1:40);

X=10:0.1:40;
pu=plot(ax,X,yu,'-b');
ps=plot(ax,X,ys,'-r');

hu.FaceColor=navy;
hs.FaceColor=red;
pu.Color=navy;
ps.Color=red;
pu.LineWidth=1.5;
ps.LineWidth=1.5;

ax.FontSize=16;
xlabel('Fiber Width (nm)')
ylabel('Number of Fibers')