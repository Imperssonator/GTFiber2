function [him, hp, ho, hax, hf] = contour_plot(img,xy_sel)

hf=figure;
him=imshow(img);
hax=gca;
hold on

for i = 1:length(xy_sel)
    xy=xy_sel{i};
    hp = plot(xy(1,:),xy(2,:),'-b','LineWidth',2);
    ho = plot(xy(1,:),xy(2,:),'ob','MarkerSize',12,'LineWidth',2);
end

hax.Position = [0 0 1 1];
him.AlphaData=0.35;

end