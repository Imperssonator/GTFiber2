function AngIm = AngleColorMap(AngMap, BW)

sat = 0.88;

AngMapNaN = AngMap;
AngMapNaN(~BW) = NaN;
AngMapNaN(AngMapNaN<0) = AngMapNaN(AngMapNaN<0)+180;
[m, n] = size(AngMap);
AngIm = zeros(m,n,3);
for i = 1:m
    for j = 1:n
        if isnan(AngMapNaN(i,j))
            AngIm(i,j,:) = [0 0 0];
        else
            AngIm(i,j,:) = hsv2rgb([2*AngMapNaN(i,j)/360, 1, sat]);
        end
    end
end

AngImDiskFilt = imfilter(AngIm,fspecial('disk',2),'replicate');
AngImHSV = rgb2hsv(AngImDiskFilt);
AngImHSV(:,:,3) = ones(size(AngImHSV,1),size(AngImHSV,2)).*AngImHSV(:,:,3)>0.1;

amf = figure; imshow(imgaussian(hsv2rgb(AngImHSV),0.7));
ampos = amf.Position;
legpos = [ampos(1)+ampos(3), ampos(2)];
Angle_Legend(legpos);

ODist(AngMapNaN,5,'');

end

function [] = Angle_Legend(pos)

h = 370;
wl = -425; wr = 370;

[X Y] = meshgrid((wl:wr),(h:-1:0));
R = sqrt((X.^2+Y.^2));
H = atan2d(Y,X)*2/360;
S = ones(size(X));
V = ones(size(X));

rr = 220;
Black = R>rr;
Color = R<=rr;

V(Black) = 0;
V(Color) = 0.84;

HSV = cat(3,H,S,V);

rgb = hsv2rgb(HSV);

rgbt1 = insertText(rgb,[0 250],['180', char(176)],'BoxColor','black','TextColor','white','FontSize',80,'BoxOpacity',0);
rgbt2 = insertText(rgbt1,[640 250],['0', char(176)],'BoxColor','black','TextColor','white','FontSize',80,'BoxOpacity',0);
rgbt3 = insertText(rgbt2,[350 30],['90', char(176)],'BoxColor','black','TextColor','white','FontSize',80,'BoxOpacity',0);

legfig = figure; imshow(rgbt3)
legfig.Position(1:2) = pos;


end