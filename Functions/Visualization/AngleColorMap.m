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
