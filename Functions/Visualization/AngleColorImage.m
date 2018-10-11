function AngIm = AngleColorImage(AngMap, BW)

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

end