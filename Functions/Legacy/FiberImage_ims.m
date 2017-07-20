function [] = FiberImage_ims(IMS)

FL = IMS.FiberLabels;
[m n] = size(FL);
PlotIm = zeros(m,n,3);          % It's gonna be an image

Fibers = IMS.Fibers;
f = length(Fibers);

BGColor = [255 255 255];
BG = double(FL==0);
BGR = BG.*BGColor(1); PlotIm(:,:,1) = BGR;
BGG = BG.*BGColor(2); PlotIm(:,:,2) = BGG;
BGB = BG.*BGColor(3); PlotIm(:,:,3) = BGB;

ColorList = ...
    [217 83 25;...
     237 177 32;...
     126 47 142;...
     119 172 48;...
     77 190 238;...
     162 20 47;...
     0 114 189];

nColors = size(ColorList,1);
count = 1;

for i = 1:f
    FiberPix = double(FL==i);
    C = randi(255,1,3); %ColorList(count,:);
    FibR = FiberPix.*C(1);
    FibG = FiberPix.*C(2);
    FibB = FiberPix.*C(3);
    PlotIm(:,:,1) = PlotIm(:,:,1) + FibR;
    PlotIm(:,:,2) = PlotIm(:,:,2) + FibG;
    PlotIm(:,:,3) = PlotIm(:,:,3) + FibB;
%     if count == nColors
%         count = 1;
%     else
%         count = count+1;
%     end
end 

PlotIm = uint8(PlotIm);

figure
imshow(PlotIm)

end