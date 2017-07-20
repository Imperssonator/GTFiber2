function out = RemoveNonLines(IM)

SegLabels = bwlabel(IM,8);                % Create an image where their labels correspond to the order regionprops found them in
NumSegs = max(SegLabels(:));

% Generate Endpoints
EndsImg = bwmorph(IM,'endpoints');        % Also create a binary matrix of the segment endpoints

% Remove Perfect Rings
out = SegLabels;
for s = 1:NumSegs;
    SegEnds = (SegLabels==s) & EndsImg;
    if sum(SegEnds(:))<2
        out(SegLabels==s)=0;
    end
end

out=~~out;

end