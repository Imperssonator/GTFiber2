function out = RemoveHoles(IM)

tic
L = bwlabel(IM);
RP = regionprops(IM,'EulerNumber');
out = zeros(size(L));

for i = 1:length(RP)
    if RP(i).EulerNumber > 0
        out = out+double(L==i);
    end
end

out = out==1;
toc

end